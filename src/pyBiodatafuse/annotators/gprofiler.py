# coding: utf-8

"""Python file for adding enrichment analysis using g:Profiler Python package. All the pathways and annotations are being added despite not being significance."""

import datetime
from typing import Any, Dict

import pandas as pd
import requests
from gprofiler import GProfiler

import pyBiodatafuse.constants as Cons
from pyBiodatafuse.utils import get_identifier_of_interest


def get_data_versions(species: str = "hsapiens") -> dict:
    """Get version of g:Profiler.

    :param species: The species to retrieve the version information for.
    :returns: a dictionary containing the version information
    """
    params = {Cons.ORGANISM: species}
    try:
        response = requests.get(Cons.GPROFILER_VERSION_ENDPOINT, params=params)

        if response.status_code == 200:
            gprofiler_version = response.json()
            return gprofiler_version
        else:
            raise Exception(f"Failed to retrieve data: {response.status_code}")
    except Exception as e:
        return {"error": str(e)}


def gene_enrichment_gprofiler(
    bridgedb_df: pd.DataFrame,
    padj_colname: str = "padj",
    padj_filter: float = 0.05,
    species: str = "hsapiens",
) -> tuple[pd.DataFrame, dict]:
    """Enrichment analysis using g:Profiler for input gene list.

    :param bridgedb_df: DataFrame containing gene data with columns "identifier" and a significance column.
    :param padj_colname: Name of the column used to filter significant genes (default is "padj").
    :param padj_filter: Significance threshold to filter genes (default is 0.05).
    :param species: Species code for g:Profiler query, e.g., "hsapiens" for human (default is "hsapiens").
    :returns: A tuple containing:
              - DataFrame containing g:Profiler analysis results.
              - Dictionary with g:Profiler version information.
        :raises RuntimeError: If the g:Profiler query fails.
    :raises RuntimeError: If the g:Profiler query fails.
    """
    # Extract the "target" values in bridgedb_df
    data_df = get_identifier_of_interest(bridgedb_df, Cons.GPROFILER_GENE_INPUT_ID)
    query_genes = (
        data_df[data_df[f"{padj_colname}_dea"] <= padj_filter][Cons.TARGET_COL].unique().tolist()
    )
    background_genes = data_df[Cons.TARGET_COL].unique().tolist()

    # Record the start time
    start_time = datetime.datetime.now()

    gp = GProfiler(return_dataframe=True)
    try:
        gprofiler_df = gp.profile(
            organism=species,
            all_results=True,
            query=query_genes,
            background=background_genes,
            no_evidences=False,
            significance_threshold_method="fdr",
            user_threshold=0.05,
        )
    except Exception as e:
        raise RuntimeError(f"g:Profiler query failed: {e}")

    # Record the end time
    end_time = datetime.datetime.now()

    """Metadata details"""
    # Get the current date and time
    current_date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    # Calculate the time elapsed
    time_elapsed = str(end_time - start_time)

    # Add version, datasource, query, query time, and the date to metadata
    # Retrieve g:Profiler version information for the specified species
    gprofiler_version = get_data_versions(species)
    gprofiler_metadata: Dict[str, Any] = {
        "datasource": Cons.GPROFILER,
        "metadata": gprofiler_version,
        "query": {
            "size": len(query_genes),
            "time": time_elapsed,
            "date": current_date,
        },
    }

    gprofiler_df.rename(columns={"native": Cons.GPROFILER_ID}, inplace=True)
    gprofiler_df[Cons.DATASOURCE] = Cons.GPROFILER

    gprofiler_df = gprofiler_df[
        ~gprofiler_df["parents"].apply(lambda x: x == [])
    ]  # rm the root terms

    return gprofiler_df, gprofiler_metadata


def _create_path_info(row, gprofiler_df: pd.DataFrame):
    """Inner helper function to create a dictionary of path info for each row.

    :param row: Row of the DataFrame.
    :param gprofiler_df: DataFrame containing raw g:Profiler analysis results.
    :returns: Dictionary of path info for each row.
    """
    path_info = {
        col: row[col]
        for col in gprofiler_df.columns
        if col not in [Cons.GPROFILER_INTERSECTIONS, Cons.SOURCE_COL]
    }
    return path_info


def process_gprofiler_data(gprofiler_df: pd.DataFrame) -> pd.DataFrame:
    """Process raw g:Profiler results.

    :param gprofiler_df: DataFrame containing raw g:Profiler analysis results.
    :returns: Processed DataFrame structured by intersections and sources.
    """
    # Create a 'gprofiler' column using the helper function
    gprofiler_df[Cons.GPROFILER_RESULT_COL] = gprofiler_df.apply(
        _create_path_info, axis=1, args=(gprofiler_df,)
    )

    # Drop all columns except for "source", "id", "intersections", and "gprofiler"
    cols_to_drop = [col for col in gprofiler_df.columns if col not in Cons.GPROFILER_COLS_TO_KEEP]
    gprofiler_df = gprofiler_df.drop(columns=cols_to_drop)

    # Explode the 'intersections' column so that each intersection gets its own row
    gprofiler_df = gprofiler_df.explode(Cons.GPROFILER_INTERSECTIONS).reset_index(drop=True)

    # Prepare a final DataFrame with unique intersections
    unique_sources = sorted(gprofiler_df[Cons.SOURCE_COL].unique())
    intermediate_df = pd.DataFrame()
    intermediate_df[Cons.GPROFILER_INTERSECTIONS] = gprofiler_df[
        Cons.GPROFILER_INTERSECTIONS
    ].unique()

    # For each unique source, group data by intersections and map to the final DataFrame
    for source in unique_sources:
        source_subset = gprofiler_df[gprofiler_df[Cons.SOURCE_COL] == source]
        source_dictionaries = (
            source_subset.groupby(Cons.GPROFILER_INTERSECTIONS)[Cons.GPROFILER_RESULT_COL]
            .apply(list)
            .to_dict()
        )
        column_name = f"{Cons.GPROFILER}_{source.lower()}"
        intermediate_df[column_name] = intermediate_df[Cons.GPROFILER_INTERSECTIONS].map(
            source_dictionaries
        )

    return intermediate_df


def get_gene_enrichment(
    bridgedb_df: pd.DataFrame,
    species: str = "hsapiens",
    padj_colname: str = "padj",
    padj_filter: float = 0.05,
) -> tuple[pd.DataFrame, dict]:
    """Enrichment analysis using g:Profiler and retrieve version info.

    :param bridgedb_df: DataFrame containing gene data with columns "identifier"
                        and a significance column.
    :param species: species for both version retrieval and g:Profiler query
                     (default is "hsapiens").
    :param padj_colname: Name of the column used to filter significant genes
                       (default is "padj").
    :param padj_filter: Significance threshold to filter genes (default is 0.05).
    :returns: A tuple containing:
              - Processed DataFrame from g:Profiler analysis.
              - Dictionary with g:Profiler version information.
    """
    # Perform gene enrichment analysis to get raw g:Profiler results
    gprofiler_df, gprofiler_metadata = gene_enrichment_gprofiler(
        bridgedb_df=bridgedb_df, padj_colname=padj_colname, padj_filter=padj_filter, species=species
    )

    # Process the raw g:Profiler results
    intermediate_df = process_gprofiler_data(gprofiler_df)
    intermediate_df = intermediate_df.rename(columns={"intersections": Cons.TARGET_COL})

    # Merge the processed DataFrame with the original bridgedb_df
    data_df = bridgedb_df[bridgedb_df[Cons.TARGET_SOURCE_COL] == Cons.GPROFILER_GENE_INPUT_ID]

    merged_df = (
        pd.merge(  # TODO: check if we can modify collapse_data_sources to handle multiple columns
            data_df,
            intermediate_df,
            on=Cons.TARGET_COL,
            how="left",
        )
    )

    return merged_df, gprofiler_metadata
