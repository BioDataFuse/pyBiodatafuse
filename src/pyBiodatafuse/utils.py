# coding: utf-8

"""Python utils file for global functions."""

import warnings
from typing import List

import pandas as pd
from collections import defaultdict

from pyBiodatafuse.id_mapper import read_resource_files
from pyBiodatafuse.annotators import bgee, disgenet, minerva, molmedb, opentargets, pubchem, stringdb, wikipathways


def get_identifier_of_interest(bridgedb_df: pd.DataFrame, db_source: str) -> pd.DataFrame:
    """Get identifier of interest from BridgeDb output file.

    :param bridgedb_df: DataFrame containing the output from BridgeDb
    :param db_source: identifier of interest from BridgeDB (e.g. "NCBI Gene")
    :returns: a DataFrame containing the identifiers of interest
    """
    # Load identifier options
    identifier_options = read_resource_files()["source"].tolist()

    # Check if source is in identifier options
    assert db_source in identifier_options, f"Source {db_source} is not in identifier options"

    # Filter rows where "target.source" is specific datasource for eg. "NCBI Gene"
    return bridgedb_df[bridgedb_df["target.source"] == db_source]


def create_or_append_to_metadata(data: dict, prev_entry: List[dict]) -> List[dict]:
    """Create and/or append data to a metadata file.

    :param data: dictionary of data to be saved to the metadata file.
    :param prev_entry: list of dictionaries containing the previous data
        The metatdata file has the following schema:
        {
            "datasource": name_of_datasource,
            "metadata": {
                "source_version": {source_version_info},
                "data_version": {data_version_info} (Optional)
            },
            "query": {
            "size": number_of_results_queried,
            "time": time_taken_to_run_the_query,  (using datetime.datetime.now())
            "date": date_of_query,
            "url": url_of_query,
            "request_string": post_request_string (Optional)

            }
        }
    :returns: a metadata dictionary
    """
    # Create a metadata file if it doesn't exist
    prev_sources = [data["datasource"] for data in prev_entry if "datasource" in data.keys()]

    assert isinstance(data, dict), "Unsupported data type. Only dict is supported."

    if data["datasource"] not in prev_sources:
        prev_entry.append(data)

    return prev_entry


def collapse_data_sources(
    data_df: pd.DataFrame,
    source_namespace: str,
    target_df: pd.DataFrame,
    common_cols: list,
    target_specific_cols: list,
    col_name: str,
) -> pd.DataFrame:
    """Collapse data sources into a single column.

    :param data_df: BridegDb dataFrame containing idenfitiers from all sources
    :param source_namespace: identifier of interest from BridgeDB (e.g. "NCBI Gene")
    :param target_df: DataFrame containing data from a external source
    :param common_cols: list of columns that are common to both dataframes and can be used to merge
    :param target_specific_cols: list of columns that are specific to the external source
    :param col_name: name of the new column to be created
    :returns: a DataFrame containing the new data columns for a new resource
    """
    data_df = data_df[data_df["target.source"] == source_namespace]

    if target_df.empty:
        # If the target_df is empty, then return the data_df as is
        data_df[col_name] = None
        data_df.reset_index(inplace=True, drop=True)
        return data_df

    merged_df = pd.merge(data_df, target_df, on=common_cols, how="left")

    # Create a new source column with values from selected columns as a list
    merged_df[col_name] = merged_df[target_specific_cols].apply(lambda row: row.to_dict(), axis=1)
    # Convert source column from string to a list of strings
    merged_df[col_name] = merged_df[col_name].apply(lambda x: [x])

    # Group by the first 4 columns and aggregate the values into a list
    cols_of_interest = data_df.columns.tolist()
    merged_df = merged_df.groupby(cols_of_interest)[col_name].sum().reset_index()

    return merged_df


def combine_sources(df_list: List[pd.DataFrame]) -> pd.DataFrame:
    """Combine multiple dataframes into a single dataframe.

    :param df_list: list of dataframes to be combined
    :returns: a single dataframe containing from a list of dataframes
    """
    m = pd.concat(df_list, axis=1)
    m = m.loc[:, ~m.columns.duplicated()]  # remove duplicate columns

    return m


def check_columns_against_constants(
    data_df: pd.DataFrame, output_dict: dict, check_values_in: list
):
    """Check if columns in the data source output DataFrame match expected types and values from a dictionary of constants.

    :param data_df: DataFrame to check.
    :param output_dict: Dictionary containing expected types for columns.
    :param check_values_in: List of column names to check values against constants.
    """
    for col, expected_type in output_dict.items():
        if col not in data_df.columns:
            warnings.warn(f"Column '{col}' is missing in the DataFrame.", stacklevel=2)
            continue

        if not data_df[col].dropna().apply(type).eq(expected_type).all():
            warnings.warn(
                f"Not all values in column '{col}' have the correct type '{expected_type}'.",
                stacklevel=2,
            )
        if col in check_values_in:
            exec(f"from pyBiodatafuse.constants import {col.upper()}")  # noqa: S102
            starts_with = locals()[col.upper()]
            if not data_df[col].apply(type).eq(int).all():
                prefixes = starts_with.split("|")
                if (
                    not data_df[col]
                    .dropna()
                    .apply(
                        lambda value, prefixes=prefixes: any(
                            value.startswith(prefix) for prefix in prefixes
                        )
                    )
                    .all()
                ):
                    warnings.warn(
                        f"All values in column '{col}' do not start with '{starts_with}'.",
                        stacklevel=2,
                    )


def process_selected_sources(
    bridgedb_df: pd.DataFrame, selected_sources_list: list
) -> pd.DataFrame:
    """query the selected databases and convert the output to a dataframe.

    @param bridgedb_df: BridgeDb output for creating the list of gene ids to query
    @param selected_sources_list: list of selected databases
    """
    # Initialize variables
    combined_data = pd.DataFrame()
    combined_metadata = defaultdict(lambda: defaultdict(str))
    # Dictionary to map the datasource names to their corresponding functions
    data_source_functions = {
        ## TODO: "bgee": bgee.get_gene_expression,
        "disgenet": disgenet.get_gene_disease,
        ## TODO: "minerva": minera.get,
        ## TODO: "molmedb": molmedb.get_mol_gene_inhibitor,
        ## TODO: "pubchem"
        "opentarget.gene_ontology": opentargets.get_gene_go_process,
        "opentarget.reactome": opentargets.get_gene_reactome_pathways,
        "opentarget.drug_interactions": opentargets.get_gene_compound_interactions,
        "opentarget.disease_associations": opentargets.get_gene_disease_interactions,
        "string": stringdb.get_ppi,
        "wikipathways": wikipathways.get_gene_wikipathway
        ## TODO: "wikidata"
    }
    warnings = []  # Initialize empty list for warnings

    for source in selected_sources_list:
        if source in data_source_functions:
            print(source)
            tmp_data, tmp_metadata = data_source_functions[source](bridgedb_df)
            combined_metadata[source] = tmp_metadata
            if tmp_data.empty:
                warnings.append(f"No annotation available for {source}")
            if not tmp_data.empty:
                combined_data = combine_sources([combined_data, tmp_data])

    return combined_data, combined_metadata
