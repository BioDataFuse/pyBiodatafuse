# coding: utf-8

"""Python utils file for global functions."""

import warnings
from typing import List, Optional

import pandas as pd

from pyBiodatafuse.id_mapper import read_resource_files


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
    subset_df = bridgedb_df[bridgedb_df["target.source"] == db_source]

    return subset_df.reset_index(drop=True)


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


def combine_sources(bridgedb_df: pd.DataFrame, df_list: List[pd.DataFrame]) -> pd.DataFrame:
    """Combine multiple dataframes into a single dataframe.

    :param bridgedb_df: BridgeDb output.
    :param df_list: list of dataframes to be combined.
    :returns: a single dataframe containing from a list of dataframes
    """
    m = bridgedb_df[
        (bridgedb_df["target.source"] == "Ensembl")
        | (bridgedb_df["target.source"] == "PubChem-compound")
    ]
    for df in df_list:
        if not df.empty:
            m = pd.merge(
                m,
                df.drop(columns=["target.source", "identifier.source", "target"], errors="ignore"),
                on="identifier",
                how="outer",
            )

    m = m.loc[:, ~m.columns.duplicated()]  # remove duplicate columns

    # Ensure "Uniprot-TrEMBL" column is present
    # if bridgedb_df["target.source"].eq("Uniprot-TrEMBL").any():
    #    uniprot_trembl_df = bridgedb_df[bridgedb_df["target.source"] == "Uniprot-TrEMBL"]
    #    uniprot_trembl_df = uniprot_trembl_df.groupby("identifier")["target"].apply(list).reset_index()
    #    uniprot_trembl_df.rename(columns={"target": "Uniprot-TrEMBL"}, inplace=True)
    #    m = pd.merge(m, uniprot_trembl_df, on="identifier", how="left")
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


def create_harmonized_input_file(
    annotated_df: pd.DataFrame,
    target_col: str,
    target_source: str,
    identifier_source: Optional[str] = None,
) -> pd.DataFrame:
    """Create a harmonized input DataFrame by extracting specific identifiers from a complex nested structure within a target column.

    :param annotated_df: DataFrame containing the initial data with nested dictionaries.
    :param target_col: Name of the column containing the nested dictionaries.
    :param target_source: The specific identifier source to extract (e.g., 'EFO', 'OMIM').
    :param identifier_source: The main identifier in the output.
    :returns: A DataFrame with original identifiers and the extracted target identifiers.
    """
    harmonized_data = []

    for _i, row in annotated_df.iterrows():
        # Extract the identifier
        if identifier_source is None:
            id = row["identifier"]
            id_source = row["identifier.source"]

        # Extract the the target column
        target_data = row[target_col]

        # Loop through each dictionary in the target data
        for entry in target_data:
            source_idx = entry.get(identifier_source)
            target_idx = entry.get(target_source)

            if source_idx is None or target_idx in [None, ""]:
                continue

            if pd.isna(target_idx) or pd.isna(source_idx):
                continue

            if source_idx.split(":")[1] == "" or target_idx.split(":")[1] == "":
                continue

            id = source_idx.replace(":", "_")
            id_source = identifier_source

            # Extract the specific target identifiers based on the target_source
            for target in target_idx.split(", "):
                # Add a new row to the harmonized data list
                harmonized_data.append(
                    {
                        "identifier": id,
                        "identifier.source": id_source,
                        "target": target.replace(":", "_"),
                        "target.source": target_source,
                    }
                )

    harmonized_df = pd.DataFrame(harmonized_data)

    return harmonized_df.drop_duplicates()
