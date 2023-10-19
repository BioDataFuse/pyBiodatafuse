# coding: utf-8

"""Python utils file for global functions."""

import json
import os
from typing import List

import pandas as pd

from pyBiodatafuse.constants import COMBINED_DIR, DATA_DIR, RESOURCES_DIR


def get_identifier_of_interest(bridgedb_df: pd.DataFrame, source: str) -> pd.DataFrame:
    """Get identifier of interest from BridgeDb output file.

    :param bridgedb_df: DataFrame containing the output from BridgeDb
    :param source: identifier of interest from BridgeDB (e.g. "NCBI Gene")
    :returns: a DataFrame containing the identifiers of interest
    """
    # Load identifier options
    identifier_options = pd.read_csv(f"{RESOURCES_DIR}/datasources.Csv")["source"].tolist()

    # Check if source is in identifier options
    assert source in identifier_options, f"Source {source} is not in identifier options"

    # Filter rows where "target.source" is specific datasource "NCBI Gene"
    return bridgedb_df[bridgedb_df["target.source"] == source]


def create_or_append_to_metadata(data: dict) -> None:
    """Create and/or append data to a metadata file.

    :param data: dictionary of data to be saved to the metadata file.
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
    """
    metadata_file_path = f"{DATA_DIR}/metadata.json"

    if os.path.exists(metadata_file_path):
        existing_data = json.load(open(metadata_file_path))
        prev_sources = [data["datasource"] for data in existing_data]
    else:
        # Initialize the existing data to an empty list
        existing_data = []
        prev_sources = []

    assert isinstance(data, dict), "Unsupported data type. Only dict is supported."

    if data["datasource"] not in prev_sources:
        existing_data.append(data)

    # Save the updated data back to the file
    with open(metadata_file_path, "w") as file:
        json.dump(list(existing_data), file, indent=4, ensure_ascii=False)


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
    merged_df = (
        merged_df.groupby(["identifier", "identifier.source", "target", "target.source"])[col_name]
        .agg(sum)
        .reset_index()
    )

    return merged_df


def combine_sources(df_list: List[pd.DataFrame]) -> pd.DataFrame:
    """Combine multiple dataframes into a single dataframe.

    :param df_list: list of dataframes to be combined
    :returns: a single dataframe containing from a list of dataframes
    """
    m = pd.concat(df_list, axis=1)
    m = m.loc[:, ~m.columns.duplicated()]  # remove duplicate columns

    # Create the directory to save combined data
    os.makedirs(COMBINED_DIR, exist_ok=True)

    m.to_csv(f"{COMBINED_DIR}/combined_annotated_data.tsv", index=False, sep="\t")
    return m
