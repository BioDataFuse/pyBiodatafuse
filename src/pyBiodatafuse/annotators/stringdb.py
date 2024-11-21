#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Python file for querying StringDB (https://string-db.org/)."""

import datetime
import logging
import warnings

import numpy as np
import pandas as pd
import requests
from requests.exceptions import RequestException

from pyBiodatafuse.constants import (
    STRING,
    STRING_ENDPOINT,
    STRING_GENE_INPUT_ID,
    STRING_OUTPUT_DICT,
    STRING_PPI_COL,
)
from pyBiodatafuse.utils import check_columns_against_constants, get_identifier_of_interest

logger = logging.getLogger("stringdb")

TIMEOUT = 10  # Timeout for requests in seconds


def check_endpoint_stringdb() -> bool:
    """Check the availability of the STRING Db endpoint.

    :returns: True if the endpoint is available, False otherwise.
    """
    try:
        response = requests.get(f"{STRING_ENDPOINT}/json/version", timeout=TIMEOUT)
        return response.status_code == 200
    except RequestException as e:
        logger.error("Error checking STRING Db endpoint: %s", e)
        return False


def get_version_stringdb() -> dict:
    """Get version of STRING-DB API.

    :returns: a dictionary containing the version information
    """
    try:
        version_call = requests.get(f"{STRING_ENDPOINT}/json/version", timeout=TIMEOUT).json()
        return {"source_version": version_call[0]["string_version"]}
    except RequestException as e:
        logger.error("Error getting STRING Db version: %s", e)
        return {"source_version": "unknown"}


def _format_data(row, network_df):
    """Reformat STRING-DB response (Helper function).

    :param row: input_df row
    :param network_df: STRING-DB response annotation DataFrame
    :returns: StringDB reformatted annotation.
    """
    gene_ppi_links = []
    target_links_set = set()

    for _, row_arr in network_df.iterrows():
        if row_arr["preferredName_A"] == row["identifier"] and row_arr["preferredName_B"] not in target_links_set:
            gene_ppi_links.append({
                "stringdb_link_to": row_arr["preferredName_B"],
                STRING_GENE_INPUT_ID: f"{STRING_GENE_INPUT_ID}:{row_arr['stringId_B'].split('.')[1]}",
                "score": row_arr["score"],
            })
            target_links_set.add(row_arr["preferredName_B"])

        elif row_arr["preferredName_B"] == row["identifier"] and row_arr["preferredName_A"] not in target_links_set:
            gene_ppi_links.append({
                "stringdb_link_to": row_arr["preferredName_A"],
                STRING_GENE_INPUT_ID: f"{STRING_GENE_INPUT_ID}:{row_arr['stringId_A'].split('.')[1]}",
                "score": row_arr["score"],
            })
            target_links_set.add(row_arr["preferredName_A"])

    return gene_ppi_links


def get_string_ids(gene_list: list):
    """Get the String identifiers of the gene list."""
    params = {
        "identifiers": "\r".join(gene_list),  # your protein list
        "species": 9606,  # species NCBI identifier
        "limit": 1,  # only one (best) identifier per input protein
        "caller_identity": "github.com",  # your app name
    }

    try:
        results = requests.post(f"{STRING_ENDPOINT}/json/get_string_ids",
                                data=params, timeout=TIMEOUT).json()
        return results
    except RequestException as e:
        logger.error("Error getting STRING IDs: %s", e)
        return []


def _get_ppi_data(gene_ids: list) -> pd.DataFrame:
    """Get the String PPI interactions of the gene list."""
    params = {
        "identifiers": "%0d".join(gene_ids),  # your protein
        "species": 9606,  # species NCBI identifier
        "caller_identity": "github.com",  # your app name
    }

    try:
        response = requests.post(f"{STRING_ENDPOINT}/json/network", data=params, timeout=TIMEOUT).json()
        return pd.DataFrame(response)
    except RequestException as e:
        logger.error("Error getting PPI data: %s", e)
        return pd.DataFrame()


def get_ppi(bridgedb_df: pd.DataFrame):
    """Annotate genes with protein-protein interactions from STRING-DB.

    :param bridgedb_df: BridgeDb output for creating the list of gene ids to query
    :returns: a DataFrame containing the StringDB output and dictionary of the metadata.
    """
    # Check if the endpoint is available
    if not check_endpoint_stringdb():
        warnings.warn(f"{STRING} endpoint is not available. Unable to retrieve data.", stacklevel=2)
        return pd.DataFrame(), {}

    string_version = get_version_stringdb()

    # Record the start time
    start_time = datetime.datetime.now()

    data_df = get_identifier_of_interest(bridgedb_df, STRING_GENE_INPUT_ID).reset_index(drop=True)
    gene_list = list(set(data_df["target"].tolist()))

    # Return empty dataframe when only one input submitted
    if len(gene_list) == 1:
        warnings.warn(
            f"There is only one input gene/protein. Provide at least two input to extract their interactions from {STRING}.",
            stacklevel=2,
        )
        return pd.DataFrame(), {}

    # Get ids
    string_ids = get_string_ids(gene_list)
    if string_ids.empty:
        return pd.DataFrame(), {}

    stringdb_ids_df = pd.DataFrame(string_ids)
    stringdb_ids_df.queryIndex = stringdb_ids_df.queryIndex.astype(str)

    # Get the PPI data
    network_df = _get_ppi_data(list(stringdb_ids_df.stringId.unique()))

    # Record the end time
    end_time = datetime.datetime.now()

    # Metadata details
    current_date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    time_elapsed = str(end_time - start_time)
    num_new_edges = network_df.drop_duplicates(subset=["stringId_A", "stringId_B"]).shape[0]

    if num_new_edges != len(network_df):
        warnings.warn(
            f"The network_df in {STRING} annotator should be checked, please create an issue https://github.com/BioDataFuse/pyBiodatafuse/issues/.",
            stacklevel=2,
        )

    string_metadata = {
        "datasource": STRING,
        "metadata": {"source_version": string_version},
        "query": {
            "size": len(gene_list),
            "input_type": STRING_GENE_INPUT_ID,
            "number_of_added_edges": num_new_edges,
            "time": time_elapsed,
            "date": current_date,
            "url": STRING_ENDPOINT,
        },
    }

    if "stringId_A" not in network_df.columns:
        warnings.warn(
            f"There is no interaction between your input list based on {STRING}, {string_version}.",
            stacklevel=2,
        )
        return pd.DataFrame(), string_metadata

    # Format the data
    data_df[STRING_PPI_COL] = data_df.apply(_format_data, network_df=network_df, axis=1)
    data_df[STRING_PPI_COL] = data_df[STRING_PPI_COL].apply(
        lambda x: ([{key: np.nan for key in STRING_OUTPUT_DICT}] if len(x) == 0 else x)
    )

    # Check if all keys in df match the keys in OUTPUT_DICT
    exploded_df = data_df.explode(STRING_PPI_COL)
    ppi_df = pd.json_normalize(exploded_df[STRING_PPI_COL])
    check_columns_against_constants(
        data_df=ppi_df,
        output_dict=STRING_OUTPUT_DICT,
        check_values_in=[],
    )

    return data_df, string_metadata
