#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Module for querying StringDB (https://string-db.org/)."""

import datetime
import logging
import warnings
from time import sleep

import numpy as np
import pandas as pd
import requests
from requests.exceptions import RequestException

from pyBiodatafuse.constants import (
    NCBI_ENDPOINT,
    STRING,
    STRING_ENDPOINT,
    STRING_GENE_INPUT_ID,
    STRING_GENE_LINK_ID,
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


def _format_data(row, string_ids_df, network_df):
    """Reformat STRING-DB response to match expected output.

    For a given input row (with key 'identifier'), if the network row
    indicates that the input gene appears as preferredName_A then its partner is
    preferredName_B, and vice versa. The output dictionaries will have the following keys:
      - "stringdb_link_to": the partner gene symbol,
      - "Ensembl": the partner's Ensembl id,
      - "score": the interaction score,
      - "Ensembl_link": the input gene's Ensembl id.

    :param row: Row from the input DataFrame (with at least 'identifier' column).
    :param string_ids_df: DataFrame returned from get_string_ids (not used in this version).
    :param network_df: DataFrame returned from the network call.
    :returns: List of dictionaries describing the interactions.
    """
    gene_ppi_links = []
    target_links_set = set()
    for _, row_arr in network_df.iterrows():
        if row_arr["preferredName_A"]:
            if (
                row_arr["preferredName_A"] == row["target"]
                or row_arr["preferredName_A"] == row["identifier"]
            ) and row_arr["preferredName_B"] not in target_links_set:
                gene_ppi_links.append(
                    {
                        "stringdb_link_to": row_arr["preferredName_B"],
                        STRING_GENE_INPUT_ID: row_arr["stringId_B"].split(".")[1],
                        "score": row_arr["score"],
                        STRING_GENE_LINK_ID: row_arr["stringId_A"].split(".")[1],
                    }
                )
                target_links_set.add(row_arr["preferredName_B"])

            elif (
                row_arr["preferredName_B"] == row["target"]
                or row_arr["preferredName_B"] == row["identifier"]
            ) and row_arr["preferredName_A"] not in target_links_set:
                gene_ppi_links.append(
                    {
                        "stringdb_link_to": row_arr["preferredName_A"],
                        STRING_GENE_INPUT_ID: row_arr["stringId_A"].split(".")[1],
                        "score": row_arr["score"],
                        STRING_GENE_LINK_ID: row_arr["stringId_B"].split(".")[1],
                    }
                )
                target_links_set.add(row_arr["preferredName_A"])

    return gene_ppi_links


def get_string_ids(gene_list: list, species):
    """Get the String identifiers of the gene list."""
    params = {
        "identifiers": "\r".join(gene_list),  # your protein list
        "species": species,  # species NCBI identifier
        "limit": 1,  # only one (best) identifier per input protein
        "caller_identity": "github.com",  # your app name
    }

    try:
        results = requests.post(
            f"{STRING_ENDPOINT}/json/get_string_ids", data=params, timeout=TIMEOUT
        ).json()
        return results
    except RequestException as e:
        logger.error("Error getting STRING IDs: %s", e)
        return []


def _get_ppi_data(gene_ids: list, species) -> pd.DataFrame:
    """Get the String PPI interactions of the gene list."""
    params = {
        "identifiers": "%0d".join(gene_ids),  # your protein
        "species": species,  # species NCBI identifier
        "caller_identity": "github.com",  # your app name
    }

    try:
        response = requests.post(
            f"{STRING_ENDPOINT}/json/network", data=params, timeout=TIMEOUT
        ).json()
        return pd.DataFrame(response)
    except RequestException as e:
        logger.error("Error getting PPI data: %s", e)
        return pd.DataFrame()


def get_ppi(bridgedb_df: pd.DataFrame, species: str = "human"):
    """Annotate genes with protein-protein interactions from STRING-DB.

    :param bridgedb_df: BridgeDb output for creating the list of gene ids to query
    :param species: The species to query. (Try 'Homo sapiens' if 'human' is not working.)
    :returns: a tuple (DataFrame containing the StringDB output, metadata dictionary)
    """
    # Check if the endpoint is available
    if not check_endpoint_stringdb():
        warnings.warn(f"{STRING} endpoint is not available. Unable to retrieve data.", stacklevel=2)
        return pd.DataFrame(), {}

    string_version = get_version_stringdb()

    # Record the start time
    start_time = datetime.datetime.now()

    # Retrieve NCBI taxonomy identifier using the given species term
    params = {"db": "taxonomy", "term": species, "retmode": "json"}
    response = requests.get(f"{NCBI_ENDPOINT}/entrez/eutils/esearch.fcgi", params=params).json()
    try:
        species_id = response["esearchresult"]["idlist"][0]
    except (KeyError, IndexError):
        logger.error("NCBI taxonomy search did not return an ID for species: %s", species)
        return pd.DataFrame(), {}

    data_df = get_identifier_of_interest(
        bridgedb_df,
        STRING_GENE_INPUT_ID,
    ).reset_index(drop=True)
    gene_list = list(set(data_df["target"].tolist()))
    logger.debug("Gene list: %s", gene_list)

    # Return empty dataframe when only one input is submitted
    if len(gene_list) == 1:
        warnings.warn(
            f"There is only one input gene/protein. Provide at least two input to extract their interactions from {STRING}.",
            stacklevel=2,
        )
        return pd.DataFrame(), {}

    # Get STRING IDs
    string_ids = get_string_ids(gene_list, species_id)
    if len(string_ids) == 0:
        return pd.DataFrame(), {}

    stringdb_ids_df = pd.DataFrame(string_ids)
    stringdb_ids_df.queryIndex = stringdb_ids_df.queryIndex.astype(str)

    # Get the PPI data
    network_df = _get_ppi_data(list(stringdb_ids_df.stringId.unique()), species_id)
    logger.debug("Network DataFrame: %s", network_df)

    # Record the end time and build metadata
    end_time = datetime.datetime.now()
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
    data_df[STRING_PPI_COL] = data_df.apply(
        lambda row: _format_data(row, stringdb_ids_df, network_df), axis=1
    )
    # Drop rows where STRING_PPI_COL is an empty list
    data_df = data_df[data_df[STRING_PPI_COL].apply(bool)].reset_index(drop=True)
    # Get matching rows
    # Filter rows with non-empty STRING_PPI_COL
    # rows_with_ppi = data_df[data_df[STRING_PPI_COL].apply(bool)]
    #
    # Iterate over rows with STRING_PPI_COL values
    # for _, row in rows_with_ppi.iterrows():
    #    identifier = row["identifier"]
    #    target = row["target"]
    #    ppi_value = row[STRING_PPI_COL]
    #
    #    # Update rows where target or identifier matches
    #    data_df.loc[
    #        (data_df["identifier"] == target)
    #        | (data_df["target"] == target)
    #        | (data_df["identifier"] == identifier)
    #        | (data_df["target"] == identifier),
    #        STRING_PPI_COL,
    #    ] = data_df[STRING_PPI_COL].apply(
    #        lambda x: ppi_value if not x else (x + ppi_value if #isinstance(x, list) else ppi_value)
    #    )

    # Collect all ENSP IDs from network_df
    ensp_ids = set(network_df["stringId_A"].str.split(".").str[1]) | set(
        network_df["stringId_B"].str.split(".").str[1]
    )
    # Get UniProt mapping
    uniprot_map = ensp_to_uniprot(list(ensp_ids))
    # Add 'Uniprot-TrEMBL' and 'Uniprot-TrEMBL_link' keys to each element in data_df[STRING_PPI_COL]
    data_df[STRING_PPI_COL] = data_df[STRING_PPI_COL].apply(
        lambda lst: (
            [
                {
                    **ppi,
                    "Uniprot-TrEMBL": uniprot_map.get(ppi.get(STRING_GENE_INPUT_ID, "")),
                    "Uniprot-TrEMBL_link": uniprot_map.get(ppi.get(STRING_GENE_LINK_ID, "")),
                }
                for ppi in lst
            ]
            if lst
            else []
        )
    )
    return data_df, string_metadata


def ensp_to_uniprot(ensp_ids):
    """Retrieve UniProt IDs from Ensembl protein IDs (ENSP).

    :param ensp_ids: List of Ensembl protein IDs (ENSP)
    :return: Dictionary mapping ENSP IDs to UniProt IDs
    """
    ensp_to_uniprot_map = {}
    url = "https://rest.uniprot.org/idmapping/run"
    headers = {"Content-Type": "application/x-www-form-urlencoded"}
    data = {"from": "Ensembl_Protein", "to": "UniProtKB", "ids": ",".join(ensp_ids)}

    try:
        # Submit the ID mapping request
        response = requests.post(url, headers=headers, data=data)
        response.raise_for_status()
        job_id = response.json()["jobId"]

        # Check the status of the job
        status_url = f"https://rest.uniprot.org/idmapping/status/{job_id}"
        while True:
            try:
                status_response = requests.get(status_url)
                status_response.raise_for_status()
                status_data = status_response.json()
                if status_data.get("results"):
                    break  # Exit loop if the job is finished
                else:
                    sleep(10)
            except requests.HTTPError as e:
                logger.error("HTTP error occurred: %s", e)
                break
            except Exception as e:
                logger.error("An unexpected error occurred: %s", e)
                break

        # Retrieve the results
        result_url = f"https://rest.uniprot.org/idmapping/results/{job_id}"
        result_response = requests.get(result_url)
        result_response.raise_for_status()
        results = result_response.json()["results"]

        # Process the results
        for result in results:
            ensp_id = result["from"]
            uniprot_id = result["to"]
            ensp_to_uniprot_map[ensp_id] = uniprot_id

    except Exception as e:
        logger.error("Error during ID mapping: %s", e)

    return ensp_to_uniprot_map
