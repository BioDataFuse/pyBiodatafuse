#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Module for querying StringDB (https://string-db.org/)."""

import datetime
import logging
import warnings
from time import sleep
from typing import Any, Dict, List, Tuple

import pandas as pd
import requests
from requests.exceptions import RequestException

import pyBiodatafuse.constants as Cons
from pyBiodatafuse.utils import get_identifier_of_interest, give_annotator_warning

logger = logging.getLogger("stringdb")

TIMEOUT = 10  # Timeout for requests in seconds


def check_endpoint_stringdb() -> bool:
    """Check the availability of the STRING Db endpoint.

    :returns: True if the endpoint is available, False otherwise.
    """
    try:
        response = requests.get(f"{Cons.STRING_ENDPOINT}/json/version", timeout=TIMEOUT)
        return response.status_code == 200
    except RequestException as e:
        logger.error("Error checking STRING Db endpoint: %s", e)
        return False


def get_version_stringdb() -> dict:
    """Get version of STRING-DB API.

    :returns: a dictionary containing the version information
    """
    try:
        version_call = requests.get(f"{Cons.STRING_ENDPOINT}/json/version", timeout=TIMEOUT).json()
        return {"source_version": version_call[0]["string_version"]}
    except RequestException as e:
        logger.error("Error getting STRING Db version: %s", e)
        return {"source_version": "unknown"}


def _format_data(row, string_ids_df, network_df) -> List[Dict[str, Any]]:
    """Reformat STRING-DB response to match expected output.

    For a given input row (with key 'identifier'), if the network row
    indicates that the input gene appears as preferredName_A then its partner is
    preferredName_B, and vice versa.

    :param row: Row from the input DataFrame (with at least 'identifier' column).
    :param string_ids_df: DataFrame returned from get_string_ids (not used in this version).
    :param network_df: DataFrame returned from the network call.
    :returns: List of dictionaries describing the interactions.
    """
    gene_ppi_links = []
    target_links_set = set()

    target = row[Cons.TARGET_COL]
    identifier = row[Cons.IDENTIFIER_COL]

    for _, row_arr in network_df.iterrows():
        prot_a = row_arr[Cons.STRING_PREFERRED_NAME_A]
        prot_b = row_arr[Cons.STRING_PREFERRED_NAME_B]
        if (prot_a == target or prot_a == identifier) and prot_b not in target_links_set:
            gene_ppi_links.append(
                {
                    Cons.STRING_PPI_INTERACTS_WITH: row_arr[Cons.STRING_PREFERRED_NAME_B],
                    Cons.STRING_GENE_INPUT_ID: row_arr["stringId_B"],
                    Cons.STRING_PPI_SCORE: row_arr["score"],
                    Cons.UNIPROT_TREMBL: row_arr["Uniprot-TrEMBL_B"],
                }
            )
            target_links_set.add(row_arr[Cons.STRING_PREFERRED_NAME_B])

        elif (prot_b == target or prot_b == identifier) and prot_a not in target_links_set:
            gene_ppi_links.append(
                {
                    Cons.STRING_PPI_INTERACTS_WITH: row_arr[Cons.STRING_PREFERRED_NAME_A],
                    Cons.STRING_GENE_INPUT_ID: row_arr["stringId_A"],
                    Cons.STRING_PPI_SCORE: row_arr["score"],
                    Cons.UNIPROT_TREMBL: row_arr["Uniprot-TrEMBL_A"],
                }
            )
            target_links_set.add(row_arr[Cons.STRING_PREFERRED_NAME_A])

    return gene_ppi_links


def get_string_ids(gene_list: List[str], species: str) -> List[Dict[str, Any]]:
    """Get the String identifiers of the gene list.

    :param gene_list: List of gene identifiers
    :param species: Species identifier
    :returns: List of String identifiers
    """
    params = {
        "identifiers": "\r".join(gene_list),  # your protein list
        "species": species,  # species NCBI identifier
        "limit": 1,  # only one (best) identifier per input protein
        "caller_identity": "github.com",  # your app name
    }

    try:
        results = requests.post(
            f"{Cons.STRING_ENDPOINT}/json/get_string_ids", data=params, timeout=TIMEOUT
        ).json()
        return results
    except RequestException as e:
        logger.error(f"Error getting {Cons.STRING} IDs: %s", e)
        return []


def _get_ppi_data(gene_ids: list, species: str) -> pd.DataFrame:
    """Get the String PPI interactions of the gene list.

    :param gene_ids: List of gene identifiers
    :param species: Species identifier
    :returns: DataFrame containing the String PPI interactions
    """
    params = {
        "identifiers": "%0d".join(gene_ids),  # your protein
        "species": species,  # species NCBI identifier
        "caller_identity": "github.com",  # your app name
    }

    try:
        response = requests.post(
            f"{Cons.STRING_ENDPOINT}/json/network", data=params, timeout=TIMEOUT
        ).json()
        return pd.DataFrame(response)
    except RequestException as e:
        logger.error(f"Error getting {Cons.STRING} PPI data: %s", e)
        return pd.DataFrame()


def ensp_to_uniprot(ensp_ids: List[str]) -> Dict[str, str]:
    """Retrieve UniProt IDs from Ensembl protein IDs (ENSP).

    :param ensp_ids: List of Ensembl protein IDs (ENSP)
    :return: Dictionary mapping ENSP IDs to UniProt IDs
    """
    ensp_to_uniprot_map = {}
    headers = {"Content-Type": "application/x-www-form-urlencoded"}
    data = {"from": "Ensembl_Protein", "to": "UniProtKB", "ids": ",".join(ensp_ids)}

    try:
        # Submit the ID mapping request
        response = requests.post(
            f"{Cons.UNIPROT_ID_MAPPER_ENDPOINT}/run", headers=headers, data=data
        )
        response.raise_for_status()
        job_id = response.json()["jobId"]

        # Check the status of the job
        status_url = f"{Cons.UNIPROT_ID_MAPPER_ENDPOINT}/status/{job_id}"
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
        result_url = f"{Cons.UNIPROT_ID_MAPPER_ENDPOINT}/results/{job_id}"
        result_response = requests.get(result_url)
        result_response.raise_for_status()
        results = result_response.json()["results"]

        # Process the results
        for result in results:
            ensp_id = result["from"]
            uniprot_id = result["to"]
            ensp_to_uniprot_map[ensp_id] = uniprot_id

    except Exception as e:
        logger.error(f"Error during {Cons.STRING} Uniprot ID mapping: %s", e)

    return ensp_to_uniprot_map


def get_ppi(
    bridgedb_df: pd.DataFrame, species: str = "human"
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """Annotate genes with protein-protein interactions from STRING-DB.

    :param bridgedb_df: BridgeDb output for creating the list of gene ids to query
    :param species: The species to query. (Try 'Homo sapiens' if 'human' is not working.)
    :returns: a tuple (DataFrame containing the StringDB output, metadata dictionary)
    """
    # Check if the endpoint is available
    if not check_endpoint_stringdb():
        warnings.warn(
            f"{Cons.STRING} endpoint is not available. Unable to retrieve data.", stacklevel=2
        )
        return pd.DataFrame(), {}

    string_version = get_version_stringdb()

    # Record the start time
    start_time = datetime.datetime.now()

    # Retrieve NCBI taxonomy identifier using the given species term
    params = {"db": "taxonomy", "term": species, "retmode": "json"}
    response = requests.get(
        f"{Cons.NCBI_ENDPOINT}/entrez/eutils/esearch.fcgi", params=params
    ).json()
    try:
        species_id = response["esearchresult"]["idlist"][0]
    except (KeyError, IndexError):
        logger.error("NCBI taxonomy search did not return an ID for species: %s", species)
        return pd.DataFrame(), {}

    data_df = get_identifier_of_interest(
        bridgedb_df,
        Cons.STRING_GENE_INPUT_ID,
    ).reset_index(drop=True)
    gene_list = data_df[Cons.TARGET_COL].unique().tolist()
    logger.debug("Gene list: %s", gene_list)

    # Return empty dataframe when only one input is submitted
    if len(gene_list) == 1:
        warnings.warn(
            f"There is only one input gene/protein. Provide at least two input to extract their interactions from {Cons.STRING}.",
            stacklevel=2,
        )
        return pd.DataFrame(), {}

    # Get STRING IDs
    string_ids = get_string_ids(gene_list, species_id)
    if len(string_ids) == 0:
        warnings.warn(
            f"No {Cons.STRING} IDs found for the input genes.",
            stacklevel=2,
        )
        return pd.DataFrame(), {}

    stringdb_ids_df = pd.DataFrame(string_ids)
    stringdb_ids_df.queryIndex = stringdb_ids_df.queryIndex.astype(str)

    # Get the PPI data
    network_df = _get_ppi_data(stringdb_ids_df.stringId.unique().tolist(), species_id)
    logger.debug("Network DataFrame: %s", network_df)

    # Record the end time and build metadata
    end_time = datetime.datetime.now()
    current_date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    time_elapsed = str(end_time - start_time)
    num_new_edges = network_df.drop_duplicates(subset=["stringId_A", "stringId_B"]).shape[0]

    stringdb_metadata: Dict[str, Any] = {
        "datasource": Cons.STRING,
        "metadata": {"source_version": string_version},
        "query": {
            "size": len(gene_list),
            "input_type": Cons.STRING_GENE_INPUT_ID,
            "time": time_elapsed,
            "date": current_date,
            "url": Cons.STRING_ENDPOINT,
        },
    }

    if "stringId_A" not in network_df.columns:
        warnings.warn(
            f"There is no interaction between your input list based on {Cons.STRING}.",
            stacklevel=2,
        )
        return pd.DataFrame(), stringdb_metadata

    # Clean up the network_df
    network_df["stringId_A"] = network_df["stringId_A"].str.split(".").str[1]
    network_df["stringId_B"] = network_df["stringId_B"].str.split(".").str[1]

    # Get UniProt mapping
    ensp_ids = set(network_df["stringId_A"].unique()) | set(network_df["stringId_B"].unique())
    uniprot_map = ensp_to_uniprot(list(ensp_ids))

    # add 'Uniprot-TrEMBL' and 'Uniprot-TrEMBL_link' to network_df
    network_df[Cons.UNIPROT_TREMBL_A] = network_df["stringId_A"].map(uniprot_map)
    network_df[Cons.UNIPROT_TREMBL_B] = network_df["stringId_B"].map(uniprot_map)
    network_df[Cons.UNIPROT_TREMBL_A] = (
        f"{Cons.UNIPROT_TREMBL}:" + network_df[Cons.UNIPROT_TREMBL_A]
    )
    network_df[Cons.UNIPROT_TREMBL_B] = (
        f"{Cons.UNIPROT_TREMBL}:" + network_df[Cons.UNIPROT_TREMBL_B]
    )

    # Format the data
    data_df[Cons.STRING_INTERACT_COL] = data_df.apply(
        lambda row: _format_data(row, stringdb_ids_df, network_df), axis=1
    )

    # Drop rows with no interactions
    data_df = data_df[data_df[Cons.STRING_INTERACT_COL].apply(bool)].reset_index(drop=True)

    if data_df.empty:
        warnings.warn(
            f"There is no interaction between your input list based on {Cons.STRING}.",
            stacklevel=2,
        )
        return pd.DataFrame(), stringdb_metadata

    # Check if the number of new edges is equal to the number of edges in the network_df
    if num_new_edges != len(network_df):
        give_annotator_warning(Cons.STRING_INTERACT_COL)

    stringdb_metadata[Cons.QUERY][Cons.NUM_NODES] = len(gene_list)
    stringdb_metadata[Cons.QUERY][Cons.NUM_EDGES] = num_new_edges

    return data_df, stringdb_metadata
