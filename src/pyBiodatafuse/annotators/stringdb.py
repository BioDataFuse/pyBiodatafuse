#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Python file for querying StringDB (https://string-db.org/)."""

from collections import defaultdict
import datetime
import logging
import warnings
import traceback
from time import time, sleep

import numpy as np
import pandas as pd
import requests
from requests.exceptions import RequestException

from pyBiodatafuse.constants import (
    NCBI_ENDPOINT,
    STRING,
    STRING_ENDPOINT,
    STRING_GENE_INPUT_ID,
    STRING_OUTPUT_DICT,
    STRING_PPI_COL,
    STRING_GENE_LINK_ID,
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
    """Reformat STRING-DB response (Helper function).

    :param row: input_df row
    :param string_ids_df: STRING-DB response identifier DataFrame
    :param network_df: STRING-DB response annotation DataFrame
    :returns: StringDB reformatted annotation.
    """
    gene_ppi_links = []
    target_links_set = set()

    # Create a mapping from preferredName to all corresponding queryItem values
    preferredName_to_queryItems = defaultdict(list)
    for _, r in string_ids_df.iterrows():
        preferredName_to_queryItems[r["preferredName"]].append(r["queryItem"])

    for _, row_arr in network_df.iterrows():
        # Check for condition where row_arr["preferredName_A"] matches and identifier fits one possibility
        if row_arr["preferredName_A"] in preferredName_to_queryItems:
            for query_item in preferredName_to_queryItems[row_arr["preferredName_A"]]:
                if row["identifier"] == query_item:
                    # For each possible link corresponding to preferredName_B
                    if row_arr["preferredName_B"] in preferredName_to_queryItems:
                        for link in preferredName_to_queryItems[row_arr["preferredName_B"]]:
                            if row_arr["preferredName_B"] not in target_links_set:
                                gene_ppi_links.append(
                                    {
                                        "stringdb_link_to": link,
                                        STRING_GENE_INPUT_ID: row_arr["stringId_B"].split('.')[1],
                                        "score": row_arr["score"],
                                        "Uniprot-TrEMBL": row_arr["preferredName_A"],
                                    }
                                )
                                target_links_set.add(row_arr["preferredName_B"])
                    # Once a matching identifier is found for A, no need to check other possibilities for this row_arr
                    break

        # Check for condition where row_arr["preferredName_B"] matches and identifier fits one possibility
        if row_arr["preferredName_B"] in preferredName_to_queryItems:
            for query_item in preferredName_to_queryItems[row_arr["preferredName_B"]]:
                if row["identifier"] == query_item:
                    # For each possible link corresponding to preferredName_A
                    if row_arr["preferredName_A"] in preferredName_to_queryItems:
                        for link in preferredName_to_queryItems[row_arr["preferredName_A"]]:
                            if row_arr["preferredName_A"] not in target_links_set:
                                gene_ppi_links.append(
                                    {
                                        "stringdb_link_to": link,
                                        STRING_GENE_INPUT_ID: row_arr["stringId_A"].split('.')[1],
                                        "score": row_arr["score"],
                                        "Uniprot-TrEMBL": row_arr["preferredName_B"],
                                    }
                                )
                                target_links_set.add(row_arr["preferredName_A"])
                    # Once a matching identifier is found for B, break out of loop for this row_arr
                    break

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
    :param species: The species to query. All species that are supported by both NCBI and STRINGDB can be used.
    :returns: a DataFrame containing the StringDB output and dictionary of the metadata.
    """
    # Check if the endpoint is available
    if not check_endpoint_stringdb():
        warnings.warn(f"{STRING} endpoint is not available. Unable to retrieve data.", stacklevel=2)
        return pd.DataFrame(), {}

    string_version = get_version_stringdb()

    # Record the start time
    start_time = datetime.datetime.now()

    # Retrieve NCBI taxonomy identifier
    params = {"db": "taxonomy", "term": species, "retmode": "json"}
    response = requests.get(f"{NCBI_ENDPOINT}/entrez/eutils/esearch.fcgi", params=params).json()
    species_id = response["esearchresult"]["idlist"][0]

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
    string_ids = get_string_ids(gene_list, species_id)
    if len(string_ids) == 0:
        return pd.DataFrame(), {}

    stringdb_ids_df = pd.DataFrame(string_ids)
    stringdb_ids_df.queryIndex = stringdb_ids_df.queryIndex.astype(str)

    # Chunk the string IDs into groups of 2000
    unique_string_ids = list(stringdb_ids_df.stringId.unique())
    all_network_dfs = []
    chunk_size = 1500

    # Get the PPI data
    for i in range(0, len(unique_string_ids), chunk_size):
        chunk = unique_string_ids[i:i + chunk_size]
        chunk_network_df = _get_ppi_data(chunk, species_id)
        if not chunk_network_df.empty:
            all_network_dfs.append(chunk_network_df)
    
    # Combine all chunked results
    if all_network_dfs:
        network_df = pd.concat(all_network_dfs, ignore_index=True)
    else:
        network_df = pd.DataFrame()
                                   
    # network_df = _get_ppi_data(list(stringdb_ids_df.stringId.unique()), species_id)

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
    data_df[STRING_PPI_COL] = data_df.apply(
        lambda row: _format_data(row, stringdb_ids_df, network_df), axis=1
    )
    data_df[STRING_PPI_COL] = data_df[STRING_PPI_COL].apply(
        lambda x: ([{key: np.nan for key in STRING_OUTPUT_DICT}] if len(x) == 0 else x)
    )
    # Collect all ENSP IDs from network_df
    ensp_ids = set(network_df['stringId_A'].str.split('.').str[1]) | set(network_df['stringId_B'].str.split('.').str[1])
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
                if lst
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
    data = {
        "from": "Ensembl_Protein",
        "to": "UniProtKB",
        "ids": ",".join(ensp_ids)
    }

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
