#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Python file for queriying StringDB (https://string-db.org/)."""

import datetime
import io
import logging
import warnings

import pandas as pd
import requests

from pyBiodatafuse.constants import STRING, STRING_ENDPOINT
from pyBiodatafuse.utils import get_identifier_of_interest

logger = logging.getLogger("stringdb")


def check_endpoint_stringdb() -> bool:
    """Check the availability of the STRING Db endpoint.

    :returns: True if the endpoint is available, False otherwise.
    """
    response = requests.get(f"{STRING_ENDPOINT}/json/version")

    # Check if API is down
    if response.status_code == 200:
        return True
    else:
        return False


def get_version_stringdb() -> dict:
    """Get version of STRING-DB API.

    :returns: a dictionary containing the version information
    """
    version_call = requests.get(f"{STRING_ENDPOINT}/json/version")
    return version_call.json()


def _format_data(row, network_df):
    """Reformat STRING-DB response (Helper function).

    :param row: input_df row
    :param network_df: STRING-DB response annotation DataFrame
    :returns: StringDB reformatted annotation.
    """
    gene_ppi_links = list()

    target_links_set = set()

    for _i, row_arr in network_df.iterrows():
        if row_arr["preferredName_A"] == row["identifier"]:
            if row_arr["preferredName_B"] not in target_links_set:
                gene_ppi_links.append(
                    {"stringdb_link_to": row_arr["preferredName_B"], "score": row_arr["score"]}
                )
                target_links_set.add(row_arr["preferredName_B"])

        elif row_arr["preferredName_B"] == row["identifier"]:
            if row_arr["preferredName_A"] not in target_links_set:
                gene_ppi_links.append(
                    {"stringdb_link_to": row_arr["preferredName_A"], "score": row_arr["score"]}
                )
                target_links_set.add(row_arr["preferredName_A"])

    return gene_ppi_links


def get_ppi(bridgedb_df: pd.DataFrame):
    """Annotate genes with protein-protein interactions from STRING-DB.

    :param bridgedb_df: BridgeDb output for creating the list of gene ids to query
    :returns: a DataFrame containing the StringDB output and dictionary of the metadata.
    """
    # Check if the IDSM endpoint is available
    api_available = check_endpoint_stringdb()
    if not api_available:
        warnings.warn(f"{STRING} endpoint is not available. Unable to retrieve data.", stacklevel=2)
        return pd.DataFrame(), {}

    # Record the start time
    start_time = datetime.datetime.now()

    data_df = get_identifier_of_interest(bridgedb_df, "Ensembl")
    data_df = data_df.reset_index(drop=True)

    gene_list = list(set(data_df["target"].tolist()))

    # --------- Get the String identifiers of the gene list --------#
    output_format = "tsv"
    method = "get_string_ids"

    params = {
        "identifiers": "\r".join(gene_list),  # your protein list
        "species": 9606,  # species NCBI identifier
        "limit": 1,  # only one (best) identifier per input protein
        "caller_identity": "github.com",  # your app name
    }

    request_url = "/".join([STRING_ENDPOINT, output_format, method])

    results = requests.post(request_url, data=params)
    results_contents = results.content.decode("utf-8")

    stringdb_ids_df = pd.read_csv(io.StringIO(results_contents), sep="\t")
    stringdb_ids_df.queryIndex = stringdb_ids_df.queryIndex.astype(str)

    # ---------- Get String PPI network using the String identifiers ---------------#

    method = "network"
    request_url = "/".join([STRING_ENDPOINT, output_format, method])

    params = {
        "identifiers": "%0d".join(list(stringdb_ids_df.stringId.unique())),  # your protein
        "species": 9606,  # species NCBI identifier
        "caller_identity": "github.com",  # your app name
    }

    response = requests.post(request_url, data=params)
    results_contents = response.content.decode("utf-8")

    network_df = pd.read_csv(io.StringIO(results_contents), sep="\t")

    # ---------- Add the interactions of each protein (row) to a new column ('stringdb') ---------------#

    data_df[STRING] = data_df.apply(_format_data, network_df=network_df, axis=1)

    # Record the end time
    end_time = datetime.datetime.now()

    # Metdata details
    # Get the current date and time
    current_date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    # Calculate the time elapsed
    time_elapsed = str(end_time - start_time)
    # Add version to metadata file

    string_version = get_version_stringdb()

    # Add the datasource, query, query time, and the date to metadata
    string_metadata = {
        "datasource": STRING,
        "metadata": {"source_version": string_version},
        "query": {
            "size": len(gene_list),
            "time": time_elapsed,
            "date": current_date,
            "url": STRING_ENDPOINT,
        },
    }

    return data_df, string_metadata
