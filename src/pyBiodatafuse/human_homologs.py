#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Python file for queriying Ensembl to get human homologs for mouse genes."""

import datetime
import warnings

import numpy as np
import pandas as pd
import requests

import pyBiodatafuse.constants as Cons
from pyBiodatafuse.utils import get_identifier_of_interest


def check_endpoint_ensembl() -> bool:
    """Check if the endpoint of the Ensembl API is available.

    :returns: A True statement if the endpoint is available, else return False
    """
    response = requests.get(f"{Cons.ENSEMBL_ENDPOINT}/info/ping")
    # Check if API is down
    if response.status_code == 200:
        return True
    else:
        return False


def check_version_ensembl() -> str:
    """Check the current version of the REST API.

    :returns: A True statement if the endpoint is available, else return False
    """
    response = requests.get(
        f"{Cons.ENSEMBL_ENDPOINT}/info/rest", headers={"Content-Type": "application/json"}
    )
    # Check if API is down
    return response.text


def get_human_homologs(row):
    """Retrieve human homologs for mouse genes using Ensembl API.

    :param row: row from input dataframe.
    :returns: dictionary mapping mouse genes to human homologs.
    """
    response = requests.get(
        f"{Cons.ENSEMBL_ENDPOINT}/homology/id/mouse/{row['target']}",
        headers={"Content-Type": "application/json"},
        params={"target_species": "homo_sapiens"},
    )

    if response.status_code != 200:
        return [{Cons.ENSEMBL_HOMOLOG_MAIN_LABEL: np.nan}]

    data = response.json()
    if "data" in data and len(data["data"]) > 0:
        for homology in data["data"][0].get("homologies", []):
            if homology["target"]["species"] == "homo_sapiens":
                homolog = homology["target"]["id"]
                return [{Cons.ENSEMBL_HOMOLOG_MAIN_LABEL: homolog}]

    return [{Cons.ENSEMBL_HOMOLOG_MAIN_LABEL: np.nan}]


def get_homologs(bridgedb_df):
    """Retrieve homologs for input DataFrame.

    :param bridgedb_df: input dataframe.
    :returns: dataframe including the human homologs as well as the metadata.
    """
    api_available = check_endpoint_ensembl()
    if not api_available:
        warnings.warn(
            f"{Cons.ENSEMBL} endpoint is not available. Unable to retrieve data.", stacklevel=2
        )
        return pd.DataFrame(), {}

    ensembl_version = check_version_ensembl()

    # Record the start time
    start_time = datetime.datetime.now()

    data_df = get_identifier_of_interest(bridgedb_df, Cons.ENSEMBL_GENE_INPUT_ID)
    data_df = data_df.reset_index(drop=True)
    gene_list = list(set(data_df[Cons.TARGET_COL].tolist()))

    # Get the human homologs
    data_df[Cons.ENSEMBL_HOMOLOG_COL] = data_df.apply(lambda row: get_human_homologs(row), axis=1)

    # Record the end time
    end_time = datetime.datetime.now()

    """Metadata details"""
    # Get the current date and time
    current_date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    # Calculate the time elapsed
    time_elapsed = str(end_time - start_time)
    # Calculate new edges
    num_new_edges = data_df.shape[0]

    # Add the datasource, query, query time, and the date to metadata
    kegg_metadata = {
        "datasource": Cons.ENSEMBL,
        "metadata": {"source_version": ensembl_version},
        "query": {
            "size": len(gene_list),
            "input_type": Cons.ENSEMBL_GENE_INPUT_ID,
            "number_of_added_edges": num_new_edges,
            "time": time_elapsed,
            "date": current_date,
            "url": Cons.ENSEMBL_ENDPOINT,
        },
    }

    return data_df, kegg_metadata
