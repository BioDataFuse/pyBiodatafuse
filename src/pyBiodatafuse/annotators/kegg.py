#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Python file for queriying StringDB (https://rest.kegg.jp/)."""

import datetime
import logging
import warnings

import numpy as np
import pandas as pd
import requests
import time

from pyBiodatafuse.constants import (
    KEGG,
    KEGG_ENDPOINT,
    KEGG_GENE_INPUT_ID,
    KEGG_COL_NAME
)

from pyBiodatafuse.utils import check_columns_against_constants, get_identifier_of_interest

def check_endpoint_kegg() -> dict:
    """Check if the endpoint of the KEGG API is available.

    :returns: a dictionary containing the version information
    """
    response = requests.get(f"{KEGG_ENDPOINT}/info/kegg")
    # Check if API is down
    if response.status_code == 200:
        return True
    else:
        return False
    

def get_kegg_ids(row):
    """Get the KEGG identifiers of the gene list.
    """
    print(row)
    results = requests.get(f"{KEGG_ENDPOINT}/conv/genes/ncbi-geneid:{row['target']}")
    kegg_id = results.text.split()
    return {"KEGG_id": kegg_id[1]}


def get_pathway_link(row):
    kegg_dict = row[KEGG_COL_NAME]
    results = requests.get(f"{KEGG_ENDPOINT}/link/pathway/{kegg_dict.get('KEGG_id')}")

    pathways = []
    for line in results.text.strip().split("\n"):
        pathway_info = {}
        parts = line.split("\t")
        pathway_id = parts[1]
        pathway_info["pathway_id"] = pathway_id 

        results_kgml = requests.get(f"{KEGG_ENDPOINT}/get/{pathway_id}/kgml")
        title_start = results_kgml.text.find('title="') + len('title="')
        title_end = results_kgml.text.find('"', title_start)

        # Extract the title substring
        pathway_title = results_kgml.text[title_start:title_end]
        pathway_info["pathway_label"] = pathway_title
        pathways.append(pathway_info)


    kegg_dict["pathways"] = pathways

    return kegg_dict
    

def get_pathways(bridgedb_df):    
    
    api_available = check_endpoint_kegg()
    if not api_available:
        warnings.warn(f"{KEGG} endpoint is not available. Unable to retrieve data.", stacklevel=2)
        return pd.DataFrame(), {}

    # Record the start time
    # start_time = datetime.datetime.now()

    data_df = get_identifier_of_interest(bridgedb_df, KEGG_GENE_INPUT_ID)
    data_df = data_df.reset_index(drop=True)

    gene_list = list(set(data_df["target"].tolist()))

    # Get the KEGG identifiers
    data_df[KEGG_COL_NAME] = data_df.apply(lambda row: get_kegg_ids(row), axis=1)
    print(data_df)

    # Get the links for the KEGG pathways
    data_df[KEGG_COL_NAME] = data_df.apply(lambda row: get_pathway_link(row), axis=1)
    print(data_df)

    return data_df
    

    
    
