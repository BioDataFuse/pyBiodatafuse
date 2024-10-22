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
    return {kegg_id[1]}


def get_pathway_link(kegg_id_dict):
    for key in kegg_id_dict:
        print(kegg_id_dict[key]) # debug
        results = requests.get(f"{KEGG_ENDPOINT}/link/pathway/{kegg_id_dict[key]}")
        print(results.text) # debug
    

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


    data_df[KEGG_COL_NAME] = data_df.apply(lambda row: get_kegg_ids(row), axis=1)

    # Get the KEGG identifiers
    # for gene in gene_list:
    #     kegg_id = get_kegg_ids(gene)
    #     data_df[KEGG_COL_NAME]

    # Get the links for the KEGG pathways
    # get_pathway_link(kegg_id_dict)
    

    
    
