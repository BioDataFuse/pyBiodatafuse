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
    """Get version of KEGG API.

    :returns: a dictionary containing the version information
    """
    response = requests.get(f"{KEGG_ENDPOINT}/info/kegg")
    # Check if API is down
    if response.status_code == 200:
        return True
    else:
        return False
    
def get_kegg_ids(gene_list):
    """Get the String identifiers of the gene list.
    """
    for gene in gene_list:
        print(gene)
        results = requests.get(f"https://rest.kegg.jp/conv/genes/ncbi-geneid:{gene}")
        # pathway_data = requests.get(f"{KEGG_ENDPOINT}/get/{gene}")

        print(results.text)
    
    
def get_info(bridgedb_df):    
    
    api_available = check_endpoint_kegg()
    if not api_available:
        warnings.warn(f"{KEGG} endpoint is not available. Unable to retrieve data.", stacklevel=2)
        return pd.DataFrame(), {}

    # Record the start time
    # start_time = datetime.datetime.now()

    data_df = get_identifier_of_interest(bridgedb_df, KEGG_GENE_INPUT_ID)
    data_df = data_df.reset_index(drop=True)

    gene_list = list(set(data_df["target"].tolist()))

    # data_df[KEGG_COL_NAME] = data_df.apply(lambda row: _format_data(row, stringdb_ids_df, network_df, to_uniprot), axis=1)

    pathway_data = get_kegg_ids(gene_list)
    print(pathway_data)
    
    
