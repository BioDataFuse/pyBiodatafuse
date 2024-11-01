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

    :returns: A True statement if the endpoint is available, else return False
    """
    response = requests.get(f"{KEGG_ENDPOINT}/info/kegg")
    # Check if API is down
    if response.status_code == 200:
        return True
    else:
        return False


def check_version_kegg() -> dict:
    """Check the current version of the KEGG databass

    :returns: a dictionary containing the version information
    """
    response = requests.get(f"{KEGG_ENDPOINT}/info/kegg")
    for line in response.text.splitlines():
        if "Release" in line:
            release_version = line.split()[2]
            release_version = release_version.rstrip(",")
            print("KEGG Release Version:", release_version)
            break



def get_kegg_ids(row):
    """Get the KEGG identifiers of the gene list.

    :param row: input_df row
    :returns: a dictionary containing the KEGG identifier
    """
    print(row)
    results = requests.get(f"{KEGG_ENDPOINT}/conv/genes/ncbi-geneid:{row['target']}")
    kegg_id = results.text.split()
    return {"KEGG_id": kegg_id[1]}


def get_pathway_link(row):
    """Get pathway information for the input genes.
    
    :param row: input_df row
    :returns: Dictionary containing pathway IDs and labels.
    """
    kegg_dict = row[KEGG_COL_NAME]
    results = requests.get(f"{KEGG_ENDPOINT}/link/pathway/{kegg_dict.get('KEGG_id')}")

    pathways = []
    for line in results.text.strip().split("\n"):
        pathway_info = {}
        parts = line.split("\t")
        pathway_id = parts[1]
        pathway_info["pathway_id"] = pathway_id 

        # Get KGML file from KEGG API
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
    
    check_version_kegg()

    # Record the start time
    start_time = datetime.datetime.now()

    data_df = get_identifier_of_interest(bridgedb_df, KEGG_GENE_INPUT_ID)
    data_df = data_df.reset_index(drop=True)
    gene_list = list(set(data_df["target"].tolist()))

    # Get the KEGG identifiers
    data_df[KEGG_COL_NAME] = data_df.apply(lambda row: get_kegg_ids(row), axis=1)
    print(data_df)

    # Get the links for the KEGG pathways
    data_df[KEGG_COL_NAME] = data_df.apply(lambda row: get_pathway_link(row), axis=1)
    print(data_df)

    # Record the end time
    end_time = datetime.datetime.now()

    """Metadata details"""
    # Get the current date and time
    current_date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    # Calculate the time elapsed
    time_elapsed = str(end_time - start_time)

    # Add the datasource, query, query time, and the date to metadata
    # string_metadata = {
    #     "datasource": KEGG,
    #     "metadata": {"source_version": string_version},
    #     "query": {
    #         "size": len(gene_list),
    #         "input_type": KEGG_GENE_INPUT_ID,
    #         "number_of_added_edges": num_new_edges,
    #         "time": time_elapsed,
    #         "date": current_date,
    #         "url": KEGG_ENDPOINT,
    #     },
    # }

    return data_df
    

    
    
