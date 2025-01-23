#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Python file for queriying KEGG (https://rest.kegg.jp/)."""

import datetime
import json
import logging
import time
import warnings

import numpy as np
import pandas as pd
import requests

from pyBiodatafuse.constants import KEGG, KEGG_COL, KEGG_ENDPOINT, KEGG_GENE_INPUT_ID
from pyBiodatafuse.utils import check_columns_against_constants, get_identifier_of_interest


def check_endpoint_kegg() -> bool:
    """Check if the endpoint of the KEGG API is available.

    :returns: A True statement if the endpoint is available, else return False
    """
    response = requests.get(f"{KEGG_ENDPOINT}/info/kegg")
    # Check if API is down
    if response.status_code == 200:
        return True
    else:
        return False


def check_version_kegg() -> str:
    """Check the current version of the KEGG database.

    :returns: a dictionary containing the version information
    """
    response = requests.get(f"{KEGG_ENDPOINT}/info/kegg")
    for line in response.text.splitlines():
        if "Release" in line:
            release_version = line.split()[2]
            release_version = release_version.rstrip(",")
            return release_version
        else:
            release_version = "Error: Release version not found."


def get_kegg_ids(row) -> dict:
    """Get the KEGG identifiers of the gene list.

    :param row: input_df row
    :returns: a dictionary containing the KEGG identifier
    """
    results = requests.get(f"{KEGG_ENDPOINT}/conv/genes/ncbi-geneid:{row['target']}")
    kegg_id = results.text.split()
    if len(kegg_id) > 1:
        return {"KEGG_id": kegg_id[1]}
    else:
        return {"KEGG_id": pd.nan}


def get_compound_genes(pathway_info, results_entry):
    """Get compounds and gene counts from a pathway.

    :param pathway_info: Dictionary containing all information of the pathway
    :param results_entry: KGML file from which further information gets extracted
    :returns: Dictionary containing compounds and gene count
    """
    genes = []
    compounds = []  # Initialize an empty list to hold compound dictionaries
    section = None

    for line in results_entry.text.splitlines():
        current_compound = {}

        if line.startswith("GENE"):
            section = "GENE"
        elif line.startswith("COMPOUND"):
            section = "COMPOUND"
        elif line.startswith("REFERENCE"):
            section = None

        if section == "GENE":
            parts = line.split()
            if len(parts) > 0 and parts[0].isdigit():  # Gene identifier is numeric
                genes.append(parts[0])

        elif section == "COMPOUND":
            parts = line.split()
            for part in parts:
                if (
                    part.startswith("C") and part[1:].isdigit()
                ):  # KEGG compound identifiers start with C
                    current_compound["KEGG_identifier"] = part

                    # Initialize "name" to None in case it's not found
                    current_compound["name"] = None

                    # Query KEGG API for the compound details
                    results = requests.get(f"{KEGG_ENDPOINT}/get/{part}")

                    # Parse the compound name
                    for line in results.text.splitlines():
                        if line.startswith("NAME"):
                            parts = line.split()
                            current_compound["name"] = parts[1] if len(parts) > 1 else None
                            # Remove any trailing semicolons from the name
                            if current_compound["name"]:
                                current_compound["name"] = current_compound["name"].rstrip(";")
                            break

                    compounds.append(current_compound)

    # Set a default structure if no compounds were found
    if not compounds:
        compounds = [{"KEGG_identifier": None, "name": None}]

    pathway_info["gene_count"] = len(genes)
    pathway_info["compounds"] = compounds

    return pathway_info


def get_pathway_info(row):
    """Get pathway information for the input genes.

    :param row: input_df row
    :returns: Dictionary containing pathway IDs and labels.
    """
    kegg_dict = row[KEGG_COL]
    if kegg_dict is None or not isinstance(kegg_dict, dict) or kegg_dict.get("KEGG_id") is None:
        # Return default structure with np.nan
        return {
            "KEGG_id": np.nan,
            "pathways": [
                {
                    "pathway_id": np.nan,
                    "pathway_label": np.nan,
                    "gene_count": np.nan,
                    "compounds": [{"KEGG_identifier": None, "name": None}],
                }
            ],
        }

    results = requests.get(f"{KEGG_ENDPOINT}/link/pathway/{kegg_dict.get('KEGG_id')}")
    if len(results.text) <= 1:
        # Return default structure with np.nan when no pathways are found
        kegg_dict["pathways"] = [
            {
                "pathway_id": np.nan,
                "pathway_label": np.nan,
                "gene_count": np.nan,
                "compounds": [{"KEGG_identifier": None, "name": None}],
            }
        ]
        return kegg_dict

    pathways = []

    for line in results.text.strip().split("\n"):
        pathway_info = {}
        parts = line.split("\t")
        pathway_id = parts[1]
        pathway_info["pathway_id"] = pathway_id

        # Get entry from KEGG API
        results_entry = requests.get(f"{KEGG_ENDPOINT}/get/{pathway_id}")

        for line in results_entry.text.splitlines():
            if line.startswith("NAME"):
                pathway_info["pathway_label"] = line.split("  ", 1)[1].strip()
                break

        pathway_info = get_compound_genes(pathway_info, results_entry)
        pathways.append(pathway_info)

    kegg_dict["pathways"] = (
        pathways
        if pathways
        else [
            {
                "pathway_id": np.nan,
                "pathway_label": np.nan,
                "gene_count": np.nan,
                "compounds": [{"KEGG_identifier": None, "name": None}],
            }
        ]
    )
    return kegg_dict


def get_pathways(bridgedb_df):
    """Annotate genes with KEGG pathway information."""
    api_available = check_endpoint_kegg()
    if not api_available:
        warnings.warn(f"{KEGG} endpoint is not available. Unable to retrieve data.", stacklevel=2)
        return pd.DataFrame(), {}

    kegg_version = check_version_kegg()

    # Record the start time
    start_time = datetime.datetime.now()

    data_df = get_identifier_of_interest(bridgedb_df, KEGG_GENE_INPUT_ID)
    data_df = data_df.reset_index(drop=True)
    gene_list = list(set(data_df["target"].tolist()))

    # Get the KEGG identifiers
    data_df[KEGG_COL] = data_df.apply(lambda row: get_kegg_ids(row), axis=1)

    # Get the links for the KEGG pathways
    data_df[KEGG_COL] = data_df.apply(lambda row: get_pathway_info(row), axis=1)

    data_df["KEGG_pathways"] = data_df["KEGG_pathways"].apply(
        lambda x: x["pathways"] if isinstance(x, dict) and "pathways" in x else []
    )

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
        "datasource": KEGG,
        "metadata": {"source_version": kegg_version},
        "query": {
            "size": len(gene_list),
            "input_type": KEGG_GENE_INPUT_ID,
            "number_of_added_edges": num_new_edges,
            "time": time_elapsed,
            "date": current_date,
            "url": KEGG_ENDPOINT,
        },
    }

    return data_df, kegg_metadata
