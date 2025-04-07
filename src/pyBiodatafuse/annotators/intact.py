#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Python file for querying IntAct (https://www.ebi.ac.uk/intact/).
"""

import datetime
import json
import logging
import numpy as np
import pandas as pd
import requests
import warnings
from time import sleep

from pyBiodatafuse.constants import (
    INTACT,
    INTACT_GENE_INPUT_ID,
    INTACT_ENDPOINT,
)
from pyBiodatafuse.utils import get_identifier_of_interest

def check_endpoint_intact() -> bool:
    """Check if the IntAct API is reachable by making a test request.

    :returns: True if the endpoint is available, else False
    """
    response = requests.get(f"{INTACT_ENDPOINT}/findInteractions/P53")
    return response.status_code == 200

def check_version_intact() -> str:
    """Check the current version of the IntAct database.
    
    :returns: A string containing the version information
    """
    response = requests.get(f"{INTACT_ENDPOINT}/version")
    if response.status_code == 200:
        return response.json().get("version", "Unknown")
    return "Error: Version not found."

def clean_id(identifier: str) -> str:
    """Strip the source suffix (e.g., ' (uniprotkb)') from an identifier string."""
    if identifier and isinstance(identifier, str):
        return identifier.split(" ")[0]
    return identifier

def get_intact_interactions(gene_id: str):
    """Retrieve protein interactions for a given gene from IntAct.

    :param gene_id: Gene identifier
    :returns: List of interactions for the given gene
    """
    response = requests.get(f"{INTACT_ENDPOINT}/findInteractions/{gene_id}")

    if response.status_code == 200:
        data = response.json()
        if not data.get("content"):
            return []

        interactions = [
            {
                "interaction_id": item.get("ac", np.nan),
                "interactor_id_A": item.get("acA", np.nan),
                "interactor_id_B": item.get("acB", np.nan),
                "binary_interaction_id": item.get("binaryInteractionId", np.nan),
                "confidence_values": item.get("confidenceValues", []),
                "intact_score": item.get("intactMiscore", []),
                "biological_role_A": item.get("biologicalRoleA", np.nan),
                "biological_role_B": item.get("biologicalRoleB", np.nan),
                "type": item.get("type", np.nan),
                "stoichiometry_A": item.get("stoichiometryA", np.nan),
                "stoichiometry_B": item.get("stoichiometryB", np.nan),
                "detection_method": item.get("detectionMethod", np.nan),
                "detection_method_id": item.get("detectionMethodMIIdentifier", np.nan),
                "host_organism": item.get("hostOrganism", np.nan),
                "interactor_A_name": item.get("intactNameA", np.nan),
                "interactor_B_name": item.get("intactNameB", np.nan),
                "interactor_A_species": item.get("speciesA", np.nan),
                "interactor_B_species": item.get("speciesB", np.nan),
                "molecule_A": item.get("moleculeA", np.nan),
                "molecule_B": item.get("moleculeB", np.nan),
                "id_A": clean_id(item.get("idA", np.nan)),
                "id_B": clean_id(item.get("idB", np.nan)),
                "pubmed_publication_id": item.get("publicationPubmedIdentifier", []),
            }
            for item in data["content"]
        ]
        return interactions

    return []

def get_filtered_interactions(gene_id: str, uniprot_map: dict) -> list:
    """Retrieve IntAct interactions for a gene, filtered to include only interactions with UniProt IDs in uniprot_map."""
    cleaned_gene_id = clean_id(gene_id)
    interactions = get_intact_interactions(gene_id)
    filtered = []

    
    for interaction in interactions:
        id_a = interaction.get("id_A")
        id_b = interaction.get("id_B")

        is_valid_interaction = False
        for ensembl_id, uniprot_ids in uniprot_map.items():
            if (id_a in uniprot_ids) or (id_b in uniprot_ids):
                is_valid_interaction = True
                break
        
        if is_valid_interaction:
            filtered.append(interaction)
    
    return filtered


def get_interactions(bridgedb_df: pd.DataFrame):
    """Annotate genes with interaction data from IntAct, filtering with UniProt IDs.

    :param bridgedb_df: BridgeDb output for creating the list of gene ids to query
    :returns: a tuple (DataFrame containing the filtered IntAct output, metadata dictionary)
    """
    api_available = check_endpoint_intact()
    if not api_available:
        warnings.warn("IntAct API endpoint is unavailable. Cannot retrieve data.", stacklevel=2)
        return pd.DataFrame(), {}

    intact_version = check_version_intact()
    start_time = datetime.datetime.now()

    # Get identifiers of interest
    data_df = get_identifier_of_interest(bridgedb_df, INTACT_GENE_INPUT_ID)
    uniprot_data_df = get_identifier_of_interest(bridgedb_df, "Uniprot-TrEMBL")

    if isinstance(data_df, tuple):
        data_df = data_df[0]

    data_df = data_df.reset_index(drop=True)
    ensembl_gene_list = set(data_df["target"].tolist())

    # Build a map of UniProt IDs from the uniprot_data_df (mapping Ensembl ID to list of UniProt IDs)
    uniprot_map = {}
    for _, row in uniprot_data_df.iterrows():
        ensembl_id = row["identifier"]
        uniprot_id = row["target"]
        if ensembl_id not in uniprot_map:
            uniprot_map[ensembl_id] = []
        uniprot_map[ensembl_id].append(uniprot_id)

    # Retrieve interactions from IntAct using Ensembl IDs
    intact_interactions = []
    for ensembl_id in ensembl_gene_list:
        intact_interactions.extend(get_intact_interactions(ensembl_id))

    # Filter interactions based on UniProt IDs
    filtered_interactions = []
    for interaction in intact_interactions:
        id_A = interaction.get("id_A")
        id_B = interaction.get("id_B")

        # Filter the interaction based on UniProt mapping
        if id_A in uniprot_map and id_B in uniprot_map:
            filtered_interactions.append(interaction)

    # Adding filtered interactions to the DataFrame
    data_df["IntAct_interactions"] = data_df["target"].apply(lambda gene_id: get_filtered_interactions(gene_id, uniprot_map))

    end_time = datetime.datetime.now()
    time_elapsed = str(end_time - start_time)
    current_date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    num_new_edges = sum(data_df["IntAct_interactions"].apply(len))

    intact_metadata = {
        "datasource": INTACT,
        "metadata": {"source_version": intact_version},
        "query": {
            "size": len(ensembl_gene_list),
            "input_type": INTACT_GENE_INPUT_ID,
            "number_of_added_edges": num_new_edges,
            "time": time_elapsed,
            "date": current_date,
            "url": INTACT_ENDPOINT,
        },
    }

    return data_df, intact_metadata
