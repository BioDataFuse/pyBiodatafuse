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
import datetime
from time import sleep

from pyBiodatafuse.constants import (
    INTACT,
    INTACT_GENE_INPUT_ID,
    INTACT_ENDPOINT,
)
from pyBiodatafuse.utils import check_columns_against_constants, get_identifier_of_interest


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

def get_filtered_interactions(gene_id: str, gene_list: set) -> list:
    """Retrieve IntAct interactions for a gene, filtered to include only interactions with genes in gene_list."""
    cleaned_gene_id = clean_id(gene_id)  # Clean the input for comparison
    interactions = get_intact_interactions(gene_id)
    filtered = []
    for interaction in interactions:
        id_a = interaction.get("id_A")
        id_b = interaction.get("id_B")
        partner_id = id_b if id_a == cleaned_gene_id else id_a
        if partner_id in gene_list:
            filtered.append(interaction)
    return filtered



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


def get_interactions(bridgedb_df: pd.DataFrame):
    """Annotate genes with interaction data from IntAct, filtering with Ensembl IDs and mapping to UniProt IDs.

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

    if isinstance(data_df, tuple):
        data_df = data_df[0]

    data_df = data_df.reset_index(drop=True)
    ensembl_gene_list = set(data_df["target"].tolist())

    # Build a list of Ensembl IDs
    ensembl_ids = list(ensembl_gene_list)

    # Retrieve interactions from IntAct using Ensembl IDs
    intact_interactions = get_intact_interactions(ensembl_ids)

    # Map Ensembl IDs to UniProt IDs
    uniprot_map = ensp_to_uniprot(ensembl_ids)

    # Process and filter the interactions
    filtered_interactions = []
    for interaction in intact_interactions:
        # Map Ensembl IDs to UniProt for both interactors
        id_A = interaction.get("id_A")
        id_B = interaction.get("id_B")

        # Filter the interaction based on UniProt mapping
        if id_A in uniprot_map and id_B in uniprot_map:
            filtered_interactions.append(interaction)

    # Adding UniProt info to the DataFrame
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

def ensp_to_uniprot(ensp_ids):
    """Retrieve UniProt IDs from Ensembl protein IDs (ENSP).

    :param ensp_ids: List of Ensembl protein IDs (ENSP)
    :return: Dictionary mapping ENSP IDs to UniProt IDs
    """
    ensp_to_uniprot_map = {}
    url = "https://rest.uniprot.org/idmapping/run"
    headers = {"Content-Type": "application/x-www-form-urlencoded"}
    data = {"from": "Ensembl_Protein", "to": "UniProtKB", "ids": ",".join(ensp_ids)}

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