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
    response = requests.get(f"{INTACT_ENDPOINT}/ws/interaction/findInteractions/P53")
    return response.status_code == 200


def check_version_intact() -> dict:
    """Get version of IntAct API.

    :returns: a dictionary containing the version information
    """
    try:
        version_call = requests.get(f"{INTACT_ENDPOINT}/version", timeout=10)
        version_call.raise_for_status()
        version_json = version_call.json()
        return {"source_version": version_json.get("version", "unknown")}
    except (requests.exceptions.RequestException, json.JSONDecodeError) as e:
        logging.error("Error getting IntAct version")
        return {"source_version": "unknown"}


def get_compound_related_interactions():
    # Use 'CHEBI:' prefix to search for compound interactions
    response = requests.get(f"{INTACT_ENDPOINT}/ws/interaction/search/CHEBI", timeout=10)
    if response.status_code != 200:
        print("Failed to query IntAct.")
        return []

    data = response.json()
    return data.get("content", [])


def clean_id(identifier: str) -> str:
    """Strip the source suffix (e.g., ' (uniprotkb)') from an identifier string."""
    if identifier and isinstance(identifier, str):
        return identifier.split(" ")[0]
    return identifier


def get_intact_interactions(gene_id: str):
    """Retrieve protein interactions for a given gene from IntAct.

    :param gene_id: Gene identifier (Ensembl ID)
    :returns: List of interactions for the given gene
    """
    response = requests.get(f"{INTACT_ENDPOINT}/ws/interaction/findInteractions/{gene_id}")

    if response.status_code == 200:
        data = response.json()
        if not data.get("content"):
            return []

        interactions = [
            {
                "interaction_id": item.get("ac", ""),
                "interactor_id_A": item.get("acA", ""),
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
                "pubmed_publication_id": item.get("publicationPubmedIdentifier", np.nan),
                "ensembl": gene_id,
            }
            for item in data["content"]
        ]
        return interactions

    return []


def get_protein_intact_acs(ensembl_id: str) -> list:
    """Get all IntAct ACs for protein interactors from a given Ensembl ID."""
    url = f"{INTACT_ENDPOINT}/ws/interactor/findInteractor/{ensembl_id}"
    try:
        response = requests.get(url, timeout=10)
        response.raise_for_status()
        data = response.json()
        protein_acs = [
            item["interactorAc"]
            for item in data.get("content", [])
            if item.get("interactorType") == "protein"
        ]
        return protein_acs
    except requests.exceptions.RequestException as e:
        logging.warning(f"Failed to get interactors for {ensembl_id}: {e}")
        return []


def get_filtered_interactions(gene_id: str, valid_intact_acs: set, intact_ac_to_ensembl: dict) -> list:
    """Get IntAct interactions for a gene, filtered to include only protein-protein interactions between input genes.
       Adds the Ensembl ID of the partner gene to each interaction.
    """
    interactions = get_intact_interactions(gene_id)
    filtered = []

    for interaction in interactions:
        id_a = interaction.get("interactor_id_A")
        id_b = interaction.get("interactor_id_B")

        if id_a in valid_intact_acs and id_b in valid_intact_acs:
            if intact_ac_to_ensembl.get(id_a) == gene_id:
                partner_gene = intact_ac_to_ensembl.get(id_b, None)
            else:
                partner_gene = intact_ac_to_ensembl.get(id_a, None)

            interaction["intact_link_to"] = partner_gene
            filtered.append(interaction)

    return filtered


def get_compound_filtered_interactions(gene_id: str) -> list:
    """
    Get IntAct interactions for a gene that involve chemical compounds (e.g., ChEBI).

    Returns interactions where at least one interactor has a ChEBI ID.
    """
    interactions = get_intact_interactions(gene_id)
    filtered = []

    for interaction in interactions:
        id_a = interaction.get("id_A", "")
        id_b = interaction.get("id_B", "")

        if isinstance(id_a, str) and "CHEBI:" in id_a or isinstance(id_b, str) and "CHEBI:" in id_b:
            filtered.append(interaction)

    return filtered


def get_interactions(bridgedb_df: pd.DataFrame):
    """Annotate genes with interaction data from IntAct, filtering with UniProt IDs."""

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

    uniprot_map = {}
    for _, row in uniprot_data_df.iterrows():
        ensembl_id = row["identifier"]
        uniprot_id = row["target"]
        if ensembl_id not in uniprot_map:
            uniprot_map[ensembl_id] = []
        uniprot_map[ensembl_id].append(uniprot_id)

    ensembl_to_intact_map = {}
    for ensembl_id in ensembl_gene_list:
        intact_acs = get_protein_intact_acs(ensembl_id)
        ensembl_to_intact_map[ensembl_id] = intact_acs
        print(ensembl_id, ensembl_gene_list)

    intact_ac_to_ensembl = {
        ac: ensembl for ensembl, acs in ensembl_to_intact_map.items() for ac in acs
    }

    valid_intact_acs = {ac for acs in ensembl_to_intact_map.values() for ac in acs}

    # Retrieve interactions from IntAct using input IDs
    intact_interactions = []
    for ensembl_id in ensembl_gene_list:
        intact_interactions.extend(get_intact_interactions(ensembl_id))

    data_df["IntAct_interactions"] = data_df["target"].apply(
    lambda gene_id: get_filtered_interactions(gene_id, valid_intact_acs, intact_ac_to_ensembl)
    )

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


def get_compound_interactions(bridgedb_df: pd.DataFrame):
    """Annotate genes with compound-related interaction data from IntAct."""

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

    # Fetch compound interactions for each gene
    data_df["IntAct_compound_interactions"] = data_df["target"].apply(
        lambda gene_id: get_compound_filtered_interactions(gene_id)
    )

    end_time = datetime.datetime.now()
    time_elapsed = str(end_time - start_time)
    current_date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    num_new_edges = sum(data_df["IntAct_compound_interactions"].apply(len))

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