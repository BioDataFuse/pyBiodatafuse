#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Python file for querying IntAct (https://www.ebi.ac.uk/intact/)."""

import datetime
import json
import logging
import warnings
from time import sleep

import numpy as np
import pandas as pd
import requests

from pyBiodatafuse.constants import (
    INTACT,
    INTACT_COMPOUND_INPUT_ID,
    INTACT_COMPOUND_INTERACT_COL,
    INTACT_ENDPOINT,
    INTACT_GENE_INPUT_ID,
    INTACT_INTERACT_COL,
)
from pyBiodatafuse.utils import get_identifier_of_interest


def check_endpoint_intact() -> bool:
    """Check if the IntAct API is reachable by making a test request.

    :returns: True if the endpoint is available, False otherwise.
    """
    response = requests.get(f"{INTACT_ENDPOINT}/ws/interaction/findInteractions/P53")
    return response.status_code == 200


# def check_version_intact() -> dict:
#     """Get version of IntAct API.
#
#     :returns: a dictionary containing the version information
#     """
#     try:
#         version_call = requests.get(f"{INTACT_ENDPOINT}/version", timeout=10)
#         version_call.raise_for_status()
#         version_json = version_call.json()
#         return {"source_version": version_json.get("version", "unknown")}
#     except (requests.exceptions.RequestException, json.JSONDecodeError) as e:
#         logging.error("Error getting IntAct version")
#         return {"source_version": "unknown"}


def clean_id(identifier: str) -> str:
    """Strip the source suffix (e.g., ' (uniprotkb)') from an identifier string.

    :param identifier: The identifier string to clean.
    :returns: The cleaned identifier without the suffix.
    """
    return identifier.split(" ")[0] if identifier else identifier


def get_intact_interactions(gene_id: str):
    """Retrieve protein interactions for a given gene from IntAct.

    :param gene_id: Gene identifier
    :returns: List of interactions for the given gene
    """
    response = requests.get(f"{INTACT_ENDPOINT}/ws/interaction/findInteractions/{gene_id}")

    if response.status_code == 200:
        data = response.json()
        content = data.get("content", [])
        if not content:
            return []

        interactions = [
            {
                "interaction_id": item.get("ac", ""),
                "interactor_id_A": item.get("acA", ""),
                "interactor_id_B": item.get("acB", np.nan),
                "binary_interaction_id": item.get("binaryInteractionId", np.nan),
                "confidence_values": item.get("confidenceValues", []),
                "score": item.get("intactMiscore", []),
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
                "altIdsA": item.get("altIdsA", np.nan),
                "altIdsB": item.get("altIdsB", np.nan),
            }
            for item in data["content"]
        ]
        return interactions

    return []


def get_protein_intact_acs(id_of_interest: str) -> list:
    """Get all IntAct ACs for protein interactors from a given Ensembl ID.

    :param id_of_interest: input gene Ensembl identifier.
    :returns: Interactor information if possible, empty list if not.
    """
    url = f"{INTACT_ENDPOINT}/ws/interactor/findInteractor/{id_of_interest}?pageSize=100"
    try:
        response = requests.get(url, timeout=120)
        response.raise_for_status()
        data = response.json()

        content = data.get("content", [])

        protein_acs = []
        for item in content:
            interactor_type = item.get("interactorType")
            interactor_ac = item.get("interactorAc")

            if interactor_type == "protein":
                protein_acs.append(interactor_ac)

        return protein_acs

    except requests.exceptions.RequestException as e:
        logging.warning(f"Failed to get interactors for {id_of_interest}: {e}")
        return []


def get_filtered_interactions(
    query_id: str,
    valid_intact_acs: set,
    intact_ac_to_entity: dict,
    entity_to_input_id: dict,
    interaction_type: str = "gene_gene",
) -> list:
    """Filter interactions based on data type.

    :param query_id: The identifier to filter interactions for.
    :param valid_intact_acs: Set of valid IntAct ACs.
    :param intact_ac_to_entity: Dictionary mapping IntAct ACs to entity.
    :param entity_to_input_id: Dictionary mapping entities to input IDs.
    :param interaction_type: Either 'gene_gene', 'gene_compound' or 'both'.
    :returns: A list of filtered interactions.
    """
    interactions = get_intact_interactions(query_id)

    if len(interactions) == 1 and pd.isna(interactions[0].get("interaction_id", None)):
        return interactions

    filtered = []

    for interaction in interactions:
        id_a = interaction.get("interactor_id_A")
        id_b = interaction.get("interactor_id_B")
        alt_ids_a = interaction.get("altIdsA", []) or []
        alt_ids_b = interaction.get("altIdsB", []) or []

        has_uniprot_a = any("uniprotkb" in x.lower() for x in alt_ids_a)
        has_uniprot_b = any("uniprotkb" in x.lower() for x in alt_ids_b)
        has_chebi_a = any("chebi" in x.lower() for x in alt_ids_a)
        has_chebi_b = any("chebi" in x.lower() for x in alt_ids_b)

        keep_interaction = False

        if interaction_type == "gene_gene":
            if has_uniprot_a and has_uniprot_b:
                if id_a in valid_intact_acs and id_b in valid_intact_acs:
                    keep_interaction = True

        elif interaction_type == "gene_compound":
            if has_chebi_a or has_chebi_b:
                keep_interaction = True

        elif interaction_type == "compound_compound":
            if has_chebi_a and has_chebi_b:
                if id_a in valid_intact_acs and id_b in valid_intact_acs:
                    keep_interaction = True

        elif interaction_type == "compound_gene":
            if (has_chebi_a and has_uniprot_b) or (has_chebi_b and has_uniprot_a):
                keep_interaction = True

        elif interaction_type == "both":
            is_gene_gene = (
                has_uniprot_a
                and has_uniprot_b
                and id_a in valid_intact_acs
                and id_b in valid_intact_acs
            )
            is_gene_compound = (has_chebi_a and has_uniprot_b) or (has_chebi_b and has_uniprot_a)
            is_compound_compound = (
                has_chebi_a
                and has_chebi_b
                and id_a in valid_intact_acs
                and id_b in valid_intact_acs
            )
            if is_gene_gene or is_gene_compound or is_compound_compound:
                keep_interaction = True

        if not keep_interaction:
            continue

        if id_a in valid_intact_acs and intact_ac_to_entity.get(id_a) == query_id:
            partner_id = intact_ac_to_entity.get(id_b)
        else:
            partner_id = intact_ac_to_entity.get(id_a)

        partner_display_id = entity_to_input_id.get(partner_id, partner_id)
        interaction["intact_link_to"] = partner_display_id

        filtered.append(interaction)

    if not filtered:
        return [
            {
                "interaction_id": np.nan,
                "interactor_id_A": np.nan,
                "interactor_id_B": np.nan,
                "binary_interaction_id": np.nan,
                "confidence_values": [],
                "score": np.nan,
                "biological_role_A": np.nan,
                "biological_role_B": np.nan,
                "type": np.nan,
                "stoichiometry_A": np.nan,
                "stoichiometry_B": np.nan,
                "detection_method": np.nan,
                "detection_method_id": np.nan,
                "host_organism": np.nan,
                "interactor_A_name": np.nan,
                "interactor_B_name": np.nan,
                "interactor_A_species": np.nan,
                "interactor_B_species": np.nan,
                "molecule_A": np.nan,
                "molecule_B": np.nan,
                "id_A": np.nan,
                "id_B": np.nan,
                "pubmed_publication_id": np.nan,
                "ensembl": query_id,
                "altIdsA": np.nan,
                "altIdsB": np.nan,
            }
        ]

    return filtered


def get_gene_interactions(bridgedb_df: pd.DataFrame, interaction_type: str = "both"):
    """Annotate genes with interaction data from IntAct.

    :param bridgedb_df: BridgeDb output for creating the list of gene ids to query.
    :param interaction_type: Either 'gene_gene', 'gene_compound' or 'both'.
    :raises ValueError: If an invalid interaction_type is provided.
    :returns: a tuple (DataFrame containing the IntAct output, metadata dictionary)
    """
    api_available = check_endpoint_intact()
    if not api_available:
        warnings.warn("IntAct API endpoint is unavailable. Cannot retrieve data.", stacklevel=2)
        return pd.DataFrame(), {}

    start_time = datetime.datetime.now()
    data_df = get_identifier_of_interest(
        bridgedb_df,
        INTACT_GENE_INPUT_ID,
    ).reset_index(drop=True)

    if interaction_type not in ["gene_gene", "gene_compound", "both"]:
        raise ValueError(
            f"Invalid interaction_type: {interaction_type}. Must be 'gene' or 'compound'."
        )

    if isinstance(data_df, tuple):
        data_df = data_df[0]

    data_df = data_df.reset_index(drop=True)
    ensembl_gene_list = set(data_df["target"].tolist())

    ensembl_to_input_id = {}
    for _, row in data_df.iterrows():
        id_of_interest = row["target"]
        input_id = row["identifier"]
        ensembl_to_input_id[id_of_interest] = input_id

    ensembl_to_intact_map = {
        id_of_interest: get_protein_intact_acs(id_of_interest)
        for id_of_interest in ensembl_gene_list
    }

    intact_ac_to_ensembl = {
        ac: ensembl for ensembl, acs in ensembl_to_intact_map.items() for ac in acs
    }

    valid_intact_acs = {ac for acs in ensembl_to_intact_map.values() for ac in acs}

    data_df[INTACT_INTERACT_COL] = data_df["target"].apply(
        lambda gene_id: get_filtered_interactions(
            gene_id,
            valid_intact_acs,
            intact_ac_to_ensembl,
            ensembl_to_input_id,
            interaction_type=interaction_type,
        )
    )

    end_time = datetime.datetime.now()
    time_elapsed = str(end_time - start_time)
    current_date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    num_new_edges = sum(data_df[INTACT_INTERACT_COL].apply(len))

    intact_metadata = {
        "datasource": INTACT,
        "metadata": {"source_version": "unknown version"},
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


def get_compound_interactions(bridgedb_df: pd.DataFrame, interaction_type: str = "both"):
    """Annotate compounds with interaction data from IntAct.

    :param bridgedb_df: BridgeDb output for creating the list of compound ids to query.
    :param interaction_type: Either 'compound_compound', 'compound_gene' or 'both'.
    :raises ValueError: If an invalid interaction_type is provided.
    :returns: a tuple (DataFrame containing the IntAct output, metadata dictionary)
    """
    api_available = check_endpoint_intact()
    if not api_available:
        warnings.warn("IntAct API endpoint is unavailable. Cannot retrieve data.", stacklevel=2)
        return pd.DataFrame(), {}

    start_time = datetime.datetime.now()
    data_df = get_identifier_of_interest(
        bridgedb_df,
        INTACT_COMPOUND_INPUT_ID,
    ).reset_index(drop=True)

    data_df = data_df[data_df["target"].str.startswith("CHEBI:")].reset_index(drop=True)

    if interaction_type not in ["compound_compound", "compound_gene", "both"]:
        raise ValueError(
            f"Invalid interaction_type: {interaction_type}. Must be 'compound_compound', 'compound_gene' or 'both'."
        )

    if isinstance(data_df, tuple):
        data_df = data_df[0]

    data_df = data_df.reset_index(drop=True)
    chebi_list = set(data_df["target"].tolist())

    chebi_to_input_id = {row["target"]: row["identifier"] for _, row in data_df.iterrows()}

    intact_ac_to_chebi = {chebi_id: chebi_id for chebi_id in chebi_list}

    valid_intact_acs = chebi_list

    data_df[INTACT_COMPOUND_INTERACT_COL] = data_df["target"].apply(
        lambda compound_id: get_filtered_interactions(
            compound_id,
            valid_intact_acs,
            intact_ac_to_chebi,
            chebi_to_input_id,
            interaction_type=interaction_type,
        )
    )

    end_time = datetime.datetime.now()
    time_elapsed = str(end_time - start_time)
    current_date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    num_new_edges = sum(data_df[INTACT_COMPOUND_INTERACT_COL].apply(len))

    intact_metadata = {
        "datasource": INTACT,
        "metadata": {"source_version": "unknown version"},
        "query": {
            "size": len(chebi_list),
            "input_type": INTACT_COMPOUND_INPUT_ID,
            "number_of_added_edges": num_new_edges,
            "time": time_elapsed,
            "date": current_date,
            "url": INTACT_ENDPOINT,
        },
    }

    return data_df, intact_metadata
