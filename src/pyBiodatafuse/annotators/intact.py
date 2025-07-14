#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Python file for querying IntAct (https://www.ebi.ac.uk/intact/)."""

import datetime
import json
import logging
import urllib.parse
import warnings
from typing import Dict, List

import numpy as np
import pandas as pd
import requests
from tqdm import tqdm

import pyBiodatafuse.constants as Cons
from pyBiodatafuse.utils import get_identifier_of_interest, give_annotator_warning


def check_endpoint_intact() -> bool:
    """Check if the IntAct API is reachable by making a test request.

    :returns: True if the endpoint is available, False otherwise.
    """
    response = requests.get(f"{Cons.INTACT_ENDPOINT}/ws/interaction/findInteractions/P53")
    return response.status_code == 200


# TODO: Wait for this function to be implemented in the IntAct API
def check_version_intact() -> dict:
    """Get version of IntAct API.

    :returns: a dictionary containing the version information
    """
    try:
        version_call = requests.get(f"{Cons.INTACT_ENDPOINT}/version", timeout=10)
        version_call.raise_for_status()
        version_json = version_call.json()
        return {"source_version": version_json.get("version", "unknown")}
    except (requests.exceptions.RequestException, json.JSONDecodeError) as e:
        logging.error(f"Error getting IntAct version: {e}")
        return {"source_version": "unknown"}


def clean_id(identifier: str) -> str:
    """Strip the source suffix (e.g., ' (uniprotkb)') from an identifier string.

    :param identifier: The identifier string to clean.
    :returns: The cleaned identifier without the suffix.
    """
    return identifier.split(" ")[0] if identifier else identifier


def get_intact_interactions(gene_ids: List[str]) -> List[dict]:
    """Retrieve protein interactions for a list of genes from IntAct.

    :param gene_ids: List of gene identifiers.
    :returns: List of interactions for the given genes.
    """
    if not gene_ids:
        return []

    joined_ids = " - ".join(gene_ids)
    encoded_ids = urllib.parse.quote(joined_ids)
    url = f"{Cons.INTACT_ENDPOINT}/ws/interaction/findInteractions/{encoded_ids}?pageSize=200"

    try:
        response = requests.get(url, timeout=60)
        data = response.json()

        content = data.get("content", [])
        if not content:
            return []

        interation_info = {
            Cons.INTACT_INTERACTION_ID: "ac",
            Cons.INTACT_INTERACTOR_ID_A: "acA",
            Cons.INTACT_INTERACTOR_ID_B: "acB",
            Cons.INTACT_SCORE: "intactMiscore",
            Cons.INTACT_BIOLOGICAL_ROLE_A: "biologicalRoleA",
            Cons.INTACT_BIOLOGICAL_ROLE_B: "biologicalRoleB",
            Cons.INTACT_TYPE: "type",
            Cons.INTACT_DETECTION_METHOD: "detectionMethod",
            Cons.INTACT_HOST_ORGANISM: "hostOrganism",
            Cons.INTACT_INTERACTOR_A_NAME: "intactNameA",
            Cons.INTACT_INTERACTOR_B_NAME: "intactNameB",
            Cons.INTACT_INTERACTOR_A_SPECIES: "speciesA",
            Cons.INTACT_INTERACTOR_B_SPECIES: "speciesB",
            Cons.INTACT_MOLECULE_A: "moleculeA",
            Cons.INTACT_MOLECULE_B: "moleculeB",
            Cons.INTACT_ID_A: "idA",
            Cons.INTACT_ID_B: "idB",
            Cons.INTACT_PUBMED_PUBLICATION_ID: "publicationPubmedIdentifier",
        }

        interactions = [
            {key: item.get(value, np.nan) for key, value in interation_info.items()}
            for item in content
        ]

        # cleanup the alternative ids
        for interaction in interactions:
            ids_a = interaction[Cons.INTACT_ID_A]
            ids_b = interaction[Cons.INTACT_ID_B]

            if ":" in ids_a:
                interaction[Cons.INTACT_ID_A] = ids_a.split(" ")[0]  # stays the same
            else:
                idx = ids_a.split(" ")[0]
                namespace = ids_a.split(" ")[1].replace("(", "").replace(")", "")
                interaction[Cons.INTACT_ID_A] = f"{namespace}:{idx}"

            if ":" in ids_b:
                interaction[Cons.INTACT_ID_B] = ids_b.split(" ")[0]  # stays the same
            else:
                idx = ids_b.split(" ")[0]
                namespace = ids_b.split(" ")[1].replace("(", "").replace(")", "")
                interaction[Cons.INTACT_ID_B] = f"{namespace}:{idx}"

        return interactions

    except requests.RequestException as e:
        logging.warning(f"Batch request failed for genes {gene_ids}: {e}")
        return []


def get_protein_intact_acs(id_of_interest: str) -> List[str]:
    """Get all IntAct ACs for protein interactors from a given Ensembl ID.

    :param id_of_interest: input gene Ensembl identifier.
    :returns: Interactor information if possible, empty list if not.
    """
    url = f"{Cons.INTACT_ENDPOINT}/ws/interactor/findInteractor/{id_of_interest}?pageSize=100"
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
    batch_ids: List[str],
    valid_intact_acs: set,
    intact_ac_to_entity: dict,
    entity_to_input_id: dict,
    interaction_type: str = "gene_gene",
) -> Dict[str, List[dict]]:
    """Filter interactions based on data type.

    :param batch_ids: List of input IDs.
    :param valid_intact_acs: Set of valid IntAct ACs.
    :param intact_ac_to_entity: Dictionary mapping IntAct ACs to entity.
    :param entity_to_input_id: Dictionary mapping entities to input IDs.
    :param interaction_type: Either 'gene_gene', 'gene_compound', 'compound_compound', 'compound_gene', or 'both'.
    :returns: A dictionary of filtered interactions per input ID.
    """
    results: Dict[str, List[dict]] = {idx: [] for idx in batch_ids}
    interactions = get_intact_interactions(batch_ids)

    for interaction in interactions:
        if interaction_type in Cons.INTACT_GENE_INTERACTION_TYPES:
            id_a = interaction.get(Cons.INTACT_INTERACTOR_ID_A)
            id_b = interaction.get(Cons.INTACT_INTERACTOR_ID_B)
            alt_ids_a = interaction.get(Cons.INTACT_ID_A)
            alt_ids_b = interaction.get(Cons.INTACT_ID_B)
        else:
            id_a = interaction.get(Cons.INTACT_ID_A)
            id_b = interaction.get(Cons.INTACT_ID_B)
            alt_ids_a = interaction.get(Cons.INTACT_INTERACTOR_ID_A)
            alt_ids_b = interaction.get(Cons.INTACT_INTERACTOR_ID_B)

        has_uniprot_a = any("uniprotkb" in x.lower() for x in [id_a, alt_ids_a])  # type: ignore
        has_uniprot_b = any("uniprotkb" in x.lower() for x in [id_b, alt_ids_b])  # type: ignore
        has_chebi_a = any("chebi" in x.lower() for x in [id_a, alt_ids_a])  # type: ignore
        has_chebi_b = any("chebi" in x.lower() for x in [id_b, alt_ids_b])  # type: ignore

        keep_interaction = False

        if interaction_type == "gene_gene":
            if (
                has_uniprot_a
                and has_uniprot_b
                and id_a in valid_intact_acs
                and id_b in valid_intact_acs
            ):
                keep_interaction = True

        elif interaction_type == "gene_compound":
            if (has_chebi_a and has_uniprot_b) or (has_chebi_b and has_uniprot_a):
                keep_interaction = True

        elif interaction_type == "compound_compound":
            if (
                has_chebi_a
                and has_chebi_b
                and id_a in valid_intact_acs
                and id_b in valid_intact_acs
            ):
                keep_interaction = True

        elif interaction_type == "compound_gene":
            if (has_chebi_a and has_uniprot_b) or (has_chebi_b and has_uniprot_a):
                keep_interaction = True

        elif "both" in interaction_type:
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

        for idx in batch_ids:
            if id_a in valid_intact_acs and intact_ac_to_entity.get(id_a) == idx:
                partner_id = intact_ac_to_entity.get(id_b)
            elif id_b in valid_intact_acs and intact_ac_to_entity.get(id_b) == idx:
                partner_id = intact_ac_to_entity.get(id_a)
            else:
                continue

            partner_display_id = entity_to_input_id.get(partner_id, partner_id)
            interaction_copy = dict(interaction)
            interaction_copy["intact_link_to"] = partner_display_id
            results[idx].append(interaction_copy)

    for gene_id in batch_ids:
        if not results[gene_id]:
            results[gene_id] = [{key: np.nan for key in Cons.INTACT_OUTPUT_DICT}]

    return results


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
    data_df = get_identifier_of_interest(bridgedb_df, Cons.INTACT_GENE_INPUT_ID).reset_index(
        drop=True
    )

    if interaction_type not in Cons.INTACT_GENE_INTERACTION_TYPES:
        raise ValueError(
            f"Invalid interaction_type: {interaction_type}. Must be {Cons.INTACT_GENE_INTERACTION_TYPES}."
        )

    ensembl_gene_list = list(set(data_df[Cons.TARGET_COL].tolist()))

    ensembl_to_input_id = {
        row[Cons.TARGET_COL]: row[Cons.IDENTIFIER_COL] for _, row in data_df.iterrows()
    }

    ensembl_to_intact_map = {
        gene_id: get_protein_intact_acs(gene_id) for gene_id in ensembl_gene_list
    }

    intact_ac_to_ensembl = {ac: gene for gene, acs in ensembl_to_intact_map.items() for ac in acs}

    all_results = {}
    batch_size = 10
    for i in tqdm(range(0, len(ensembl_gene_list), batch_size), desc="Querying IntAct for genes"):
        batch = ensembl_gene_list[i : i + batch_size]
        batch_results = get_filtered_interactions(
            batch,
            set(intact_ac_to_ensembl.keys()),
            intact_ac_to_ensembl,
            ensembl_to_input_id,
            interaction_type=interaction_type,
        )
        all_results.update(batch_results)

    data_df[Cons.INTACT_INTERACT_COL] = data_df[Cons.TARGET_COL].map(all_results)

    end_time = datetime.datetime.now()
    time_elapsed = str(end_time - start_time)
    current_date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    intact_metadata = {
        "datasource": Cons.INTACT,
        "metadata": {"source_version": "unknown"},
        "query": {
            "size": len(ensembl_gene_list),
            "input_type": Cons.INTACT_GENE_INPUT_ID,
            "time": time_elapsed,
            "date": current_date,
            "url": Cons.INTACT_ENDPOINT,
        },
    }

    """Update metadata"""
    # Calculate the number of new nodes
    num_new_nodes = 0  # TODO: Implement this

    # Calculate the number of new edges
    num_new_edges = len(all_results)

    # Check the intermediate_df
    if num_new_edges != len(data_df):
        give_annotator_warning(Cons.INTACT)

    # Add the number of new nodes and edges to metadata
    intact_metadata[Cons.QUERY][Cons.NUM_NODES] = num_new_nodes  # type: ignore
    intact_metadata[Cons.QUERY][Cons.NUM_EDGES] = num_new_edges  # type: ignore

    return data_df, intact_metadata


def get_compound_interactions(bridgedb_df: pd.DataFrame, interaction_type: str = "both_compounds"):
    """Annotate compounds with interaction data from IntAct.

    :param bridgedb_df: BridgeDb output for creating the list of compound ids to query.
    :param interaction_type: Either 'compound_compound', 'compound_gene' or 'both_compounds'.
    :raises ValueError: If an invalid interaction_type is provided.
    :returns: a tuple (DataFrame containing the IntAct output, metadata dictionary)
    """
    api_available = check_endpoint_intact()
    if not api_available:
        warnings.warn("IntAct API endpoint is unavailable. Cannot retrieve data.", stacklevel=2)
        return pd.DataFrame(), {}

    start_time = datetime.datetime.now()
    data_df = get_identifier_of_interest(bridgedb_df, Cons.INTACT_COMPOUND_INPUT_ID).reset_index(
        drop=True
    )
    data_df = data_df[data_df[Cons.TARGET_COL].str.startswith("CHEBI:")].reset_index(drop=True)

    if interaction_type not in Cons.INTACT_COMPOUND_INTERACTION_TYPES:
        raise ValueError(
            f"Invalid interaction_type: {interaction_type}. Must be {Cons.INTACT_COMPOUND_INTERACTION_TYPES}."
        )

    chebi_list = list(set(data_df[Cons.TARGET_COL].tolist()))

    chebi_to_input_id = {
        row[Cons.TARGET_COL]: row[Cons.IDENTIFIER_COL] for _, row in data_df.iterrows()
    }

    intact_ac_to_chebi = {
        chebi_id: chebi_id for chebi_id in chebi_list
    }  # intact id is same as chebi id

    all_results = {}
    batch_size = 10
    for i in tqdm(range(0, len(chebi_list), batch_size), desc="Querying IntAct for compounds"):
        batch = chebi_list[i : i + batch_size]
        batch_results = get_filtered_interactions(
            batch_ids=batch,
            valid_intact_acs=set(chebi_list),
            intact_ac_to_entity=intact_ac_to_chebi,
            entity_to_input_id=chebi_to_input_id,
            interaction_type=interaction_type,
        )
        all_results.update(batch_results)

    data_df[Cons.INTACT_COMPOUND_INTERACT_COL] = data_df[Cons.TARGET_COL].map(all_results)

    end_time = datetime.datetime.now()
    time_elapsed = str(end_time - start_time)
    current_date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    intact_metadata = {
        "datasource": Cons.INTACT,
        "metadata": {"source_version": "unknown"},
        "query": {
            "size": len(chebi_list),
            "input_type": Cons.INTACT_COMPOUND_INPUT_ID,
            "time": time_elapsed,
            "date": current_date,
            "url": Cons.INTACT_ENDPOINT,
        },
    }

    """Update metadata"""
    # Calculate the number of new nodes
    num_new_nodes = 0  # TODO: Implement this

    # Calculate the number of new edges
    num_new_edges = sum(data_df[Cons.INTACT_COMPOUND_INTERACT_COL].apply(len))

    # Check the intermediate_df
    if num_new_edges != len(data_df):
        give_annotator_warning(Cons.INTACT)

    # Add the number of new nodes and edges to metadata
    intact_metadata[Cons.QUERY][Cons.NUM_NODES] = num_new_nodes  # type: ignore
    intact_metadata[Cons.QUERY][Cons.NUM_EDGES] = num_new_edges  # type: ignore

    return data_df, intact_metadata
