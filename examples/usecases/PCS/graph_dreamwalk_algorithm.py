#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Codes to run Dreamwalk algorithm on BDF graph."""

import logging
import os

import pandas as pd

import pyBiodatafuse.algorithms.DREAMwalk.generate_files as gen
from pyBiodatafuse.algorithms.DREAMwalk.calculate_drug_scores import find_candidates
from pyBiodatafuse.algorithms.DREAMwalk.constant import (
    DISEASE_SIM_FILE,
    DRUG_HIERACHY_FILE,
    DRUG_SIM_FILE,
    EMBEDDING_FILE,
    GRAPH_FILE,
    NODE_TYPE_FILE,
    SIM_FILE,
)
from pyBiodatafuse.algorithms.DREAMwalk.generate_embeddings import save_embedding_files
from pyBiodatafuse.algorithms.DREAMwalk.generate_similarity_net import save_sim_graph
from pyBiodatafuse.algorithms.DREAMwalk.predict_associations import predict_dda
from pyBiodatafuse.analyzer.summarize import BioGraph

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def preprocess_files() -> None:
    """Generating base files for the algorithm."""
    graph_obj = BioGraph(graph_path="PCS_graph.gml", graph_format="gml")
    gen.create_files(graph_obj=graph_obj, output_dir="dreamwalk_data")
    logger.info("Files generated successfully.")


def generate_sim(
    drug_hierarchy_file: str,
    dis_sim_file: str,
    drug_sim_file: str,
    network_file: str,
    output_file: str,
) -> None:
    """Generate similarity graph.

    :param drug_hierarchy_file: Drug hierarchy file path
    :param dis_sim_file: Disease similarity file path
    :param drug_sim_file: Drug similarity file path
    :param network_file: Network file path
    :param output_file: Final similarity file path
    """
    if not os.path.exists(drug_sim_file):
        save_sim_graph(
            networkf=network_file, hierf=drug_hierarchy_file, outputf=drug_sim_file, cutoff=0.4
        )
        logger.info("Similarity graph generated successfully.")

    # create similarity graph merging disease and drugs
    drug_sim = pd.read_csv(drug_sim_file, sep="\t")
    dis_sim = pd.read_csv(dis_sim_file, sep="\t")

    similarty_graph = pd.concat([dis_sim, drug_sim], ignore_index=True)
    similarty_graph["edge_id"] = range(0, len(similarty_graph))  # reindexing the counts
    similarty_graph.to_csv(output_file, sep="\t", index=False)


def generate_node_embeddings(
    similarity_file: str,
    node_type_file: str,
    network_file: str,
    output_file: str,
) -> None:
    """Generate node embeddings."""

    save_embedding_files(
        netf=network_file,
        sim_netf=similarity_file,
        outputf=output_file,
        nodetypef=node_type_file,
        tp_factor=0.3,
    )
    logger.info("Node embeddings generated successfully.")


def train_model(embedding_file: str, pair_file: str, model_file: str) -> None:
    """Predict drug-disease association.

    :param embedding_file: Embedding file path
    :param pair_file: Pair file path
    :param model_file: Model file path
    """
    predict_dda(embeddingf=embedding_file, pairf=pair_file, modelf=model_file)
    logger.info("Drug-disease association predicted successfully.")


def find_potential_drugs(
    disease_id: str,
    output_dir: str,
    embeddingf: str,
    model_folder: str,
    network_file: str,
    node_type_file: str,
) -> None:
    """Calculate drug scores.

    :param disease_id: Disease ID for which drug scores are to be calculated
    :param output_dir: Output directory
    :param embeddingf: Embedding file path
    :param model_folder: Folder where all models are stored
    :param network_file: Network file path
    :param node_type_file: Node type file path
    :return: None
    """
    all_results = find_candidates(
        networkf=network_file,
        nodetypef=node_type_file,
        embeddingf=embeddingf,
        model_folder=model_folder,
        query_disease=disease_id,
    )  # By deault all the drugs are considered; to limit the drugs, use candidate_count parameter
    all_results.to_csv(f"{output_dir}/model_predictions.tsv", sep="\t", index=False)
    logger.info("Drug scores calculated successfully.")


if __name__ == "__main__":
    networkf = f"dreamwalk_data/{GRAPH_FILE}"
    hierarchyf = f"dreamwalk_data/{DRUG_HIERACHY_FILE}"
    drug_simf = f"dreamwalk_data/{DRUG_SIM_FILE}"
    disease_simf = f"dreamwalk_data/{DISEASE_SIM_FILE}"
    nodetypef = f"dreamwalk_data/{NODE_TYPE_FILE}"
    embeddingf = f"dreamwalk_data/{EMBEDDING_FILE}"
    simf = f"dreamwalk_data/{SIM_FILE}"

    preprocess_files()
    generate_sim(
        drug_hierarchy_file=hierarchyf,
        dis_sim_file=disease_simf,
        drug_sim_file=drug_simf,
        network_file=networkf,
        output_file=simf,
    )
    generate_node_embeddings(
        similarity_file=simf,
        node_type_file=nodetypef,
        network_file=networkf,
        output_file=embeddingf,
    )

    MODEL_DIR = "dreamwalk_data/dda"

    for i in range(1, 6):
        pairf = f"{MODEL_DIR}/dda_{i}.tsv"
        modelf = f"{MODEL_DIR}/clf_{i}.pkl"
        train_model(embedding_file=embeddingf, pair_file=pairf, model_file=modelf)

    find_potential_drugs(
        disease_id="UMLS:C5433293",  # Long COVID node id
        output_dir="dreamwalk_data",
        embeddingf=embeddingf,
        model_folder=MODEL_DIR,
        network_file=networkf,
        node_type_file=nodetypef,
    )
