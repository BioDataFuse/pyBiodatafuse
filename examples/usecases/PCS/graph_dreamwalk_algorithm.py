#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Codes to run Dreamwalk algorithm on BDF graph."""

import logging
import os

import pandas as pd

import pyBiodatafuse.algorithms.DREAMwalk.generate_files as gen
from pyBiodatafuse.algorithms.DREAMwalk.calculate_drug_scores import find_candidates
from pyBiodatafuse.algorithms.DREAMwalk.constant import (
    DDA_DIRECTORY,
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


def preprocess_files():
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
):
    """Generate similarity graph.
    :param drug_hierarchy_file: Drug hierarchy file path
    :param dis_sim_file: Disease similarity file path
    :param drug_sim_file: Drug similarity file path
    :param network_file: Network file path
    :param output_file: Final similarity file path
    """
    if os.path.exists(drug_sim_file):
        drug_sim = pd.read_csv(drug_sim_file, sep="\t", header=None)
    else:
        save_sim_graph(
            networkf=network_file, hierf=drug_hierarchy_file, outputf=drug_sim_file, cutoff=0.4
        )
        logger.info("Similarity graph generated successfully.")

    # create similarity graph merging disease and drugs
    drug_sim = pd.read_csv(drug_sim_file, sep="\t", header=None)
    dis_sim = pd.read_csv(dis_sim_file, sep="\t", header=None)

    similarty_graph = pd.concat([dis_sim, drug_sim], ignore_index=True)
    similarty_graph[4] = range(0, len(similarty_graph))  # reindexing the counts
    similarty_graph.to_csv(output_file, sep="\t", index=False, header=False)


def generate_node_embeddings(
    similarity_file: str,
    node_type_file: str,
    network_file: str,
    output_file: str,
):
    """Generate node embeddings."""

    save_embedding_files(
        netf=network_file,
        sim_netf=similarity_file,
        outputf=output_file,
        nodetypef=node_type_file,
        tp_factor=0.3,
    )
    logger.info("Node embeddings generated successfully.")


def predict_dda():
    """Predict drug-disease association."""
    pairf = "dreamwalk_data/dda1.tsv"  # TODO: change to actual file
    modelf = "dreamwalk_data/clf1.pkl"
    embeddingf = "dreamwalk_data/embedding_file.pkl"

    predict_dda(embeddingf=embeddingf, pairf=pairf, modelf=modelf)
    logger.info("Drug-disease association predicted successfully.")


def find_candidates():
    """Calculate drug scores."""
    embeddingf = "dreamwalk_data/embedding_file.pkl"
    model_folder = "dreamwalk_data"
    query_disease = "C00000"  # TODO: change to actual disease code
    kgfile = "dreamwalk_data/preprocessed_graph.csv"  # TODO: change to actual file

    results1 = find_candidates(
        kgfile, embeddingf, model_folder, query_disease, candidates_count=400
    )
    results1.to_csv("dreamwalk_data/results.csv")
    logger.info("Drug scores calculated successfully.")


if __name__ == "__main__":
    networkf = f"dreamwalk_data/{GRAPH_FILE}"
    hierarchyf = f"dreamwalk_data/{DRUG_HIERACHY_FILE}"
    drug_simf = f"dreamwalk_data/{DRUG_SIM_FILE}"
    disease_simf = f"dreamwalk_data/{DISEASE_SIM_FILE}"
    nodetypef = f"dreamwalk_data/{NODE_TYPE_FILE}"
    embeddingf = f"dreamwalk_data/{EMBEDDING_FILE}"
    simf = f"dreamwalk_data/{SIM_FILE}"

    # preprocess_files()
    # generate_sim(
    #     drug_hierarchy_file=hierarchyf,
    #     dis_sim_file=disease_simf,
    #     drug_sim_file=drug_simf,
    #     network_file=networkf,
    #     output_file=simf,
    # )
    generate_node_embeddings(
        similarity_file=simf,
        node_type_file=nodetypef,
        network_file=networkf,
        output_file=embeddingf,
    )
    # predict_dda()
    # find_candidates()
