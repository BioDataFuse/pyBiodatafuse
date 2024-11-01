"""Codes taken from DreamWalk repository: https://github.com/eugenebang/DREAMwalk."""

import logging
import os
import pickle
from typing import Dict

import pandas as pd

from pyBiodatafuse.algorithms.DREAMwalk.utils import read_graph

logger = logging.getLogger(__name__)


def process_drugs(drugs: set, disease: str, embeddingf: str, model_list: list) -> pd.DataFrame:
    """Process drugs to find candidates.

    :param drugs: Set of drugs to process
    :param disease: Query disease node
    :param embeddingf: Embedding file path
    :param model_list: List of models
    :returns: Dataframe of candidate drugs
    """
    with open(embeddingf, "rb") as fin:
        embedding_dict = pickle.load(fin)

    results = []

    for drug in drugs:
        prob_sum = []

        for i in range(len(model_list)):
            model_prob_sum = model_list[i].predict_proba(
                [embedding_dict[drug] - embedding_dict[disease]]
            )[:, 1][0]
            prob_sum.append(model_prob_sum)

        prob_avg = sum(prob_sum) / len(prob_sum)
        results.append({"drug": drug, "avg_prob": prob_avg})

    result_df = pd.DataFrame(results)
    candidates = result_df.sort_values(by="avg_prob", ascending=False)

    return candidates


def find_candidates(
    networkf: str,
    nodetypef: str,
    embeddingf: str,
    model_folder: str,
    query_disease: str,
    candidate_count: int = 0,
    weighted: bool = True,
    directed: bool = False,
    net_delimiter: str = "\t",
) -> pd.DataFrame:
    """Find predictions for non-associated drugs.

    :param networkf: Network file path
    :param nodetypef: Node type file path
    :param embeddingf: Embedding file path
    :param model_folder: Folder where all models are stored
    :param query_disease: Query disease node
    :param candidate_count: Number of candidate drugs to find
    :param weighted: Whether the network is weighted or not
    :param directed: Whether the network is directed or not
    :param net_delimiter: Delimiter for the network file
    :returns: Dataframe of candidate drugs
    """
    graph = read_graph(networkf, weighted=weighted, directed=directed, delimiter=net_delimiter)
    node_types = pd.read_csv(nodetypef, sep="\t")

    drug_nodes = node_types[node_types["type"] == "drug"]["node"].values
    disease_nodes = node_types[node_types["type"] == "disease"]["node"].values

    assert query_disease in disease_nodes, "Query disease not found in the network."

    # Find connected and non-connected drugs
    drug_disease_links = {
        "connected": set(),
        "non_connected": set(),
    }  # type: Dict[str, set]

    for drug in drug_nodes:
        if graph.has_edge(drug, query_disease):
            drug_disease_links["connected"].add(drug)
        else:
            drug_disease_links["non_connected"].add(drug)

    logging.info(
        f"Found {len(drug_disease_links['connected'])} connected drugs and {len(drug_disease_links['non_connected'])} unconnected drugs."
    )

    # Load all the different models
    model_list = []

    model_files = os.listdir(model_folder)
    for file in model_files:
        if not file.endswith(".pkl"):
            continue

        file_path = os.path.join(model_folder, file)
        with open(file_path, "rb") as f:
            model = pickle.load(f)

        # Process content as needed
        model_list.append(model)

    logger.info(f"Found {len(model_list)} loaded models.")

    candidate_drugs = process_drugs(
        drugs=drug_disease_links["non_connected"],
        disease=query_disease,
        embeddingf=embeddingf,
        model_list=model_list,
    )

    if candidate_count > 0:
        return candidate_drugs.head(candidate_count)

    return candidate_drugs
