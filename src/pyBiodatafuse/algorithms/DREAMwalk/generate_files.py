"""Codes taken from DreamWalk repository: https://github.com/eugenebang/DREAMwalk."""

import logging
import os
from itertools import combinations
from typing import Any, Dict, List

import networkx as nx
import pandas as pd
from tqdm import tqdm

from pyBiodatafuse.algorithms.DREAMwalk.constant import (
    ATC_HIERARCHY_FILE,
    DDA_DIRECTORY,
    DISEASE_SIM_FILE,
    DRUG_HIERACHY_FILE,
    GRAPH_FILE,
    NODE_TYPE_FILE,
    SUBGRAPH_FILE,
)
from pyBiodatafuse.algorithms.DREAMwalk.utils import read_graph
from pyBiodatafuse.analyzer.summarize import BioGraph

logger = logging.getLogger(__name__)


def jaccard_similarity(set1: set, set2: set) -> float:
    """Calculate Jaccard similarity between two lists.

    :param set1: Set 1
    :param set2: Set 2
    :returns: Jaccard similarity
    """
    intersection = len(set1.intersection(set2))
    union = len(set1.union(set2))
    return intersection / union


def get_drug_disease_file(output_dir: str) -> None:
    """Generate drug-disease association (DDA) files.

    :param output_dir: output directory to save the files
    """
    graph = read_graph(f"{output_dir}/{GRAPH_FILE}", weighted=True, directed=False, delimiter="\t")
    nodetype_df = pd.read_csv(f"{output_dir}/{NODE_TYPE_FILE}", sep="\t")

    drug_nodes = nodetype_df[nodetype_df["type"] == "drug"]["node"].values
    disease_nodes = nodetype_df[nodetype_df["type"] == "disease"]["node"].values

    drug_disease_edges = []

    for drug in tqdm(drug_nodes, desc="Creating Drug-Disease Association files"):
        for disease in disease_nodes:
            if graph.has_edge(drug, disease):
                drug_disease_edges.append({"drug": drug, "disease": disease, "label": 1})
            else:
                drug_disease_edges.append({"drug": drug, "disease": disease, "label": 0})

    tmp = pd.DataFrame(drug_disease_edges)

    known_dda = tmp[tmp["label"] == 1]
    unknown_dda = tmp[tmp["label"] == 0]

    logger.warning(f"Known DDA: {len(known_dda)}, Unknown DDA: {len(unknown_dda)}")

    # Sampling same number of negative samples as positive samples
    os.makedirs(f"{output_dir}/{DDA_DIRECTORY}", exist_ok=True)
    for i in range(1, 11):
        sampled = unknown_dda.sample(n=len(known_dda))
        sampled = pd.concat([known_dda, sampled])
        sampled.to_csv(f"{output_dir}/{DDA_DIRECTORY}/dda_{i}.tsv", sep="\t", index=False)


def create_nodetype_file(graph: nx.MultiDiGraph, output_dir: str) -> None:
    """Generate node type file for the algorithm.

    :param graph: networkx.Graph object
    :param output_dir: output directory to save the files
    :returns: None
    """
    tmp = []
    for node, node_label in graph.nodes(data="labels"):
        tmp.append({"node": node, "type": node_label})

    tmp_df = pd.DataFrame(tmp)
    tmp_df["type"] = tmp_df["type"].map({"Gene": "gene", "Disease": "disease", "Compound": "drug"})
    tmp_df = tmp_df.sort_values(by="type")
    tmp_df.to_csv(f"{output_dir}/{NODE_TYPE_FILE}", sep="\t", index=False)
    return None


def create_network_file(graph: nx.MultiDiGraph, output_dir: str) -> None:
    """Generate network file for the algorithm.

    :param graph: networkx.Graph object
    :param output_dir: output directory to save the files
    :returns: None
    """
    rel_to_id = {
        "activates": 1,
        "inhibits": 1,
        "active": 1,
        "inactive": 1,
        "associated_with": 2,
        "interacts_with": 3,
        "treats": 4,
    }

    graph_data = []

    counter = 1
    for source, target, edge in graph.edges(data=True):
        if edge["label"] not in rel_to_id:
            continue

        edge_id = rel_to_id[edge["label"]]
        graph_data.append(
            {
                "source": source,
                "target": target,
                "type": edge_id,
                "weight": 1,
                "edge_id": counter,
            }
        )
        counter += 1

    output_graph = pd.DataFrame(graph_data)
    output_graph.to_csv(f"{output_dir}/{GRAPH_FILE}", sep="\t", index=False)
    return None


def check_hierarchy(graph: nx.MultiDiGraph, data_dir: str) -> set:
    """Check if the node has a hierarchy.

    :param graph: Networkx object
    :param data_dir: Data directory to save the files
    :returns: Set of nodes to remove
    """
    remove_nodes = set()

    # ATC hierarchy
    atc_hierarchy = pd.read_csv(
        f"{data_dir}/{ATC_HIERARCHY_FILE}", usecols=["dbID", "atcClassification", "id"]
    )

    atc_hierarchy = atc_hierarchy.rename(
        columns={"dbID": "drug", "atcClassification": "atc", "id": "chembl_id"}
    )
    atc_hierarchy["drug"] = "DrugBank:" + atc_hierarchy["drug"]
    atc_hierarchy["chembl_id"] = "ChEMBL:" + atc_hierarchy["chembl_id"]

    drug_nodes = [node for node, label in graph.nodes(data="labels") if label == "Compound"]

    for drug in drug_nodes:
        node_info = graph.nodes[drug]
        if "drugbank_id" in node_info:
            atc_classes = atc_hierarchy[atc_hierarchy["drug"] == node_info["drugbank_id"]][
                "atc"
            ].values
        elif "chembl_id" in node_info:
            atc_classes = atc_hierarchy[atc_hierarchy["chembl_id"] == node_info["chembl_id"]][
                "atc"
            ].values
        else:
            atc_classes = []

        if len(atc_classes) == 0:
            remove_nodes.add(drug)

    logger.warning(f"Nodes to remove: {len(remove_nodes)}")
    return remove_nodes


def get_drug_hierarchy(graph: nx.MultiDiGraph, data_dir: str, output_dir: str) -> None:
    """Generate the drug hierarchy using ATC classification.

    :param graph: Networkx object
    :param data_dir: Data directory to save the files
    :param output_dir: output directory to save the files
    :returns: None
    """
    drug_classes = []

    atc_hierarchy = pd.read_csv(
        f"{data_dir}/{ATC_HIERARCHY_FILE}", usecols=["dbID", "atcClassification", "id"]
    )

    atc_hierarchy = atc_hierarchy.rename(
        columns={"dbID": "drug", "atcClassification": "atc", "id": "chembl_id"}
    )
    atc_hierarchy["drug"] = "DrugBank:" + atc_hierarchy["drug"]
    atc_hierarchy["chembl_id"] = "ChEMBL:" + atc_hierarchy["chembl_id"]

    # With hierarch of drugs in KG
    drug_nodes = [node for node, label in graph.nodes(data="labels") if label == "Compound"]

    for drug in drug_nodes:
        node_info = graph.nodes[drug]
        if "drugbank_id" in node_info:
            atc_classes = atc_hierarchy[atc_hierarchy["drug"] == node_info["drugbank_id"]][
                "atc"
            ].values[0]
        else:
            atc_classes = atc_hierarchy[atc_hierarchy["chembl_id"] == node_info["chembl_id"]][
                "atc"
            ].values[0]

        pnames = atc_classes.split(";")

        for atc_classification in pnames:
            atc_classification_list = atc_classification.split(",")
            drug_classes.append({"child": drug, "parent": atc_classification_list[0]})

            parents = atc_classification_list[1:]
            children = atc_classification_list[:-1]

            for child, parent in zip(children, parents):
                drug_classes.append({"child": child, "parent": parent})

    drug_class_df = pd.DataFrame(drug_classes)
    drug_class_df.drop_duplicates(inplace=True)
    drug_class_df.to_csv(f"{output_dir}/{DRUG_HIERACHY_FILE}", index=False)
    return None


def get_disease_similarity(graph: nx.MultiDiGraph, data_dir: str) -> None:
    """Generate disease similarity graph based on gene overlap.

    :param graph: Networkx graph object
    :param data_dir: Data directory to save the files
    :returns: None
    """
    disease_nodes = [node for node, label in graph.nodes(data="labels") if label == "Disease"]
    disease_pairs = list(combinations(disease_nodes, 2))

    disease_similarity_matrix = []  # type: List[Dict[str, Any]]

    for disease_1, disease_2 in tqdm(disease_pairs, desc="Calculating Disease Similarity"):
        gene_list_1 = set([g for g, d in graph.in_edges(disease_1)])
        gene_list_2 = set([g for g, d in graph.in_edges(disease_2)])

        similarity = jaccard_similarity(gene_list_1, gene_list_2)
        if similarity < 0.4:
            continue

        disease_similarity_matrix.append(
            {
                "source": disease_1,
                "target": disease_2,
                "type": 2,
                "weight": similarity,
                "edge_id": len(disease_similarity_matrix),
            }
        )

    sim_graph = pd.DataFrame(disease_similarity_matrix)
    sim_graph.to_csv(f"{data_dir}/{DISEASE_SIM_FILE}", sep="\t", index=False, header=False)
    return None


def create_files(graph_obj: BioGraph, output_dir: str = "./dreamwalk_data") -> None:
    """Generate files for the DREAMwalk algorithm from the BioGraph object.

    :param graph_obj: BioGraph object
    :param output_dir: output directory to save the files
    """
    subgraph = graph_obj.get_subgraph(node_types=["Gene", "Disease", "Compound"])
    logger.warning(f"Subgraph nodes: {len(subgraph.nodes())}, edges: {len(subgraph.edges())}")

    # Removing nodes with no heirarchy
    remove_nodes = check_hierarchy(subgraph, output_dir)
    updated_graph = nx.MultiDiGraph(subgraph)

    for node in remove_nodes:
        updated_graph.remove_node(node)

    logger.info(
        f"Updated graph nodes: {len(updated_graph.nodes())}, edges: {len(updated_graph.edges())}"
    )
    nx.write_gml(updated_graph, f"{output_dir}/{SUBGRAPH_FILE}")

    # Node type label file - TSV files with two columns: node_id, node_type
    create_nodetype_file(updated_graph, output_dir)

    # Network file - TXT file with four columns: source, target, edgetype, weight
    create_network_file(updated_graph, output_dir)

    # Hierarchy file - CSV file with two columns: child, parent
    get_drug_hierarchy(updated_graph, data_dir=output_dir, output_dir=output_dir)

    # Disease similairty file - CSV file with two columns: child, parent
    get_disease_similarity(updated_graph, output_dir)

    # Positive/negative drug-disease association file - TSV file with three columns: drug, disease, label
    get_drug_disease_file(output_dir)
