"""Codes taken from DreamWalk repository: https://github.com/eugenebang/DREAMwalk"""

import logging
import os

import networkx as nx
import pandas as pd
from tqdm import tqdm

from pyBiodatafuse.analyzer.summarize import BioGraph

logger = logging.getLogger(__name__)


def get_drug_disease_file(graph, output_dir: str):
    """Generate drug-disease association (DDA) files.
    :param graph: networkx.Graph object
    :param output_dir: output directory to save the files
    :return: None
    """
    drug_nodes = [node for node, label in graph.nodes(data="labels") if label == "Compound"]

    disease_nodes = [node for node, label in graph.nodes(data="labels") if label == "Disease"]

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
    os.makedirs(f"{output_dir}/dda", exist_ok=True)
    for i in range(1, 11):
        sampled = unknown_dda.sample(n=len(known_dda))
        sampled.to_csv(f"{output_dir}/dda/dda_{i}.tsv", sep="\t", index=False)


def create_nodetype_file(graph, output_dir: str) -> None:
    """Generate node type file for the algorithm.
    :param graph: networkx.Graph object
    :param output_dir: output directory to save the files
    :return: None
    """
    tmp = []
    for node, node_label in graph.nodes(data="labels"):
        tmp.append({"node": node, "type": node_label})

    tmp_df = pd.DataFrame(tmp)
    tmp_df["type"] = tmp_df["type"].map({"Gene": "gene", "Disease": "disease", "Compound": "drug"})
    tmp_df = tmp_df.sort_values(by="type")
    tmp_df.to_csv(f"{output_dir}/nodetypes.tsv", sep="\t", index=False)
    return None


def create_network_file(graph, output_dir: str) -> None:
    """Generate network file for the algorithm.
    :param graph: networkx.Graph object
    :param output_dir: output directory to save the files
    :return: None
    """
    rel_to_id = {"activates": 1, "inhibits": 1, "associated_with": 2, "interacts_with": 3}

    graph_data = []

    for source, target, edge in graph.edges(data=True):
        if edge["label"] not in rel_to_id:
            continue

        edge_id = rel_to_id[edge["label"]]

        graph_data.append(
            {
                "source": source,
                "target": target,
                "edgetype": edge_id,
                "weight": 1,
            }
        )

    output_graph = pd.DataFrame(graph_data)
    output_graph["edge_counter"] = range(1, len(output_graph) + 1)
    output_graph.to_csv(f"{output_dir}/graph.txt", sep="\t", index=False, header=False)
    return None


def check_hierarchy(graph: BioGraph, data_dir: str) -> set:
    """Check if the node has a hierarchy.
    :param graph: BioGraph object
    :param data_dir: Data directory to save the files
    :return: remove_nodes: Set of nodes to remove
    """
    remove_nodes = set()

    # ATC hierarchy
    atc_hierarchy = pd.read_csv(
        f"{data_dir}/atc_hierarchy.csv", usecols=["dbID", "atcClassification", "id"]
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

    # MONDO hierarchy
    mondo_hierarchy = pd.read_csv(
        f"{data_dir}/mondo_hierarchy.csv"
    )  # file created using graph nodes itself

    for disease_id, _, parents in mondo_hierarchy.values:
        if pd.isna(parents) or isinstance(parents, float):
            remove_nodes.add(disease_id)
        elif parents == "disease":
            remove_nodes.add(disease_id)

    logger.warning(f"Nodes to remove: {len(remove_nodes)}")
    return remove_nodes


def get_drug_hierarchy(graph: BioGraph, data_dir: str) -> pd.DataFrame:
    """Generating the drug hierarchy using ATC classification.
    :param graph: BioGraph object
    :param data_dir: Data directory to save the files
    :return: drug_hierarchy_df: DataFrame
    """
    drug_classes = []

    atc_hierarchy = pd.read_csv(
        f"{data_dir}/atc_hierarchy.csv", usecols=["dbID", "atcClassification", "id"]
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

            parents = atc_classification_list[:-1]
            children = atc_classification_list[1:]
            for parent, child in zip(parents, children):
                drug_classes.append({"child": child, "parent": parent})

    drug_class_df = pd.DataFrame(drug_classes)
    drug_class_df.drop_duplicates(inplace=True)

    # make "Drug" into lowercase
    drug_class_df["child"] = drug_class_df["child"].apply(lambda x: x.lower() if x == "Drug" else x)

    return drug_class_df


def get_disease_hierarchy(data_dir: str) -> pd.DataFrame:
    """Generating the disease hierarchy using MONDO classification.
    :param output_dir: output directory to save the files
    :return: disease_hierarchy_df: DataFrame
    """
    mondo_hierarchy = pd.read_csv(
        f"{data_dir}/mondo_hierarchy.csv"
    )  # file created using graph nodes itself

    disease_hierarchy = []

    for disease_id, _, parents in mondo_hierarchy.values:
        if pd.isna(parents) or isinstance(parents, float) or parents == "disease":
            continue

        all_parents = parents.split(",")
        parent = all_parents
        children = [disease_id] + list(reversed(all_parents[:-1]))
        for p, c in zip(parent, children):
            disease_hierarchy.append({"child": c, "parent": p})

    hierarchy_df = pd.DataFrame(disease_hierarchy)
    hierarchy_df.drop_duplicates(inplace=True)
    return hierarchy_df


def create_files(graph_obj: BioGraph, output_dir: str = "./dreamwalk_data") -> None:
    """Generate files for the DREAMwalk algorithm from the BioGraph object.
    :param graph_obj: BioGraph object
    :param output_dir: output directory to save the files
    :return: None
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
    nx.write_gml(updated_graph, f"{output_dir}/subgraph_graph.gml")

    # Node type label file - TSV files with two columns: node_id, node_type
    create_nodetype_file(updated_graph, output_dir)

    # Network file - TXT file with four columns: source, target, edgetype, weight
    create_network_file(updated_graph, output_dir)

    # Hierarchy file - CSV file with two columns: child, parent
    drug_hierarchy_df = get_drug_hierarchy(updated_graph, output_dir)
    disease_hierarchy_df = get_disease_hierarchy(output_dir)
    hierarchy_df = pd.concat([drug_hierarchy_df, disease_hierarchy_df], ignore_index=True)
    hierarchy_df.to_csv(f"{output_dir}/hierarchy.csv", index=False)

    # Positive/negative drug-disease association file - TSV file with three columns: drug, disease, label
    get_drug_disease_file(updated_graph, output_dir)
