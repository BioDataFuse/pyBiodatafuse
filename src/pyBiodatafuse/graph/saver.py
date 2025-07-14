# coding: utf-8

"""Python module to save the graph to multiple formats."""

import json
import os
import pickle
from logging import Logger
from typing import Any, Dict

import networkx as nx
import pandas as pd

from pyBiodatafuse.graph.generator import build_networkx_graph

logger = Logger(__name__)


def save_graph_to_graphml(g: nx.MultiDiGraph, output_path: str):
    """Convert a NetworkX graph to Neo4J graphml file.

    :param g: the NetworkX graph object.
    :param output_path: the output path of the graphml file
    """
    nx.write_graphml(g, output_path, named_key_ids=True)


def save_graph_to_edgelist(g: nx.MultiDiGraph, output_path: str):
    """Convert a NetworkX graph to edgelist file.

    :param g: the NetworkX graph object.
    :param output_path: the output path of the edgelist file
    """
    nx.write_edgelist(
        g,
        output_path,
        delimiter=" ",
        data=True,
        encoding="utf-8",
    )


def save_cytoscape_json(graph: dict, output_path: str):
    """Write Cytoscape JSON graph to file.

    :param graph: the Cytoscape JSON graph object.
    :param output_path: the output path.
    """
    with open(output_path, "w") as out:
        json.dump(graph, out)


def save_graph_to_tsv(g, output_dir="output"):
    """Save the graph to TSV files for nodes and edges.

    :param g: the input graph to save.
    :param output_dir: the directory to save the TSV files.
    """
    os.makedirs(output_dir, exist_ok=True)

    # Save nodes
    nodes_path = os.path.join(output_dir, "nodes.tsv")
    with open(nodes_path, "w") as f:
        f.write("node_id\tattributes\n")
        for node, attrs in g.nodes(data=True):
            f.write(f"{node}\t{json.dumps(attrs)}\n")

    # Save edges
    edges_path = os.path.join(output_dir, "edges.tsv")
    with open(edges_path, "w") as f:
        f.write("source\ttarget\tkey\tattributes\n")
        for u, v, k, attrs in g.edges(keys=True, data=True):
            f.write(f"{u}\t{v}\t{k}\t{json.dumps(attrs)}\n")


def save_graph(
    combined_df: pd.DataFrame,
    combined_metadata: Dict[Any, Any],
    disease_compound: pd.DataFrame = None,
    graph_name: str = "combined",
    graph_dir: str = "examples/usecases/",
):
    """Save the graph to a file.

    :param combined_df: the input df to be converted into a graph.
    :param combined_metadata: the metadata of the graph.
    :param disease_compound: the input df containing disease-compound relationships.
    :param graph_name: the name of the graph.
    :param graph_dir: the directory to save the graph.
    :returns: a NetworkX MultiDiGraph

    """
    os.makedirs(graph_dir, exist_ok=True)

    df_path = f"{graph_dir}/{graph_name}_df.pkl"
    metadata_path = f"{graph_dir}/{graph_name}_metadata.pkl"
    graph_path_pickle = f"{graph_dir}/{graph_name}_graph.pkl"
    graph_path_gml = f"{graph_dir}/{graph_name}_graph.gml"
    graph_path_edgelist = f"{graph_dir}/{graph_name}_graph.edgelist"
    # Save the combined DataFrame
    combined_df.to_pickle(df_path)
    logger.warning(f"Combined DataFrame saved in {df_path}")

    # Save the metadata
    with open(metadata_path, "wb") as file:
        pickle.dump(combined_metadata, file)
    logger.warning(f"Metadata saved in {metadata_path}")

    # Save the graph
    g = build_networkx_graph(combined_df, disease_compound)
    logger.warning("Graph is built successfully")

    with open(graph_path_pickle, "wb") as f:
        pickle.dump(g, f)

    save_graph_to_graphml(g, graph_path_gml)  # for neo4j import
    logger.warning(f"Graph saved in: \n {graph_path_pickle} \n {graph_path_gml}")

    save_graph_to_edgelist(g, graph_path_edgelist)
    logger.warning(f"Graph saved in {graph_path_edgelist}")

    return g
