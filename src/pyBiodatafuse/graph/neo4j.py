# coding: utf-8

"""Python module to export a NetworkX graph to Neo4J graphml file."""

import networkx as nx


def write_graph_to_neo4j_graphml(g: nx.MultiDiGraph, output_path: str):
    """Convert a NetworkX graph to Neo4J graphml file.

    :param g: the NetworkX graph object.
    :param output_path: the output path of the graphml file
    """
    nx.write_graphml(g, output_path, named_key_ids=True)
