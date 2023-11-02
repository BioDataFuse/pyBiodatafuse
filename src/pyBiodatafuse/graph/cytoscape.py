# coding: utf-8

"""Python module to export a NetworkX graph to cytoscape-compliant format and set the styling for cytoscape."""

import json

import networkx as nx


def replace_graph_attrs_for_cytoscape(g: nx.MultiDiGraph):
    """Adapt the node and edge attributes keys to the cytoscape json structure.

    :param g: input NetworkX graph object.
    :returns: output NetworkX graph object.
    """
    for node in list(g.nodes()):
        for k, v in list(g.nodes[node].items()):
            if k == "labels":
                nx.set_node_attributes(g, {node: {"name": v}})
                del g.nodes[node]["labels"]
            elif k == "source":
                nx.set_node_attributes(g, {node: {"datasource": v}})
                del g.nodes[node]["source"]

    for u, v, k in list(g.edges(keys=True)):
        for x, y in list(g[u][v][k].items()):
            if x == "label":
                g[u][v][k]["interaction"] = y
                del g[u][v][k]["label"]

    return g


def convert_graph_to_cytoscape_json(g: nx.MultiDiGraph):
    """Convert a NetworkX graph to cytoscape json file.

    :param g: the NetworkX graph object.
    :returns: a cytoscape network json object.
    """
    adj_g = replace_graph_attrs_for_cytoscape(g)

    cytoscape_graph = nx.cytoscape_data(adj_g)

    return cytoscape_graph


def write_cytoscape_graph_to_file(cytoscape_graph: dict, output_path: str):
    """Write cytoscape graph to json file.

    :param cytoscape_graph: the cytoscape graph object.
    :param output_path: the output path.
    """
    with open(output_path, "w") as out:
        json.dump(cytoscape_graph, out)
