# coding: utf-8

"""Python module to export a NetworkX graph to cytoscape-compliant format and set the styling for cytoscape."""

import json
from importlib import resources
from pathlib import Path

import networkx as nx
import py4cytoscape as p4c


def _replace_graph_attrs_for_cytoscape(g: nx.MultiDiGraph):
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


def convert_graph_to_cytoscapejs_json(g: nx.MultiDiGraph):
    """Convert a NetworkX graph to cytoscape json file.

    :param g: the NetworkX graph object.
    :returns: a cytoscape network json object.
    """
    adj_g = _replace_graph_attrs_for_cytoscape(g)

    cytoscape_graph = nx.cytoscape_data(adj_g)

    return cytoscape_graph


def save_graph_to_cytoscape_graphml(g: nx.MultiDiGraph, output_path: str, export_style=True):
    """Convert a NetworkX graph to cytoscape graphml file.

    :param g: the NetworkX graph object.
    :param output_path: the output path of the graphml file
    :param export_style: option to export styles.xml along the cytoscape graph file
    """
    adj_g = _replace_graph_attrs_for_cytoscape(g)

    nx.write_graphml(adj_g, output_path, named_key_ids=True)

    if export_style:
        styles = ""

        with resources.path("pyBiodatafuse.resources", "styles.xml") as style_file:
            with style_file.open(mode="r", encoding="utf-8") as style_file_stream:
                styles = style_file_stream.read()

        if styles != "":
            output_dir = Path(output_path).parent.absolute()

            with open(str(output_dir) + "/styles.xml", "w") as out:
                out.write(styles)


def load_graph_into_cytoscape(g: nx.MultiDiGraph, network_name: str):
    """Load the obtained graph into a running instance of Cytoscape.

    :param g: input NetworkX graph object.
    :param network_name: Network name to appear in Cytoscape.
    """
    adj_g = _replace_graph_attrs_for_cytoscape(g)

    # Define the visual style as a dictionary
    default = {
        "title": "BioDataFuse_style",
        "defaults": [
            {"visualProperty": "NODE_FILL_COLOR", "value": "#FF0000"},
            {"visualProperty": "EDGE_COLOR", "value": "#000000"},
        ],
        "mappings": [],
    }

    # Create the network in Cytoscape
    p4c.networks.create_network_from_networkx(
        adj_g,
        title=network_name,
        collection="BioDataFuse",
    )

    # Apply the visual style
    p4c.styles.create_visual_style(default)

    # Define node shape and color mapping
    column = "node_type"
    values = [
        "gene",
        "disease",
        "gene ontology",
        "reactome pathways",
        "drug interactions",
    ]
    shapes = ["DIAMOND", "RECTANGLE", "OCTAGON", "HEXAGON", "ELLIPSE"]
    colors = ["#AAFF88", "#B0C4DE", "pink", "yellow", "red"]

    # Apply node shape and color mappings
    p4c.set_node_color_mapping(column, values, colors, mapping_type="d", style_name="default")

    p4c.set_node_shape_mapping(column, values, shapes, style_name="default")


def save_cytoscape_json_to_file(cytoscape_graph: dict, output_path: str):
    """Write cytoscape graph to json file.

    :param cytoscape_graph: the cytoscape graph object.
    :param output_path: the output path.
    """
    with open(output_path, "w") as out:
        json.dump(cytoscape_graph, out)
