# coding: utf-8

"""Python module to export a NetworkX graph to cytoscape-compliant format and set the styling for cytoscape."""

import networkx as nx
import py4cytoscape as p4c

from pyBiodatafuse.constants import (
    ANATOMICAL_NODE_LABELS,
    AO_NODE_LABELS,
    AOP_NODE_LABELS,
    COMPOUND_NODE_LABELS,
    DISEASE_NODE_LABELS,
    GENE_NODE_LABELS,
    GO_BP_NODE_LABELS,
    GO_CC_NODE_LABELS,
    GO_MF_NODE_LABELS,
    HOMOLOG_NODE_LABELS,
    KEY_EVENT_NODE_LABELS,
    MIE_NODE_LABELS,
    PATHWAY_NODE_LABELS,
    SIDE_EFFECT_NODE_LABELS,
)


def _replace_graph_attrs(g: nx.MultiDiGraph):
    """Adapt the node and edge attributes keys to the cytoscape json structure.

    :param g: input NetworkX graph object.
    :returns: output NetworkX graph object.
    """
    for node in list(g.nodes()):
        for k, v in list(g.nodes[node].items()):
            if k == "name":
                nx.set_node_attributes(g, {node: {"label": v}})
                del g.nodes[node]["name"]
            elif k == "source":
                nx.set_node_attributes(g, {node: {"datasource": v}})
                del g.nodes[node]["source"]

    for u, v, k in list(g.edges(keys=True)):
        for x, y in list(g[u][v][k].items()):
            if x == "label":
                g[u][v][k]["interaction"] = y
                del g[u][v][k]["label"]
            elif x == "source":
                g[u][v][k]["datasource"] = y
                del g[u][v][k]["source"]

    return g


def convert_graph_to_json(g: nx.MultiDiGraph):
    """Convert a NetworkX graph to cytoscape json file.

    :param g: the NetworkX graph object.
    :returns: a cytoscape network json object.
    """
    adj_g = _replace_graph_attrs(g)

    cytoscape_graph = nx.cytoscape_data(adj_g)

    return cytoscape_graph


def load_graph(g: nx.MultiDiGraph, network_name: str):
    """Load the obtained graph into a running instance of Cytoscape.

    :param g: input NetworkX graph object.
    :param network_name: Network name to appear in Cytoscape.
    """
    adj_g = _replace_graph_attrs(g)

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
    column = "labels"
    values = [
        GENE_NODE_LABELS,
        ANATOMICAL_NODE_LABELS,
        DISEASE_NODE_LABELS,
        GO_BP_NODE_LABELS,
        GO_MF_NODE_LABELS,
        GO_CC_NODE_LABELS,
        PATHWAY_NODE_LABELS,
        COMPOUND_NODE_LABELS,
        SIDE_EFFECT_NODE_LABELS,
        HOMOLOG_NODE_LABELS,
        KEY_EVENT_NODE_LABELS,
        MIE_NODE_LABELS,
        AOP_NODE_LABELS,
        AO_NODE_LABELS,
    ]
    shapes = [
        "ELLIPSE",  # Genes
        "HEXAGON",  # Anatomical
        "VEE",  # Diseases
        "PARALLELOGRAM",  # GO BP
        "ROUND_RECTANGLE",  # GO MF
        "RECTANGLE",  # GO CC
        "OCTAGON",  # Pathways
        "DIAMOND",  # Compounds
        "TRIANGLE",  # Side Effects
        "Ellipse",  # Homologs
        "TRIANGLE",  # Key Events
        "TRIANGLE",  # MIE
        "VEE",  # AOP
        "OCTAGON",  # AO
    ]
    colors = [
        "#42d4f4",  # Cyan for Genes
        "#4363d8",  # Blue for Anatomical
        "#e6194B",  # Red for Diseases
        "#ff7b00",  # Orange for GO BP
        "#ffa652",  # Orange for GO MF
        "#ffcd90",  # Orange for GO CC
        "#3cb44b",  # Green for Pathways
        "#ffd700",  # Gold for Compounds
        "#aaffc3",  # Mint for Side Effects
        "#9b59b6",  # Purple for Homologs
        "#aaffc3",  # Mint for Key Events
        "#3cb44b",  # Green for MIE
        "#000075",  # Navy for AOP
        "#e6194B",  # Red for AO
    ]

    # Apply node shape and color mappings
    p4c.set_node_color_mapping(
        column,
        values,
        colors,
        mapping_type="d",
        style_name="default",
    )

    p4c.set_node_shape_mapping(column, values, shapes, style_name="default")
