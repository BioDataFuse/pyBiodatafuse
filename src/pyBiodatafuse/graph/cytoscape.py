# coding: utf-8

"""Python module to export a NetworkX graph to cytoscape-compliant format and set the styling for cytoscape."""

import networkx as nx
import py4cytoscape as p4c

import pyBiodatafuse.constants as Cons


def _replace_graph_attrs(g: nx.MultiDiGraph):
    """Adapt the node and edge attributes keys to the cytoscape json structure.


    :param g: input NetworkX graph object.
    :returns: output NetworkX graph object.
    """
    # Fix node attributes
    for node in list(g.nodes()):
        attrs = g.nodes[node]
        # Set 'id' explicitly
        attrs['id'] = node
        if 'name' in attrs:
            attrs['label'] = attrs['name']
            del attrs['name']
        if 'source' in attrs:
            attrs['datasource'] = attrs['source']
            del attrs['source']

    # Fix edge attributes
    for u, v, k in g.edges(keys=True):
        attrs = g[u][v][k]
        # Edge source/target must match node IDs
        attrs['source'] = u
        attrs['target'] = v
        if 'label' in attrs:
            attrs['interaction'] = attrs['label']
            del attrs['label']
        if 'source' in attrs:
            attrs['datasource'] = attrs['source']
            del attrs['source']

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
    column = Cons.LABEL
    values = [
        Cons.GENE_NODE_LABEL,
        Cons.ANATOMICAL_NODE_LABEL,
        Cons.DISEASE_NODE_LABEL,
        Cons.GO_BP_NODE_LABEL,
        Cons.GO_MF_NODE_LABEL,
        Cons.GO_CC_NODE_LABEL,
        Cons.PATHWAY_NODE_LABEL,
        Cons.COMPOUND_NODE_LABEL,
        Cons.SIDE_EFFECT_NODE_LABEL,
        Cons.HOMOLOG_NODE_LABEL,
        Cons.KEY_EVENT_NODE_LABEL,
        Cons.MIE_NODE_LABEL,
        Cons.AOP_NODE_LABEL,
        Cons.AO_NODE_LABEL,
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
