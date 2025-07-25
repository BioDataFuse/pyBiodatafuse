# coding: utf-8

"""Python module to export a NetworkX graph to cytoscape-compliant format and set the styling for cytoscape."""

import networkx as nx
import py4cytoscape as p4c
from py4cytoscape.networks import create_network_from_networkx

import pyBiodatafuse.constants as Cons


def _replace_graph_attrs(g: nx.MultiDiGraph):
    """Adapt the node and edge attributes keys to the cytoscape json structure.

    :param g: input NetworkX graph object.
    :returns: output NetworkX graph object.
    :raises ValueError: if a node is missing a name or an edge is missing a label.
    """
    # Fix node attributes
    for node, attrs in g.nodes(data=True):
        # Set 'id' explicitly
        attrs[Cons.ID] = node

        # Rename label to node type
        attrs[Cons.NODE_TYPE] = attrs[Cons.LABEL]
        del attrs[Cons.LABEL]

        # Adding a label i.e name of the node
        if Cons.NAME in attrs:
            attrs[Cons.LABEL] = attrs[Cons.NAME]
            del attrs[Cons.NAME]
        else:
            raise ValueError("Missing name ", attrs)

    # Fix edge attributes
    for u, v, k in g.edges(keys=True):
        attrs = g[u][v][k]

        # Edge source/target must match node IDs
        attrs[Cons.SOURCE] = u
        attrs[Cons.TARGET] = v
        if Cons.LABEL in attrs:
            attrs[Cons.INTERACTION] = attrs[Cons.LABEL]
        else:
            raise ValueError("Missing label ", attrs)

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
    g = g.copy()
    adj_g = _replace_graph_attrs(g)

    # Define the visual style as a dictionary
    default = {
        "title": "BioDataFuse_style",
        "layout_name": "force-directed",
        "defaults": [
            {"visualProperty": "NODE_FILL_COLOR", "value": "#FF0000"},
            {"visualProperty": "EDGE_COLOR", "value": "#000000"},
        ],
        "mappings": [],
    }

    # Create the network in Cytoscape
    create_network_from_networkx(
        adj_g,
        title=network_name,
        collection="BioDataFuse",
    )

    # Apply the visual style
    p4c.styles.create_visual_style(default)

    # Define node shape and color mapping
    column = Cons.LABEL
    values = list(Cons.ALL_NODE_LABELS.keys())
    shapes = list(Cons.ALL_NODE_LABELS.values())
    colors = list(Cons.COLOR_MAPPER.values())

    # Apply node shape and color mappings
    p4c.set_node_color_mapping(
        column,
        values,
        colors,
        mapping_type="d",
        style_name="default",
    )

    p4c.set_node_shape_mapping(column, values, shapes, style_name="default")
