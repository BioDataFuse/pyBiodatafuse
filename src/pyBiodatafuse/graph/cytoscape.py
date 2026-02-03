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

    # Create the network in Cytoscape
    create_network_from_networkx(
        adj_g,
        title=network_name,
        collection="BioDataFuse",
    )

    # Define the visual style
    style_name = "BioDataFuse_style"

    # Check if style already exists, if so delete it first
    existing_styles = p4c.styles.get_visual_style_names()
    if style_name in existing_styles:
        p4c.styles.delete_visual_style(style_name)

    # Create the visual style with defaults
    defaults = {
        "NODE_FILL_COLOR": "#808080",
        "EDGE_TARGET_ARROW_SHAPE": "ARROW",
    }
    p4c.styles.create_visual_style(style_name, defaults=defaults)

    # Define node shape and color mapping
    column = Cons.NODE_TYPE
    values = list(Cons.ALL_NODE_LABELS.keys())
    shapes = list(Cons.ALL_NODE_LABELS.values())
    colors = list(Cons.COLOR_MAPPER.values())

    # Apply node color mapping
    p4c.set_node_color_mapping(
        column,
        values,
        colors,
        mapping_type="d",
        style_name=style_name,
    )

    # Apply node shape mapping
    p4c.set_node_shape_mapping(column, values, shapes, style_name=style_name)

    # Set the node label to show the "labels" attribute (which contains the name)
    # The column name in Cytoscape is "labels" (from Cons.LABEL)
    p4c.set_node_label_mapping(Cons.LABEL, style_name=style_name)

    # Apply the visual style to the network
    p4c.styles.set_visual_style(style_name)
