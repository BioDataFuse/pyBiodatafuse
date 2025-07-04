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
    for node in list(g.nodes()):
        for k, v in list(g.nodes[node].items()):
            if k == "name":
                nx.set_node_attributes(g, {node: {"label": v}})
                del g.nodes[node]["name"]
            elif k == "source":
                nx.set_node_attributes(g, {node: {"datasource": v}})
                del g.nodes[node]["source"]

            if k == "UBERON" and "id" in g.nodes[node]:
                current_id = g.nodes[node]["id"]
                if isinstance(current_id, str) and current_id.startswith("UBERON_"):
                    corrected_id = current_id.replace("UBERON_", "UBERON:")
                    nx.set_node_attributes(g, {node: {"id": corrected_id}})
            
            if k == "UBERON": 
                if "value" in g.nodes[node]:
                    current_value = g.nodes[node]["value"]
                    if isinstance(current_value, str) and current_value.startswith("UBERON_"):
                        g.nodes[node]["value"] = current_value.replace("UBERON_", "UBERON:")
                if "name" in g.nodes[node]:
                    current_name = g.nodes[node]["name"]
                    if isinstance(current_name, str) and current_name.startswith("UBERON_"):
                        g.nodes[node]["name"] = current_name.replace("UBERON_", "UBERON:")
            
    for u, v, k in list(g.edges(keys=True)):
        source_node_id = g.nodes[u].get("id")
        target_node_id = g.nodes[v].get("id")

        if g.nodes[u].get("datasource") == "Bgee" and g[u][v][k]["source"] == source_node_id.replace(":", "_"):
            g[u][v][k]["source"] = source_node_id
        if g.nodes[v].get("datasource") == "Bgee" and g[u][v][k]["target"] == target_node_id.replace(":", "_"):
            g[u][v][k]["target"] = target_node_id
        
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