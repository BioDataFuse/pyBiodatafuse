"""Provide utils for BDF RDF."""

import logging
import os

import numpy as np
import pandas as pd
from bioregistry import normalize_curie
from rdflib import BNode, Graph, Literal, URIRef
from rdflib.namespace import RDF, RDFS, SH, XSD
from shexer.consts import SHACL_TURTLE, TURTLE
from shexer.shaper import Shaper

from pyBiodatafuse.constants import DATA_SOURCES, NAMESPACE_BINDINGS, NAMESPACE_SHAPES, NODE_TYPES

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def replace_na_none(item):
    """Replace occurrences of NA values (such as 'na', 'nan', 'none') with None.

    :param item: Item to process. Can be a string, float, list, dict, or numpy array.
    :return: Processed item with NA values replaced by None.
    """
    if isinstance(item, str) and item.lower() in ["na", "nan", "none"]:
        return None
    elif item is None or (isinstance(item, float) and pd.isna(item)):
        return None
    elif isinstance(item, list):
        return [replace_na_none(sub_item) for sub_item in item]
    elif isinstance(item, dict):
        return {key: replace_na_none(value) for key, value in item.items()}
    elif isinstance(item, np.ndarray):
        return np.array([replace_na_none(sub_item) for sub_item in item], dtype=object)
    return item


def extract_curie(prefix, identifier):
    """Generate a CURIE by normalizing a prefix and identifier.

    :param prefix: Prefix string, such as a registry identifier.
    :param identifier: Identifier to be appended to the prefix.
    :return: Normalized CURIE or None if normalization fails.
    """
    curie = normalize_curie(f"{prefix}:{identifier}")
    if not curie:
        logger.warning("Could not normalize CURIE for %s:%s", prefix, identifier)
    return curie


def construct_uri(base_uri, identifier):
    """Construct a URIRef from a base URI and an identifier.

    :param base_uri: Base URI string for the RDF resource.
    :param identifier: Identifier to append to the base URI.
    :return: A URIRef representing the constructed URI.
    """
    return URIRef(f"{base_uri}/{identifier}")


def add_data_source_node(g: Graph, source: str) -> URIRef:
    """Create and add a data source node to the RDF graph.

    :param g: RDF graph to which the data source node will be added.
    :param source: String containing the name of the source of the data.
    :return: URIRef for the created data source node.
    """
    # TODO: fix
    if "OpenTargets" in source:
        source = "OpenTargets"
    data_source_name = Literal(source, datatype=XSD.string)
    data_source_url = URIRef(DATA_SOURCES[source])
    g.add((data_source_url, RDF.type, URIRef(NODE_TYPES["data_source_node"])))
    g.add((data_source_url, RDFS.label, data_source_name))
    return data_source_url


def get_shapes(
    g,
    base_uri,
    path,
    threshold=0.000000001,
    graph_type="shex",  # "shex" or "shacl"
    uml_figure_path=None,
    print_string_output=True,
    additional_namespaces=None,
):
    """Use shexer (https://github.com/DaniFdezAlvarez/shexer) on the BDF graph to generate Shex or SHACL.

    :param g: RDF graph to generate shapes from.
    :param base_uri: The graph iri to be added to shaper namespaces.
    :param path: relative path in which the graph TTL will be saved, if provided.
    :param threshold: float between [0,1] used to accept shapes based on frequency.
    :param graph_type: "shex" or "shacl", to specify which graph type to generate.
    :param uml_figure_path: str path where the generated UML is stored.
    :param print_string_output: bool, print or not the generated TTL as a string.
    :param additional_namespaces: dictionary containing {namespace: prefix} pairs.
    :raises ValueError: If the graph type is not a valid string or not in ['shex', 'shacl'].
    :return: shaper shex or shacl graph
    """
    # Graph type: shex or shacl
    graph_type = graph_type.lower()
    if graph_type not in ["shex", "shacl"]:
        raise ValueError("Invalid graph_type specified. Choose 'shex' or 'shacl'.")

    try:
        graph_type = graph_type.lower()
    except AttributeError as exc:
        raise ValueError("graph_type must be a string.") from exc

    # Default namespaces: NAMESPACE_SHAPES
    NAMESPACE_SHAPES[base_uri] = "graph"
    # Merge with additional namespaces if provided
    if additional_namespaces:
        NAMESPACE_SHAPES.update(additional_namespaces)
    # Initialize Shaper with the given graph and namespaces
    shaper = Shaper(
        all_classes_mode=True,
        rdflib_graph=g,
        # input_format=TURTLE,
        namespaces_dict=NAMESPACE_SHAPES,
        # disable_or_statements=False, Workaround for bug
    )
    graph_result = None
    # Generate the appropriate graph (Shex or SHACL)
    rdf_png_path = os.path.join(os.getcwd(), uml_figure_path) if uml_figure_path else None
    if graph_type == "shex":
        graph_result = shaper.shex_graph(
            string_output=True, acceptance_threshold=threshold, to_uml_path=rdf_png_path
        )
    elif graph_type == "shacl":
        graph_result = shaper.shex_graph(
            string_output=True,
            acceptance_threshold=threshold,
            to_uml_path=rdf_png_path,
            output_format=SHACL_TURTLE,
        )

    # Save the output to a file if path is provided
    if path and graph_result:
        # Ensure the directory exists before saving the file
        dir_path = os.path.dirname(path)
        if dir_path and not os.path.exists(dir_path):
            os.makedirs(dir_path)  # Create the directory if it doesn't exist
        with open(path, "a", encoding="utf-8") as f:
            f.write(graph_result)

    # Optionally print the graph to the console if print_string_output is True
    if print_string_output:
        print(graph_result)

    return graph_result


def get_shacl_prefixes(namespaces, path, new_uris):
    """
    Generate SHACL prefix declarations and save them in Turtle format.

    :param namespaces: Optional dictionary of prefix to namespace URI mappings to include in the SHACL declarations.
    :param path: Optional path to a file where the Turtle data will be written. If not provided, the data is not written to disk.
    :param new_uris: Dictionary of prefix to namespace URI mappings to include in the SHACL declarations.

    :return: A RDFLib Graph containing the SHACL prefix declarations.
    """
    graph = Graph()

    def add_declarations(prefix_dict):
        """
        Add declarations for an dictionary of namespaces.

        :param prefix_dict: Dictionary of  {prefix:namespace,}.
        """
        for prefix, ns_uri in prefix_dict.items():
            declare_node = BNode()
            graph.add((declare_node, SH.prefix, Literal(prefix)))
            graph.add((declare_node, SH.namespace, Literal(ns_uri, datatype=XSD.anyURI)))
            graph.add((BNode(), SH.declare, declare_node))

    add_declarations(new_uris)
    add_declarations(NAMESPACE_BINDINGS)
    if namespaces:
        add_declarations(namespaces)

    ttl_data = graph.serialize(format="ttl")

    if path:
        try:
            with open(path, "w", encoding="UTF-8") as f:
                f.write(ttl_data)
        except IOError as e:
            logger.error("Error writing to file %s: %s", path, e)

    print(ttl_data)
    return graph


def get_node_label(g, node):
    """
    Retrieve the label of a given node from an RDF graph.

    :param g: The RDF graph containing the data.
    :param node: The node whose label is to be retrieved.

    :return: The label of the node if it exists, otherwise None.
    """
    for stmt in g.triples((node, RDFS.label, None)):
        return stmt[2]
    return None
