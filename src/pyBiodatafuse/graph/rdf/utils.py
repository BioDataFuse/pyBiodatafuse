"""Provide utils for BDF RDF."""

import logging
import os

import rdflib
from shexer.consts import TURTLE, SHACL_TURTLE
from shexer.shaper import Shaper
import numpy as np
import pandas as pd
from bioregistry import normalize_curie
from rdflib import Graph, Literal, URIRef
from rdflib.namespace import RDF, RDFS, XSD

from pyBiodatafuse.constants import DATA_SOURCES, NODE_TYPES

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
    data_source_name = Literal(source, datatype=XSD.string)
    data_source_url = URIRef(DATA_SOURCES[source])
    g.add((data_source_url, RDF.type, URIRef(NODE_TYPES["data_source_node"])))
    g.add((data_source_url, RDFS.label, data_source_name))
    return data_source_url


def matching_triples(temp_graph, subj_set, pred_set, obj_set):
    """Match triples from a graph and a set of s, p, o from another graph.

    :param g: RDF graph to which the data source node will be added.
    :param subj_set: Set containing all s nodes.
    :param pred_set: Set containing all p nodes.
    :param obj_set: Set containing all o nodes.
    :yield: Tuples (subj, pred, obj) where at least one element is found in subj_set, pred_set, or obj_set.
    """
    for subj, pred, obj in temp_graph:
        if subj in subj_set or pred in pred_set or obj in obj_set:
            yield (subj, pred, obj)


def get_shapes(
    g,
    base_uri,
    path,
    threshold=0.001,
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
    :return: shaper shex or shacl graph
    """
    graph_type = graph_type.lower()
    # Default namespaces
    namespaces_dict = {
        "http://www.w3.org/1999/02/22-rdf-syntax-ns#": "rdf",
        "http://example.org/": "ex",
        "http://weso.es/shapes/": ":",
        "http://www.w3.org/2001/XMLSchema#": "xsd",
        "http://www.w3.org/2002/07/owl#": "owl",
        f"{base_uri}": "graph",
        "http://purl.obolibrary.org/obo/": "obo",
        "http://purl.obolibrary.org/obo/so#": "so",
    }

    # Merge with additional namespaces if provided
    if additional_namespaces:
        namespaces_dict.update(additional_namespaces)

    # Initialize Shaper with the given graph and namespaces
    shaper = Shaper(
        all_classes_mode=True,
        rdflib_graph=g,
        input_format=TURTLE,
        namespaces_dict=namespaces_dict,
    )

    # Generate the appropriate graph (Shex or SHACL)
    rdf_png_path = os.path.join(os.getcwd(), uml_figure_path) if uml_figure_path else None

    if graph_type == "shex":
        graph_result = shaper.shex_graph(
            string_output=True, acceptance_threshold=threshold, to_uml_path=rdf_png_path
        )
    elif graph_type == "shacl":
        graph_result = shaper.shex_graph(
            string_output=True, acceptance_threshold=threshold, to_uml_path=rdf_png_path,
            output_format=SHACL_TURTLE
        )
    else:
        raise ValueError("Invalid graph_type specified. Choose 'shex' or 'shacl'.")

    # Save the output to a file if path is provided
    if path:
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


# Implement function to generate prefix shacl (used in SPARQL endpoint) from the graph
