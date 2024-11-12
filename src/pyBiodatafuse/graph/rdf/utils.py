import logging

import numpy as np
import pandas as pd
from bioregistry import get_iri, normalize_curie
from rdflib import Graph, Literal, URIRef
from rdflib.namespace import OWL, RDF, RDFS, XSD

from pyBiodatafuse.constants import DATA_SOURCES, NODE_TYPES

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


# Helper for NA replacement
def replace_na_none(item):
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


# Safe CURIE extraction
def extract_curie(prefix, identifier):
    curie = normalize_curie(f"{prefix}:{identifier}")
    if not curie:
        logger.warning(f"Could not normalize CURIE for {prefix}:{identifier}")
    return curie


# Generic URI constructor
def construct_uri(base_uri, identifier):
    return URIRef(f"{base_uri}/{identifier}")


# Add data source node
def add_data_source_node(g: Graph, source: str) -> URIRef:
    """Create and add a data source node to the RDF graph.

    :param g: RDF graph to which the data source node will be added.
    :param source: String containing the name of the source of the data
    :return: URIRef for the created data source node.
    """
    data_source_name = Literal(source, datatype=XSD.string)
    data_source_url = URIRef(DATA_SOURCES[source])
    g.add((data_source_url, RDF.type, URIRef(NODE_TYPES["data_source_node"])))
    g.add((data_source_url, RDFS.label, data_source_name))
    return data_source_url
