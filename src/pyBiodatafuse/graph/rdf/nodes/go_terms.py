# coding: utf-8

"""
Gene Ontology term node generation for BDF RDF graphs.

This module creates GO term nodes from OpenTargets annotation data.
"""

from typing import Optional

from bioregistry import get_iri
from rdflib import Graph, URIRef
from rdflib.namespace import RDFS

import pyBiodatafuse.constants as Cons
from pyBiodatafuse.graph.rdf.nodes.base import add_label, add_triple, link_has_source, safe_get
from pyBiodatafuse.graph.rdf.utils import add_data_source_node


def add_go_cpf(g: Graph, process_data: dict) -> Optional[URIRef]:
    """
    Create and add a Gene Ontology (GO) node to the RDF graph.

    :param g: RDF graph.
    :param process_data: Dictionary containing GO information.
    :return: URIRef for the GO node, or None if invalid data.
    """
    curie = safe_get(process_data, Cons.GO_ID)
    if not curie:
        return None

    iri = get_iri(curie)
    if not iri:
        return None

    go_node = URIRef(iri)

    # Add label
    label = safe_get(process_data, Cons.GO_NAME, curie)
    add_label(g, go_node, label)

    # Add GO type as subclass using add_triple
    go_type = safe_get(process_data, Cons.GO_TYPE)
    if go_type and go_type in Cons.GO_TYPES:
        add_triple(g, go_node, RDFS.subClassOf, URIRef(Cons.GO_TYPES[go_type]))

    # Add data source
    source_node = add_data_source_node(g, Cons.OPENTARGETS_REACTOME_COL)
    link_has_source(g, go_node, source_node)

    return go_node

