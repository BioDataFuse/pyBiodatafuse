# coding: utf-8

"""
Gene and protein node generation for BDF RDF graphs.

This module creates gene and protein nodes from BioDataFuse data.
"""

from typing import Optional

from rdflib import Graph, URIRef

import pyBiodatafuse.constants as Cons
from pyBiodatafuse.graph.rdf.nodes.base import add_label, add_type, create_node


def get_gene_node(g: Graph, row) -> Optional[URIRef]:
    """
    Create and add a gene node to the RDF graph.

    Dynamically constructs the gene node URI based on the target source.
    Supports multiple identifier sources (Ensembl, Entrez Gene, etc.)
    as defined in NODE_URI_PREFIXES.

    :param g: RDF graph.
    :param row: Data row containing gene information.
    :return: URIRef for the gene node, or None if invalid data.
    """
    target = row.get(Cons.TARGET_COL)
    target_source = row.get(Cons.TARGET_SOURCE_COL)

    if not target or not target_source:
        return None

    # Look up the URI prefix for this source
    uri_prefix = Cons.NODE_URI_PREFIXES.get(target_source)
    if not uri_prefix:
        return None

    gene_node = create_node(uri_prefix, target)
    add_type(g, gene_node, Cons.NODE_TYPES["gene_node"])
    add_label(g, gene_node, target)

    return gene_node
