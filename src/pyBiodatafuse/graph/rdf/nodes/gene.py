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

    :param g: RDF graph.
    :param row: Data row containing gene information.
    :return: URIRef for the gene node, or None if invalid data.
    """
    target = row.get(Cons.TARGET_COL)
    source = row.get(Cons.TARGET_SOURCE_COL)

    if source != Cons.ENSEMBL or not target:
        return None

    gene_node = create_node(Cons.NODE_URI_PREFIXES[Cons.ENSEMBL], target)
    add_type(g, gene_node, Cons.NODE_TYPES["gene_node"])
    add_label(g, gene_node, row.get(Cons.IDENTIFIER_COL, target))

    return gene_node

