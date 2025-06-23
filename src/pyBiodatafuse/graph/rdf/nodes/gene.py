# gene.py


"""Populate a BDF RDF graph with gene and protein nodes."""


from typing import Optional

from rdflib import Graph, Literal, URIRef
from rdflib.namespace import RDF, RDFS, XSD

import pyBiodatafuse.constants as Cons


def get_gene_node(g: Graph, row) -> tuple:
    """Create and add a gene node and associated protein node to the RDF graph.

    :param g: (Graph): RDF graph to which the gene and protein nodes are added.
    :param row: (pd.Series): Data row containing gene information.

    :return: URIRef for the gene node.
    """
    target = row.get(Cons.TARGET_COL, None)
    source = row.get(Cons.TARGET_SOURCE_COL, None)
    if source and source == Cons.ENSEMBL:
        if target:
            gene_node = URIRef(Cons.NODE_URI_PREFIXES[Cons.ENSEMBL] + target)
            g.add((gene_node, RDF.type, URIRef(Cons.NODE_TYPES["gene_node"])))
            g.add((gene_node, RDFS.label, Literal(row[Cons.IDENTIFIER_COL], datatype=XSD.string)))
            return gene_node
    return (None, None)
