# gene.py


"""Populate a BDF RDF graph with gene and protein nodes."""


from rdflib import Graph, Literal, URIRef
from rdflib.namespace import RDF, RDFS, XSD

from pyBiodatafuse.constants import NODE_TYPES


def add_gene_nodes(g: Graph, row) -> URIRef:
    """Create and add a gene node to the RDF graph.

    :param g: RDF graph to which the gene nodes are added.
    :param row: DataFrame row containing gene data.
    :return: URIRef for the gene node.
    """
    source = row.get("source")
    target = row.get("target")
    if source and source == "Ensembl":
        gene_node = URIRef(f"http://identifiers.org/ensembl#{target}")
        g.add((gene_node, RDF.type, URIRef(NODE_TYPES["gene_node"])))
        g.add((gene_node, RDFS.label, Literal(row["identifier"], datatype=XSD.string)))
    # if source and source == "HGNC":
    return gene_node
