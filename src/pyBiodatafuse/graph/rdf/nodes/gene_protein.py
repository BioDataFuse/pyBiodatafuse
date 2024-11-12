from bioregistry import get_iri, normalize_curie
from rdflib import Graph, Literal, URIRef
from rdflib.namespace import RDF, RDFS, XSD

from pyBiodatafuse.constants import NODE_TYPES, PREDICATES
from pyBiodatafuse.graph.rdf.utils import add_data_source_node


def add_gene_protein_nodes(g: Graph, row) -> tuple:
    """Creates and adds a gene node and associated protein node to the RDF graph.

    Args:
        g (Graph): RDF graph to which the gene and protein nodes are added.
        row (pd.Series): Data row containing gene information.

    Returns:
        tuple: A tuple containing URIRefs for the gene node and protein node, respectively.
    """
    target = row["target"]
    if target:
        gene_node = URIRef(f"http://identifiers.org/ensembl#{target}")
        protein_node = URIRef(f"http://identifiers.org/ensembl#{target}xProtein")

        # Add gene node to graph
        g.add((gene_node, RDF.type, URIRef(NODE_TYPES["gene_node"])))
        g.add((gene_node, RDFS.label, Literal(row["identifier"], datatype=XSD.string)))

        # Link gene to protein
        g.add((protein_node, URIRef(PREDICATES["translation_of"]), gene_node))
        g.add((gene_node, URIRef(PREDICATES["translates_to"]), protein_node))

        # Add protein node to graph
        g.add((protein_node, RDF.type, URIRef(NODE_TYPES["protein_node"])))
        g.add((protein_node, RDFS.label, Literal(f"{target}xProtein", datatype=XSD.string)))

        return gene_node, protein_node
    return None, None
