# gene_protein.py


"""Populate a BDF RDF graph with gene and protein nodes."""


from rdflib import Graph, Literal, URIRef
from rdflib.namespace import RDF, RDFS, XSD

from pyBiodatafuse.constants import NODE_TYPES, PREDICATES


def add_gene_protein_nodes(g: Graph, row) -> tuple:
    """Create and add a gene node and associated protein node to the RDF graph.

    :param g: (Graph): RDF graph to which the gene and protein nodes are added.
    :param row: (pd.Series): Data row containing gene information.

    :return: tuple containing URIRefs for the gene node and a list of URIRefs for protein nodes respectively.
    """
    target = row.get("target", None)
    if target:
        protein_nodes = []
        uniprot = row.get("Uniprot-TrEMBL", [])
        gene_node = URIRef(f"http://identifiers.org/ensembl#{target}")
        g.add((gene_node, RDF.type, URIRef(NODE_TYPES["gene_node"])))
        g.add((gene_node, RDFS.label, Literal(row["identifier"], datatype=XSD.string)))
        if uniprot:
            for prot in uniprot:
                protein_node = URIRef(f"https://www.uniprot.org/uniprot/{prot}")
                prot_label = prot
                # Link gene to protein
                g.add((protein_node, URIRef(PREDICATES["translation_of"]), gene_node))
                g.add((gene_node, URIRef(PREDICATES["translates_to"]), protein_node))

                # Add protein node to graph
                g.add((protein_node, RDF.type, URIRef(NODE_TYPES["protein_node"])))
                g.add((protein_node, RDFS.label, Literal(prot_label, datatype=XSD.string)))
                protein_nodes.append(protein_node)
        return gene_node, protein_nodes
    return None, None
