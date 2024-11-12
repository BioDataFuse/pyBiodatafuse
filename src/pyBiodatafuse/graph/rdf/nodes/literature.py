# literature.py


"""Populate a BDF RDF graph with literature-based evidence."""


from rdflib import Graph, Literal, URIRef
from rdflib.namespace import RDF, RDFS

from pyBiodatafuse.constants import NODE_TYPES, PREDICATES


def add_literature_based_data(g: Graph, entry: dict, gene_node: URIRef) -> None:
    """Add literature-based data for gene associations.

    :param g: (Graph): RDF graph to which the literature-based data is added.
    :param entry: (dict): Dictionary with literature-based association information.
    :param gene_node: (URIRef): URIRef of the gene node associated with the literature data.
    """
    source = entry.get("source", None)
    if source and "PMID" in source:
        source_id = source.split(": ")[1]
        identifier = entry["id"]
        disease_name = entry["disease_name"]
        source_url = f"https://pubmed.ncbi.nlm.nih.gov/{source_id}"

        # Create the source node and add metadata
        article_node = URIRef(source_url)
        g.add((article_node, RDF.type, URIRef(NODE_TYPES["article"])))
        g.add((article_node, URIRef(PREDICATES["sio_refers_to"]), gene_node))
        g.add(
            (
                article_node,
                URIRef(PREDICATES["sio_refers_to"]),
                URIRef(f"https://biodatafuse.org/identifiers/{identifier}"),
            )
        )

        # Add disease label
        disease_node = URIRef(f"https://biodatafuse.org/identifiers/{identifier}")
        g.add((disease_node, RDFS.label, Literal(disease_name)))
        g.add((disease_node, RDF.type, URIRef(NODE_TYPES["disease_node"])))
