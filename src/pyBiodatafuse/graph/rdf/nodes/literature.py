# literature.py


"""Populate a BDF RDF graph with literature-based evidence."""


from rdflib import Graph, Literal, URIRef
from rdflib.namespace import RDF, RDFS

from pyBiodatafuse.constants import NAMESPACE_BINDINGS, NODE_TYPES, PREDICATES


def add_literature_based_data(g: Graph, entry: dict, gene_node: URIRef) -> None:
    """Add literature-based data for gene associations.

    :param g: (Graph): RDF graph to which the literature-based data is added.
    :param entry: (dict): Dictionary with literature-based association information.
    :param gene_node: (URIRef): URIRef of the gene node associated with the literature data.
    """
    source = entry.get("source", None)
    if source and "PMID" in source:
        source_id = source.split(": ")[1]
        source_url = f"https://pubmed.ncbi.nlm.nih.gov/{source_id}"
        article_node = URIRef(source_url)
        umls = entry.get("UMLS", None)
        mondo = entry.get("MONDO", None)
        disease_name = entry["disease_name"]
        # Create the source node and add metadata
        g.add((article_node, RDF.type, URIRef(NODE_TYPES["article"])))
        g.add((article_node, URIRef(PREDICATES["sio_refers_to"]), gene_node))
        if umls:
            disease_node = URIRef(f"{NAMESPACE_BINDINGS['umls']}{umls}")
            g.add(
                (
                    article_node,
                    URIRef(PREDICATES["sio_refers_to"]),
                    disease_node,
                )
            )
            g.add((disease_node, RDFS.label, Literal(disease_name)))
            g.add((disease_node, RDF.type, URIRef(NODE_TYPES["disease_node"])))
        if mondo:
            disease_node = URIRef(f"{NAMESPACE_BINDINGS['mondo']}{mondo}")
            g.add(
                (
                    article_node,
                    URIRef(PREDICATES["sio_refers_to"]),
                    disease_node,
                )
            )
            g.add((disease_node, RDFS.label, Literal(disease_name)))
            g.add((disease_node, RDF.type, URIRef(NODE_TYPES["disease_node"])))
