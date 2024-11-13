# protein_protein.py

"""Populate a BDF RDF graph with PPI nodes."""

from rdflib import Graph, Literal, URIRef
from rdflib.namespace import OWL, RDF, RDFS, XSD

from pyBiodatafuse.constants import NODE_TYPES, PREDICATES, BASE_URLS_DBS


def add_ppi_data(  # TODO refactor using constants
    g: Graph, protein_node: URIRef, protein_name: str, entry: dict, base_uri: str, new_uris: dict
) -> URIRef:
    """Add a protein protein interaction node.

    :param g: RDFLib graph
    :param protein_node: URIRef for the target protein
    :param protein_name: the name of the target protein
    :param entry: the ppi dictionary
    :param base_uri: the base URI for the project
    :param new_uris: dictionary with project node URIs
    :return: a ppi URIRef node
    """
    stringdb_link_to = entry.get("stringdb_link_to", None)
    ensembl = entry.get("Ensembl", None).split(":")[1]
    score = entry.get("score", None)

    if score:
        score = float(score)
        # Nodes
        ppi_node = URIRef(base_uri + f"ppi/{protein_name}_{ensembl}")
        g.add(
            (
                ppi_node,
                RDF.type,
                URIRef(NODE_TYPES["ppi_node"]),
            )
        )
        g.add(
            (
                URIRef(BASE_URLS_DBS["uniprot"] + f"{stringdb_link_to}"),
                RDF.type,
                URIRef(NODE_TYPES["protein_node"]),
            )
        )
        g.add(
            (
                ppi_node,
                URIRef(PREDICATES["sio_has_part"]),
                URIRef(BASE_URLS_DBS["uniprot"] + f"{stringdb_link_to}"),
            )
        )
        g.add(
            (
                ppi_node,
                URIRef(PREDICATES["sio_has_part"]),
                protein_node,
            )
        )
        g.add(
            (
                URIRef(BASE_URLS_DBS["uniprot"] + f"{stringdb_link_to}"),
                URIRef(PREDICATES["sio_is_part_of"]),
                ppi_node,
            )
        )
        g.add(
            (
                URIRef(BASE_URLS_DBS["uniprot"] + f"{stringdb_link_to}"),
                RDF.type,
                URIRef(NODE_TYPES["protein_node"]),
            )
        )
        g.add(
            (
                URIRef(BASE_URLS_DBS["uniprot"] + f"{stringdb_link_to}"),
                RDFS.label,
                Literal(stringdb_link_to, datatype=XSD.string),
            )
        )
        g.add(
            (
                URIRef(BASE_URLS_DBS["ensembl"] + f"{ensembl}"),
                OWL.sameAs,
                URIRef(BASE_URLS_DBS["uniprot"] + f"{stringdb_link_to}"),
            )
        )
        g.add(
            (
                URIRef(BASE_URLS_DBS["ensembl"] + f"{ensembl}"),
                RDF.type,
                URIRef(NODE_TYPES["protein_node"]),
            )
        )
        score_node = URIRef(f"{new_uris['score_base_node']}/{protein_name}_{ensembl}")
        g.add((ppi_node, URIRef(PREDICATES["sio_has_measurement_value"]), score_node))
        g.add((score_node, RDF.type, URIRef(NODE_TYPES["score_node"])))
        g.add(
            (
                score_node,
                URIRef(PREDICATES["sio_has_value"]),
                Literal(score, datatype=XSD.double),
            )
        )
        return ppi_node
