# protein_protein.py

"""Populate a BDF RDF graph with PPI nodes."""

from typing import Optional

from rdflib import Graph, Literal, URIRef
from rdflib.namespace import RDF, XSD

import pyBiodatafuse.constants as Cons


def add_ppi_data(
    g: Graph, gene_node: URIRef, entry: dict, base_uri: str, new_uris: dict
) -> Optional[URIRef]:
    """Add a protein protein interaction node.

    :param g: RDFLib graph
    :param gene_node: URIRef for the target protein
    :param entry: the ppi dictionary
    :param base_uri: the base URI for the project
    :param new_uris: dictionary with project node URIs
    :return: a ppi URIRef node
    """
    gene_link = entry.get(Cons.STRING_PPI_LINK_TO, None)
    gene_link_node = URIRef(Cons.BASE_URLS_DBS["uniprot"] + gene_link)  # type: ignore
    # ensembl = entry.get("Ensembl", None).split(":")[1]
    score = entry.get(Cons.STRING_PPI_SCORE, None)
    uniprot = entry.get(Cons.UNIPROT_TREMBL_A, None)
    uniprot_link = entry.get(Cons.UNIPROT_TREMBL_B, None)
    if score:
        score = float(score)
        # Nodes
        ppi_node = URIRef(base_uri + f"ppi/{uniprot}_{gene_link_node}")
        # ensembl_node = URIRef(BASE_URLS_DBS["ensembl"] + f"{ensembl}")
        if not uniprot:
            return None

        protein_node = URIRef(Cons.BASE_URLS_DBS["uniprot"] + uniprot)  # type: ignore
        protein_link_node = URIRef(Cons.BASE_URLS_DBS["uniprot"] + uniprot_link)  # type: ignore
        score_node = URIRef(f"{new_uris['score_base_node']}/{uniprot}_{uniprot_link}")

        # Edges
        new_edges_to_add = [
            (ppi_node, URIRef(Cons.PREDICATES["sio_has_part"]), protein_node),
            (ppi_node, URIRef(Cons.PREDICATES["sio_has_part"]), protein_link_node),
            (gene_node, URIRef(Cons.PREDICATES["translates_to"]), protein_node),
            (gene_link_node, URIRef(Cons.PREDICATES["translates_to"]), protein_link_node),
            (protein_node, URIRef(Cons.PREDICATES["translation_of"]), gene_node),
            (protein_link_node, URIRef(Cons.PREDICATES["translation_of"]), gene_link_node),
            (ppi_node, RDF.type, URIRef(Cons.NODE_TYPES["ppi_node"])),
            (protein_link_node, RDF.type, URIRef(Cons.NODE_TYPES["protein_node"])),
            (protein_node, RDF.type, URIRef(Cons.NODE_TYPES["protein_node"])),
            (score_node, RDF.type, URIRef(Cons.NODE_TYPES["score_node"])),
            (
                score_node,
                URIRef(Cons.PREDICATES["sio_has_value"]),
                Literal(score, datatype=XSD.double),
            ),
            (ppi_node, URIRef(Cons.PREDICATES["sio_has_measurement_value"]), score_node),
        ]

        for edge in new_edges_to_add:
            g.add(edge)

        return ppi_node

    return None
