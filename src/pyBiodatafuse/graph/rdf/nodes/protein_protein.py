# coding: utf-8

"""
Protein-protein interaction node generation for BDF RDF graphs.

This module creates PPI nodes from StringDB interaction data.
"""

from typing import Optional

from rdflib import Graph, Literal, URIRef
from rdflib.namespace import XSD

import pyBiodatafuse.constants as Cons
from pyBiodatafuse.graph.rdf.nodes.base import add_triples, add_type, create_node, safe_get


def add_ppi_data(g: Graph, gene_node: URIRef, entry: dict) -> Optional[URIRef]:
    """
    Add a protein-protein interaction node to the RDF graph.

    :param g: RDF graph.
    :param gene_node: URIRef for the source gene.
    :param entry: Dictionary containing PPI data from StringDB.
    :return: URIRef for the PPI node, or None if invalid data.
    """
    # Extract required data
    score = safe_get(entry, Cons.STRING_PPI_SCORE)
    uniprot_a = safe_get(entry, Cons.UNIPROT_TREMBL_A)
    uniprot_b = safe_get(entry, Cons.UNIPROT_TREMBL_B)
    gene_link = safe_get(entry, Cons.STRING_PPI_LINK_TO)

    if not score or not uniprot_a:
        return None

    score = float(score)

    # Create nodes
    uniprot_base = Cons.BASE_URLS_DBS["uniprot"]
    protein_a = create_node(uniprot_base, uniprot_a)
    protein_b = create_node(uniprot_base, uniprot_b) if uniprot_b else None
    gene_link_node = create_node(uniprot_base, gene_link) if gene_link else None

    ppi_node = URIRef(f"{g.base_uri}ppi/{uniprot_a}_{uniprot_b or gene_link}")
    score_node = URIRef(f"{g.new_uris['score_base_node']}/{uniprot_a}_{uniprot_b or gene_link}")

    # Add types
    add_type(g, ppi_node, Cons.NODE_TYPES["ppi_node"])
    add_type(g, protein_a, Cons.NODE_TYPES["protein_node"])
    add_type(g, score_node, Cons.NODE_TYPES["score_node"])
    if protein_b:
        add_type(g, protein_b, Cons.NODE_TYPES["protein_node"])

    # Build relationships
    triples = [
        # PPI structure
        (ppi_node, URIRef(Cons.PREDICATES["sio_has_part"]), protein_a),
        (ppi_node, URIRef(Cons.PREDICATES["sio_has_measurement_value"]), score_node),
        # Score value
        (score_node, URIRef(Cons.PREDICATES["sio_has_value"]), Literal(score, datatype=XSD.double)),
        # Gene-protein relationships
        (gene_node, URIRef(Cons.PREDICATES["translates_to"]), protein_a),
        (protein_a, URIRef(Cons.PREDICATES["translation_of"]), gene_node),
    ]

    if protein_b:
        triples.append((ppi_node, URIRef(Cons.PREDICATES["sio_has_part"]), protein_b))

    if gene_link_node and protein_b:
        triples.extend(
            [
                (gene_link_node, URIRef(Cons.PREDICATES["translates_to"]), protein_b),
                (protein_b, URIRef(Cons.PREDICATES["translation_of"]), gene_link_node),
            ]
        )

    add_triples(g, triples)

    return ppi_node
