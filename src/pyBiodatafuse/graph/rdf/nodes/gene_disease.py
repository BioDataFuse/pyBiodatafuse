"""Populate a BDF RDF graph with gene-disease relationship data."""

from typing import Optional

from bioregistry import get_iri
from rdflib import Graph, URIRef
from rdflib.namespace import RDFS, XSD

import pyBiodatafuse.constants as Cons
from pyBiodatafuse.graph.rdf.nodes.base import (
    add_label,
    add_same_as,
    add_triple,
    add_type,
    add_value,
    create_node,
    link_refers_to,
    safe_get,
)


def add_gene_disease_associations(
    g: Graph,
    id_number: str,
    source_idx: str,
    gene_node: URIRef,
    disease_data: dict,
    new_uris: dict,
    i: int,
) -> Optional[URIRef]:
    """Process and add gene-disease association data to the RDF graph.

    :param g: RDF graph to which the associations will be added.
    :param id_number: Unique identifier for the associations.
    :param source_idx: Source index for the associations.
    :param gene_node: URIRef of the gene node associated with the disease data.
    :param disease_data: Dictionary containing disease association information.
    :param new_uris: Dictionary with updated project base URIs for the nodes.
    :param i: The index of the row.
    :return: disease node (URIRef).
    """
    gene_id = g.value(subject=gene_node, predicate=RDFS.label)
    disease_node = add_disease_node(g, disease_data)
    if not disease_node:
        return None

    umlscui = safe_get(disease_data, Cons.UMLS, "")
    assoc_node = create_node(
        f"{new_uris['gene_disease_association']}/", f"{gene_id}{umlscui}"
    )
    add_type(g, assoc_node, Cons.NODE_TYPES["gene_disease_association"])
    link_refers_to(g, assoc_node, gene_node)
    link_refers_to(g, assoc_node, disease_node)

    score = safe_get(disease_data, Cons.DISGENET_SCORE, "")
    if score:
        score_node = add_score_node(
            g, id_number, source_idx, umlscui, score, new_uris, i, str(gene_id)
        )
        add_triple(g, assoc_node, "sio_has_measurement_value", score_node)

    if safe_get(disease_data, Cons.DISGENET_EI):
        ei_node = add_evidence_node(
            g, id_number, source_idx, Cons.DISGENET_EI, disease_data, new_uris, i
        )
        if ei_node:
            add_triple(g, assoc_node, "sio_has_measurement_value", ei_node)

    if safe_get(disease_data, Cons.DISGENET_EL):
        el_node = add_evidence_node(
            g, id_number, source_idx, Cons.DISGENET_EL, disease_data, new_uris, i
        )
        if el_node:
            add_triple(g, assoc_node, "sio_has_measurement_value", el_node)

    return disease_node


def add_disease_node(g: Graph, disease_data: dict) -> Optional[URIRef]:
    """Create and add a disease node to the RDF graph.

    :param g: RDF graph to which the disease node will be added.
    :param disease_data: Dictionary containing disease information.
    :return: URIRef for the created disease node.
    """
    disease_curie = safe_get(disease_data, Cons.UMLS)
    if not disease_curie:
        return None

    disease_node = create_node(Cons.NODE_URI_PREFIXES["medgen"], disease_curie)
    add_type(g, disease_node, Cons.NODE_TYPES["disease_node"])
    add_label(g, disease_node, safe_get(disease_data, Cons.DISEASE_NAME))

    for identifier in Cons.DISEASE_IDENTIFIER_TYPES:
        if safe_get(disease_data, identifier):
            curies = safe_get(disease_data, identifier, "").split(", ")
            for curie in curies:
                source_iri = get_iri(curie) or get_iri(f"obo:{curie.split(':')[-1]}")
                if source_iri:
                    add_same_as(g, disease_node, URIRef(source_iri), bidirectional=False)

    return disease_node


def add_score_node(
    g: Graph,
    id_number: str,
    source_idx: str,
    disease_id: str,
    score: float,
    new_uris: dict,
    i: int,
    gene_id: str,
) -> URIRef:
    """Create and add a score node for gene-disease associations.

    :param g: RDF graph to which the score node will be added.
    :param id_number: Unique identifier for the score node.
    :param source_idx: Source index for the association.
    :param disease_id: Disease identifier associated with the score.
    :param score: Score value for the gene-disease association.
    :param new_uris: Dictionary with updated project base URIs for the nodes.
    :param i: Index or iteration number used for node URI construction.
    :param gene_id: String value of the gene ID.
    :return: URIRef for the created score node.
    """
    score_node = create_node(
        f"{new_uris['score_base_node']}/",
        f"{id_number}{i}{source_idx}_{disease_id}{gene_id}"
    )
    add_type(g, score_node, Cons.NODE_TYPES["score_node"])
    add_value(g, score_node, score, XSD.double)

    return score_node


def add_evidence_node(
    g: Graph,
    id_number: str,
    source_idx: str,
    evidence_type: str,
    disease_data: dict,
    new_uris: dict,
    i: int,
) -> Optional[URIRef]:
    """Create and add an evidence node (EI or EL) to the RDF graph.

    :param g: RDF graph to which the evidence node will be added.
    :param id_number: Unique identifier for the evidence node.
    :param source_idx: Source index for the association.
    :param evidence_type: Type of evidence ("ei" or "el").
    :param disease_data: Dictionary containing disease information.
    :param new_uris: Dictionary with updated project base URIs for the nodes.
    :param i: Index or iteration number used for node URI construction.
    :return: URIRef for the created evidence node.
    """
    value = safe_get(disease_data, evidence_type)
    if value is None:
        return None

    umls_id = safe_get(disease_data, Cons.UMLS, "")
    node = create_node(
        f"{new_uris['score_base_node']}/{evidence_type}/",
        f"{id_number}{i}{source_idx}_{umls_id}"
    )

    datatype = XSD.double if evidence_type != "el" else XSD.string
    add_type(g, node, Cons.NODE_TYPES[f"{evidence_type}_node"])
    add_value(g, node, value, datatype)

    return node
