"""Populate a BDF RDF graph with gene-disease relationship data."""

from typing import Optional

from bioregistry import get_iri
from rdflib import Graph, Literal, URIRef
from rdflib.namespace import OWL, RDF, RDFS, XSD

import pyBiodatafuse.constants as Cons


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
    umlscui = disease_data.get(Cons.UMLS, "")
    assoc_node = URIRef(f"{new_uris['gene_disease_association']}/{gene_id}{umlscui}")
    g.add((assoc_node, RDF.type, URIRef(Cons.NODE_TYPES["gene_disease_association"])))
    g.add((assoc_node, URIRef(Cons.PREDICATES["sio_refers_to"]), gene_node))
    g.add((assoc_node, URIRef(Cons.PREDICATES["sio_refers_to"]), disease_node))
    score = disease_data.get(Cons.DISGENET_SCORE, "")
    if disease_data.get(Cons.DISGENET_SCORE):
        score_node = add_score_node(
            g,
            id_number,
            source_idx,
            umlscui,
            score,
            new_uris,
            i,
            str(gene_id),
        )
        g.add((assoc_node, URIRef(Cons.PREDICATES["sio_has_measurement_value"]), score_node))

    if disease_data.get(Cons.DISGENET_EI):
        ei_node = add_evidence_node(
            g, id_number, source_idx, Cons.DISGENET_EI, disease_data, new_uris, i
        )
        if ei_node:
            g.add((assoc_node, URIRef(Cons.PREDICATES["sio_has_measurement_value"]), ei_node))

    if disease_data.get(Cons.DISGENET_EL):
        el_node = add_evidence_node(
            g, id_number, source_idx, Cons.DISGENET_EL, disease_data, new_uris, i
        )
        if el_node:
            g.add((assoc_node, URIRef(Cons.PREDICATES["sio_has_measurement_value"]), el_node))

    # data_source_node = add_data_source_node(g, "DISGENET")
    # g.add((assoc_node, URIRef(PREDICATES["sio_has_source"]), data_source_node))
    return disease_node


def add_disease_node(g: Graph, disease_data: dict) -> Optional[URIRef]:
    """Create and add a disease node to the RDF graph.

    :param g: RDF graph to which the disease node will be added.
    :param disease_data: Dictionary containing disease information.
    :return: URIRef for the created disease node.
    """
    disease_curie = disease_data.get(Cons.UMLS)
    if not disease_curie:
        return None

    disease_iri = Cons.NODE_URI_PREFIXES["medgen"] + disease_curie
    disease_node = URIRef(disease_iri)
    g.add((disease_node, RDF.type, URIRef(Cons.NODE_TYPES["disease_node"])))
    g.add(
        (
            disease_node,
            RDFS.label,
            Literal(disease_data.get(Cons.DISEASE_NAME), datatype=XSD.string),
        )
    )

    for identifier in Cons.DISEASE_IDENTIFIER_TYPES:
        if disease_data.get(identifier):
            curies = disease_data.get(identifier, "").split(", ")
            for curie in curies:
                source_iri = get_iri(curie) or get_iri(f"obo:{curie.split(':')[-1]}")
                if source_iri:
                    g.add((disease_node, OWL.sameAs, URIRef(source_iri)))

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
    score_node = URIRef(
        f"{new_uris['score_base_node']}/{id_number}{i}{source_idx}_{disease_id}{gene_id}"
    )
    g.add((score_node, RDF.type, URIRef(Cons.NODE_TYPES["score_node"])))
    g.add(
        (score_node, URIRef(Cons.PREDICATES["sio_has_value"]), Literal(score, datatype=XSD.double))
    )
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
    value = disease_data.get(evidence_type)
    if value is None:
        return None

    node = URIRef(
        f"{new_uris['score_base_node']}/{evidence_type}/{id_number}{i}{source_idx}_{disease_data[Cons.UMLS]}"
    )
    datatype = XSD.double
    if evidence_type == "el":
        datatype = XSD.string
    g.add((node, RDF.type, URIRef(Cons.NODE_TYPES[f"{evidence_type}_node"])))
    g.add((node, URIRef(Cons.PREDICATES["sio_has_value"]), Literal(value, datatype=datatype)))
    return node
