# gene_disease.py

"""Module to populate a BDF RDF graph with gene-disease relationship data."""


from bioregistry import get_iri
from rdflib import Graph, Literal, URIRef
from rdflib.namespace import OWL, RDF, RDFS, XSD

from pyBiodatafuse.constants import DISEASE_IDENTIFIER_TYPES, NODE_TYPES, PREDICATES
from pyBiodatafuse.graph.rdf.metadata import (  # pylint: disable=no-name-in-module
    add_data_source_node,
)


def add_gene_disease_associations(
    g: Graph,
    id_number: str,
    source_idx: str,
    gene_node: URIRef,
    disease_data: dict,
    new_uris: dict,
    i: int,
) -> None:
    """Process and add gene-disease association data to the RDF graph.

    :param g: RDF graph to which the associations will be added.
    :param id_number: Unique identifier for the associations.
    :param source_idx: Source index for the associations.
    :param gene_node: URIRef of the gene node associated with the disease data.
    :param disease_data: List of dictionaries containing disease association information.
    :param new_uris: Dictionary with updated project base URIs for the nodes.
    :param i: the index of the row
    """
    data = disease_data
    gene_id = g.value(subject=gene_node, predicate=RDFS.label)
    disease_node = add_disease_node(g, data)
    if disease_node:
        gene_disease_assoc_node = URIRef(
            f"{new_uris['gene_disease_association']}/{id_number}{i}{source_idx}"
        )
        g.add((gene_disease_assoc_node, RDF.type, URIRef(NODE_TYPES["gene_disease_association"])))
        g.add((gene_disease_assoc_node, URIRef(PREDICATES["sio_refers_to"]), gene_node))
        g.add((gene_disease_assoc_node, URIRef(PREDICATES["sio_refers_to"]), disease_node))
        if data.get("score"):
            disease_umlscui = disease_data.get("disease_umlscui", None)
            if disease_umlscui:
                score_node = add_score_node(
                    g=g,
                    id_number=id_number,
                    source_idx=source_idx,
                    score=data["score"],
                    new_uris=new_uris,
                    i=i,
                    disease_id=disease_umlscui,
                    gene_id=gene_id,
                )
                g.add(
                    (
                        gene_disease_assoc_node,
                        URIRef(PREDICATES["sio_has_measurement_value"]),
                        score_node,
                    )
                )
        data_source_node = add_data_source_node(g, "DISGENET")
        g.add((gene_disease_assoc_node, URIRef(PREDICATES["sio_has_source"]), data_source_node))


def add_evidence_idx_node(
    g: Graph, id_number: str, source_idx: str, disease_id: str, ei: float, new_uris: dict, i
) -> URIRef:
    """Create and add an evidence index (EI) node to the RDF graph.

    :param g: RDF graph to which the EI node will be added.
    :param id_number: Unique identifier for the EI node.
    :param source_idx: Source index for the association.
    :param disease_id: Disease identifier associated with the EI.
    :param ei: Evidence index value for the gene-disease association.
    :param new_uris: Dictionary with updated project base URIs for the nodes.
    :param i: Index or iteration number used for node URI construction.
    :return: URIRef for the created EI node.
    """
    evidence_idx_node = URIRef(
        f"{new_uris['score_base_node']}/{id_number}{i}{source_idx}_{disease_id}"
    )
    g.add((evidence_idx_node, RDF.type, URIRef(NODE_TYPES["evidence_idx_node"])))
    g.add(
        (
            evidence_idx_node,
            URIRef(PREDICATES["sio_has_value"]),
            Literal(ei, datatype=XSD.double),
        )
    )
    return evidence_idx_node


def add_disease_node(g: Graph, disease_data: dict) -> URIRef:
    """Create and add a disease node to the RDF graph.

    :param g: RDF graph to which the disease node will be added.
    :param disease_data: Dictionary containing disease information.
    :return: URIRef for the created disease node.
    """
    # UMLS IRIs not in Bioregistry
    disease_curie = disease_data.get("disease_umlscui")
    if disease_curie is None:
        disease_curie = disease_data.get("UMLS")
    disease_iri = f"https://www.ncbi.nlm.nih.gov/medgen/{disease_curie}"
    disease_data.get()
    disease_node = URIRef(disease_iri)
    g.add((disease_node, RDF.type, URIRef(NODE_TYPES["disease_node"])))
    g.add(
        (disease_node, RDFS.label, Literal(disease_data.get("disease_name"), datatype=XSD.string))
    )
    for identifier_type in DISEASE_IDENTIFIER_TYPES:
        curie_field = disease_data.get(identifier_type, None)
        if curie_field and "ncit" not in curie_field:
            curies = [curie_field]
            if "," in (curie_field):
                curies = [i for i in curie_field.split(", ")]
            for curie in curies:
                disease_source_iri = get_iri(curie)
                if disease_source_iri is None:
                    if ":" in curie:
                        curie = curie.split(":")[1]
                    disease_source_iri = get_iri("obo:" + curie)
                    g.add(
                        (disease_node, OWL.SameAs, URIRef(disease_source_iri))
                    )  # Some of the data does not look like a skos:exactMatch
                else:
                    g.add((disease_node, OWL.sameAs, URIRef(disease_source_iri)))
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
    :param gene_id: String value of the gene ID
    :return: URIRef for the created score node.
    """
    score_node = URIRef(
        f"{new_uris['score_base_node']}/{id_number}{i}{source_idx}_{disease_id}{gene_id}"
    )
    g.add((score_node, RDF.type, URIRef(NODE_TYPES["score_node"])))
    g.add(
        (
            score_node,
            URIRef(PREDICATES["sio_has_value"]),
            Literal(score, datatype=XSD.double),
        )
    )
    return score_node
