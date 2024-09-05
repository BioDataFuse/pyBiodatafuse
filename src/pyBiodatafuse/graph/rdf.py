# coding: utf-8

"""Python module to produce an RDF graph from the property table."""

import pandas as pd
import logging
from rdflib import Graph, URIRef, Literal, Namespace
from rdflib.namespace import RDF, RDFS, XSD, SKOS
from pyBiodatafuse.constants import (
    NAMESPACE_BINDINGS,
    URIS,
    NODE_TYPES,
    PREDICATES,
    BGEE_GENE_EXPRESSION_LEVELS_COL,
    DISGENET_DISEASE_COL
)
from bioregistry import get_iri, normalize_curie

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def create_gene_node(g: Graph, gene_base_node: str, id_number: str) -> URIRef:
    """Create and add a gene node to the RDF graph.

    :param g: RDF graph to which the gene node will be added.
    :param gene_base_node: Base URI for the gene node.
    :param id_number: Unique identifier for the gene node.
    :return: URIRef for the created gene node.
    """
    gene_node = URIRef(f"{gene_base_node}/{id_number}/")
    g.add((gene_node, RDF.type, URIRef(NODE_TYPES["gene_node"])))
    return gene_node


def create_disease_node(g: Graph, disease_data: dict) -> URIRef:
    """Create and add a disease node to the RDF graph.

    :param g: RDF graph to which the disease node will be added.
    :param disease_data: Dictionary containing disease information.
    :return: URIRef for the created disease node.
    """
    # UMLS IRIs not in Bioregistry
    disease_curie = disease_data.get("disease_umlscui")
    
    if disease_curie is None:
        logger.warning(f"Skipping disease node creation due to invalid identifier: {disease_data.get('disease_umlscui')}")
        return None
    else:
        disease_iri = f'https://www.ncbi.nlm.nih.gov/medgen/{disease_curie}'
        print(disease_iri)
        disease_node = URIRef(disease_iri)
        g.add((disease_node, RDF.type, URIRef(NODE_TYPES["disease_node"])))
        g.add((disease_node, RDFS.label, Literal(disease_data.get("disease_name", ""), datatype=XSD.string)))
        g.add((disease_node, SKOS.exactMatch, disease_node))
        return disease_node


def create_gene_disease_association_node(
    g: Graph, id_number: str, source_idx: str, disease_id: str
) -> URIRef:
    """Create and add a gene-disease association node to the RDF graph.

    :param g: RDF graph to which the association node will be added.
    :param id_number: Unique identifier for the association node.
    :param source_idx: Source index for the association.
    :param disease_id: Disease identifier associated with the gene.
    :return: URIRef for the created association node.
    """
    gene_disease_assoc_node = URIRef(
        f"{URIS['gene_disease_association']}/{id_number}/{source_idx}_assoc_{disease_id}"
    )
    g.add((gene_disease_assoc_node, RDF.type, URIRef(NODE_TYPES["gene-disease_association"])))
    return gene_disease_assoc_node


def create_score_node(
    g: Graph, id_number: str, source_idx: str, disease_id: str, score: float
) -> URIRef:
    """Create and add a score node to the RDF graph.

    :param g: RDF graph to which the score node will be added.
    :param id_number: Unique identifier for the score node.
    :param source_idx: Source index for the association.
    :param disease_id: Disease identifier associated with the score.
    :param score: Score value for the gene-disease association.
    :return: URIRef for the created score node.
    """
    score_node = URIRef(f"{URIS['score_base_node']}/{id_number}/{source_idx}_{disease_id}")
    g.add((score_node, RDF.type, URIRef(NODE_TYPES["score_node"])))
    g.add((score_node, URIRef(NAMESPACE_BINDINGS["sio"] + "has_value"), Literal(score, datatype=XSD.double)))
    return score_node


def create_evidence_source_node(
    g: Graph, id_number: str, source_idx: str, disease_id: str, evidence_source: str
) -> URIRef:
    """Create and add an evidence source node to the RDF graph.

    :param g: RDF graph to which the evidence source node will be added.
    :param id_number: Unique identifier for the evidence source node.
    :param source_idx: Source index for the association.
    :param disease_id: Disease identifier associated with the evidence source.
    :param evidence_source: Evidence source description.
    :return: URIRef for the created evidence source node.
    """
    DCAT = Namespace(NAMESPACE_BINDINGS['dcat'])
    evidence_source_node = URIRef(f"{URIS['source_base_node']}/{id_number}/{source_idx}_{disease_id}")
    g.add((evidence_source_node, RDF.type, DCAT.Dataset))
    g.add((evidence_source_node, RDFS.label, Literal(evidence_source, datatype=XSD.string)))
    return evidence_source_node


def create_data_source_node(g: Graph) -> URIRef:
    """Create and add a data source node to the RDF graph.

    :param g: RDF graph to which the data source node will be added.
    :return: URIRef for the created data source node.
    """
    data_source_name = Literal("DisGeNET", datatype=XSD.string)
    data_source_url = URIRef("https://disgenet.com/")
    g.add((data_source_url, RDF.type, URIRef(NODE_TYPES["data_source_node"])))
    g.add((data_source_url, RDFS.label, data_source_name))
    return data_source_url


def add_gene_disease_associations(
    g: Graph, id_number: str, source_idx: str, gene_node: URIRef, disease_data: list
) -> None:
    """Process and add gene-disease association data to the RDF graph.

    :param g: RDF graph to which the associations will be added.
    :param id_number: Unique identifier for the associations.
    :param source_idx: Source index for the associations.
    :param gene_node: URIRef of the gene node associated with the disease data.
    :param disease_data: List of dictionaries containing disease association information.
    """
    for data in disease_data:
        if (data['disease_umlscui']) != str:
            continue

        disease_node = create_disease_node(g, data)
        gene_disease_assoc_node = create_gene_disease_association_node(g, id_number, source_idx, data["disease_id"])
        g.add((gene_disease_assoc_node, URIRef(PREDICATES['sio_refers_to']), gene_node))
        g.add((gene_disease_assoc_node, URIRef(PREDICATES['sio_refers_to']), disease_node))

        score_node = create_score_node(g, id_number, source_idx, data["disease_id"], data["score"])
        g.add((gene_disease_assoc_node, URIRef(PREDICATES['sio_has_measurement_value']), score_node))

        evidence_source_node = create_evidence_source_node(g, id_number, source_idx, data["disease_id"], data["evidence_source"])
        g.add((gene_disease_assoc_node, URIRef(PREDICATES['sio_has_source']), evidence_source_node))

        data_source_node = create_data_source_node(g)
        g.add((gene_disease_assoc_node, URIRef(PREDICATES['sio_has_source']), data_source_node))


def add_gene_expression_data(
    g: Graph, id_number: str, source_idx: str, gene_node: URIRef, expression_data: list
) -> None:
    """Process and add gene expression data to the RDF graph.

    :param g: RDF graph to which the expression data will be added.
    :param id_number: Unique identifier for the expression data.
    :param source_idx: Source index for the expression data.
    :param gene_node: URIRef of the gene node associated with the expression data.
    :param expression_data: List of dictionaries containing gene expression information.
    """
    for data in expression_data:
        if pd.isna(data["anatomical_entity_id"]):
            continue

        exp_process_node = URIRef(f"{URIS['experimental_process_node']}/{id_number}/{source_idx}")
        developmental_stage_node = URIRef(get_iri(data["developmental_stage_id"].replace("_", ":")))
        anatomical_entity_node = URIRef(get_iri(data["anatomical_entity_id"].replace("_", ":")))
        gene_expression_value_node = URIRef(
            f"{URIS['gene_expression_value_base_node']}/{id_number}/{source_idx}_{data['anatomical_entity_id']}"
        )

        g.add((gene_node, URIRef(PREDICATES['sio_is_associated_with']), gene_expression_value_node))
        g.add((gene_node, URIRef(PREDICATES['sio_is_associated_with']), anatomical_entity_node))
        g.add((gene_node, URIRef(PREDICATES['sio_is_associated_with']), developmental_stage_node))

        g.add((gene_expression_value_node, RDF.type, URIRef(NODE_TYPES['gene_expression_value_node'])))
        g.add((gene_expression_value_node, URIRef(PREDICATES['sio_has_value']),
              Literal(data["expression_level"], datatype=XSD.double)))

        g.add((anatomical_entity_node, RDF.type, URIRef(NODE_TYPES["anatomical_entity_node"])))
        g.add((anatomical_entity_node, RDFS.label, Literal(data["anatomical_entity_name"], datatype=XSD.string)))
        g.add((developmental_stage_node, RDFS.label, Literal(data["developmental_stage_name"], datatype=XSD.string)))
        g.add((anatomical_entity_node, SKOS.exactMatch, anatomical_entity_node))
        g.add((URIRef(URIS["life_cycle_base_node"]), SKOS.exactMatch, developmental_stage_node))

        g.add((exp_process_node, RDF.type, URIRef(NODE_TYPES['gene_expression_value_node'])))
        g.add((exp_process_node, URIRef(PREDICATES['sio_has_input']), gene_node))
        g.add((exp_process_node, URIRef(PREDICATES['sio_has_output']), gene_expression_value_node))
        g.add((exp_process_node, URIRef(PREDICATES['sio_has_input']), anatomical_entity_node))


def generate_rdf(df: pd.DataFrame) -> Graph:
    """Generate RDF graph from the property table.

    :param df: DataFrame containing the data.
    :return: RDF graph.
    """
    g = Graph()
    # Get namespace bindings
    for key, value in NAMESPACE_BINDINGS.items():
        logger.info(f'Binding {key} to {value}')
        g.bind(key, value)
    for key, value in URIS.items():
        g.bind(key, value)

    for i, row in df.iterrows():
        # Unpack the relevant columns
        source_idx = row['identifier']
        source_namespace = row['identifier.source']
        target_idx = row['target']
        target_namespace = row['target.source']
        expression_data = row[BGEE_GENE_EXPRESSION_LEVELS_COL]
        disease_data = row[DISGENET_DISEASE_COL]

        # Check if any of the essential columns are NaN, and skip the iteration if so
        if pd.isna(source_idx) or pd.isna(source_namespace) or pd.isna(target_idx) or pd.isna(target_namespace):
            logger.warning(f"Skipping row {i} due to missing essential data.")
            continue

        # Normalize CURIEs and convert to IRIs
        source_curie = normalize_curie(f"{source_namespace}:{source_idx}")
        target_curie = normalize_curie(f"{target_namespace}:{target_idx}")

        # Skip the row if CURIE normalization fails (e.g., if it's None)
        if not source_curie or not target_curie:
            logger.warning(f"Skipping row {i} due to CURIE normalization failure.")
            continue

        id_number = f"{i:06d}"  # Create ID for this row
        gene_node = create_gene_node(g, URIS["gene_base_node"], id_number)

        # Add gene-disease associations
        add_gene_disease_associations(g, id_number, source_idx, gene_node, disease_data)

        # Add gene expression data
        add_gene_expression_data(g, id_number, source_idx, gene_node, expression_data)

    return g
