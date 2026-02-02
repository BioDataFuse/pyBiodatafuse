# gene_expression.py

"""Populate a BDF RDF graph with gene expression data."""

from bioregistry import get_iri
from rdflib import Graph, URIRef
from rdflib.namespace import XSD

import pyBiodatafuse.constants as Cons
from pyBiodatafuse.graph.rdf.nodes.base import (
    add_label,
    add_triple,
    add_type,
    add_value,
    create_node,
    link_has_part,
    link_has_source,
    safe_get,
)
from pyBiodatafuse.graph.rdf.nodes.experimental_process import add_experimental_process_node
from pyBiodatafuse.graph.rdf.utils import add_data_source_node


def add_gene_expression_data(
    g: Graph,
    id_number: str,
    source_idx: str,
    gene_node: URIRef,
    expression_data: list,
    experimental_process_data: list,
    new_uris: dict,
) -> None:
    """Process and add gene expression data to the RDF graph.

    :param g: RDF graph to which the expression data will be added.
    :param id_number: Unique identifier for the expression data.
    :param source_idx: Source index for the expression data.
    :param gene_node: URIRef of the gene node associated with the expression data.
    :param expression_data: List of dictionaries containing gene expression information.
    :param experimental_process_data: List of experimental process data.
    :param new_uris: Dictionary with updated project base URIs for the nodes.
    """
    for data in expression_data:
        anatomical_id = safe_get(data, Cons.ANATOMICAL_ID)
        if not anatomical_id:
            continue

        # Get developmental stage node
        dev_id = safe_get(data, Cons.DEVELOPMENTAL_ID)
        if not dev_id:
            continue
        developmental_stage_iri = get_iri(dev_id.replace("_", ":"))
        if not developmental_stage_iri:
            continue
        developmental_stage_node = URIRef(developmental_stage_iri)

        # Get anatomical entity node
        anatomical_entity_iri = get_iri(anatomical_id.replace("_", ":"))
        if not anatomical_entity_iri:
            continue
        anatomical_entity_node = URIRef(anatomical_entity_iri)

        # Create expression value node
        exp_uri = new_uris["gene_expression_value_base_node"]
        gene_expression_value_node = create_node(
            f"{exp_uri}/{id_number}/", f"{source_idx}_{anatomical_id}"
        )

        # Link gene to expression value
        add_triple(g, gene_node, "sio_has_measurement_value", gene_expression_value_node)
        add_triple(g, gene_expression_value_node, "sio_is_associated_with", anatomical_entity_node)
        add_triple(
            g, gene_expression_value_node, "sio_is_associated_with", developmental_stage_node
        )

        # Add types and labels
        add_type(g, gene_expression_value_node, Cons.NODE_TYPES["gene_expression_value_node"])
        expression_level = safe_get(data, "expression_level")
        if expression_level is not None:
            add_value(g, gene_expression_value_node, expression_level, XSD.double)

        add_type(g, anatomical_entity_node, Cons.NODE_TYPES["anatomical_entity_node"])
        add_label(g, anatomical_entity_node, safe_get(data, Cons.ANATOMICAL_NAME))

        add_type(g, developmental_stage_node, Cons.NODE_TYPES["developmental_stage_node"])
        add_label(g, developmental_stage_node, safe_get(data, Cons.DEVELOPMENTAL_STAGE_NAME))

        # Link input
        add_triple(g, gene_expression_value_node, "sio_has_input", gene_node)

        # Add experimental process nodes
        if experimental_process_data:
            for exp in experimental_process_data:
                experimental_process_node = add_experimental_process_node(g=g, data=exp)
                if experimental_process_node:
                    link_has_part(g, experimental_process_node, gene_node)
                    link_has_part(g, experimental_process_node, anatomical_entity_node)
                    data_source_node = add_data_source_node(g, "Bgee")
                    link_has_source(g, gene_expression_value_node, data_source_node)
