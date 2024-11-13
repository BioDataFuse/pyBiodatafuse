# gene_expression.py

"""Populate a BDF RDF graph with gene expression data."""


from bioregistry import get_iri
from rdflib import Graph, Literal, URIRef
from rdflib.namespace import RDF, RDFS, XSD

from pyBiodatafuse.constants import NODE_TYPES, PREDICATES
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
        if not data.get("anatomical_entity_id", None):
            continue

        gene_expression_value_node = URIRef(
            f"{new_uris['gene_expression_value_base_node']}/{id_number}/{source_idx}"
        )
        developmental_stage_node = URIRef(get_iri(data["developmental_stage_id"].replace("_", ":")))
        anatomical_entity_node = URIRef(get_iri(data["anatomical_entity_id"].replace("_", ":")))
        exp_uri = new_uris["gene_expression_value_base_node"]
        anatomical_entity = data["anatomical_entity_id"]
        gene_expression_value_node = URIRef(
            f"{exp_uri}/{id_number}/{source_idx}_{anatomical_entity}"
        )

        g.add(
            (gene_node, URIRef(PREDICATES["sio_has_measurement_value"]), gene_expression_value_node)
        )
        g.add(
            (
                gene_expression_value_node,
                URIRef(PREDICATES["sio_is_associated_with"]),
                anatomical_entity_node,
            )
        )
        g.add(
            (
                gene_expression_value_node,
                URIRef(PREDICATES["sio_is_associated_with"]),
                developmental_stage_node,
            )
        )

        # g.add((gene_node, URIRef(PREDICATES['sio_is_part_of']), cellular_component))
        # g.add((gene_node, URIRef(PREDICATES['sio_is_represented_by']), gene_symbol))

        g.add(
            (gene_expression_value_node, RDF.type, URIRef(NODE_TYPES["gene_expression_value_node"]))
        )
        g.add(
            (
                gene_expression_value_node,
                URIRef(PREDICATES["sio_has_value"]),
                Literal(data["expression_level"], datatype=XSD.double),
            )
        )

        g.add((anatomical_entity_node, RDF.type, URIRef(NODE_TYPES["anatomical_entity_node"])))
        g.add(
            (
                anatomical_entity_node,
                RDFS.label,
                Literal(data["anatomical_entity_name"], datatype=XSD.string),
            )
        )
        g.add(
            (
                developmental_stage_node,
                RDFS.label,
                Literal(data["developmental_stage_name"], datatype=XSD.string),
            )
        )
        g.add(
            (
                developmental_stage_node,
                RDF.type,
                URIRef(NODE_TYPES["developmental_stage_node"]),
            )
        )
        g.add(
            (gene_expression_value_node, RDF.type, URIRef(NODE_TYPES["gene_expression_value_node"]))
        )
        g.add((gene_expression_value_node, URIRef(PREDICATES["sio_has_input"]), gene_node))

        # Add experimental node
        if experimental_process_data:
            for exp in experimental_process_data:
                experimental_process_node = add_experimental_process_node(
                    g=g,
                    data=exp,
                )
                if experimental_process_node:
                    g.add(
                        (gene_node, URIRef(PREDICATES["sio_is_part_of"]), experimental_process_node)
                    )
                    g.add(
                        (experimental_process_node, URIRef(PREDICATES["sio_has_part"]), gene_node)
                    )
                    g.add(
                        (
                            experimental_process_node,
                            URIRef(PREDICATES["sio_has_part"]),
                            anatomical_entity_node,
                        )
                    )
                    data_source_node = add_data_source_node(g, "Bgee")
                    g.add(
                        (
                            gene_expression_value_node,
                            URIRef(PREDICATES["sio_has_source"]),
                            data_source_node,
                        )
                    )
