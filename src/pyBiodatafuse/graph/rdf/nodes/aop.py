#!/usr/bin/env python3

"""Populate a BDF RDF graph with AOP data."""

from typing import Optional

from rdflib import Graph, URIRef

from pyBiodatafuse.constants import NAMESPACE_BINDINGS, NODE_TYPES, SOURCE_NAMESPACES
from pyBiodatafuse.graph.rdf.nodes.base import (
    add_label,
    add_triple,
    add_type,
    create_node,
    link_associated_with,
    safe_get,
)


def add_aop_data(
    g: Graph, entry: dict, compound_node: Optional[URIRef] = None, gene_node: Optional[URIRef] = None
) -> None:
    """Add Adverse Outcome Pathway (AOP) data to the RDF graph.

    This function integrates AOP data into the provided RDF graph by linking
    compound and gene nodes with the specified entry data.

    :param g: The RDF graph to which the AOP data will be added.
    :param entry: The data entry containing information about the AOP.
    :param compound_node: The RDF node representing the compound, defaults to None.
    :param gene_node: The RDF node representing the gene, defaults to None.
    """
    aop_id = safe_get(entry, "aop")
    aop_node = None

    # Create AOP node
    if aop_id:
        aop_node = create_node(SOURCE_NAMESPACES["aopwiki"] + "aop/", aop_id)
        add_type(g, aop_node, NODE_TYPES["aop"])
        add_label(g, aop_node, safe_get(entry, "aop_title", ""))

    # Add upstream key event
    ke_upstream = safe_get(entry, "KE_upstream")
    if ke_upstream:
        ke_upstream_node = create_node(
            SOURCE_NAMESPACES["aopwiki"] + "aop.events/", ke_upstream
        )
        add_type(g, ke_upstream_node, NODE_TYPES["ke"])
        add_label(g, ke_upstream_node, safe_get(entry, "KE_upstream_title", ""))

        ke_upstream_organ = safe_get(entry, "KE_upstream_organ")
        if ke_upstream_organ:
            ke_upstream_organ_node = URIRef(
                NAMESPACE_BINDINGS["obo"] + ke_upstream_organ.replace(":", "/")
            )
            add_triple(g, ke_upstream_node, "occurs_in", ke_upstream_organ_node)

        if aop_node:
            add_triple(g, ke_upstream_node, "has_downstream_key_event", aop_node)

    # Add downstream key event
    ke_downstream = safe_get(entry, "KE_downstream")
    if ke_downstream:
        ke_downstream_node = create_node(
            SOURCE_NAMESPACES["aopwiki"] + "aop.events/", ke_downstream
        )
        add_type(g, ke_downstream_node, NODE_TYPES["ke"])
        add_label(g, ke_downstream_node, safe_get(entry, "KE_downstream_title", ""))

        ke_downstream_organ = safe_get(entry, "KE_downstream_organ")
        if ke_downstream_organ:
            ke_downstream_organ_node = URIRef(
                NAMESPACE_BINDINGS["obo"] + ke_downstream_organ.replace(":", "/")
            )
            add_triple(g, ke_downstream_node, "occurs_in", ke_downstream_organ_node)

        if aop_node:
            add_triple(g, ke_downstream_node, "has_upstream_key_event", aop_node)

    # Link gene and compound nodes to AOP
    if gene_node and aop_node:
        link_associated_with(g, gene_node, aop_node, bidirectional=False)
    if compound_node and aop_node:
        link_associated_with(g, compound_node, aop_node, bidirectional=False)
