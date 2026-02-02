#!/usr/bin/env python3

"""Populate a BDF RDF graph with AOP data."""

from typing import Optional

from rdflib import Graph, URIRef

from pyBiodatafuse.constants import NAMESPACE_BINDINGS, NODE_TYPES, PREDICATES, SOURCE_NAMESPACES
from pyBiodatafuse.graph.rdf.nodes.base import (
    add_label,
    add_triple,
    add_type,
    create_node,
    link_associated_with,
    safe_get,
)


def add_aop_data(
    g: Graph,
    entry: dict,
    compound_node: Optional[URIRef] = None,
    gene_node: Optional[URIRef] = None,
) -> None:
    """Add Adverse Outcome Pathway (AOP) data to the RDF graph.

    This function integrates AOP data into the provided RDF graph following
    the AOP ontology schema:
    - AOP has_key_event KE
    - AOP has_key_event_relationship KER
    - AOP has_adverse_outcome AO
    - AOP has_molecular_initiating_event MIE
    - KER has_upstream_key_event KE_upstream
    - KER has_downstream_key_event KE_downstream
    - KE OrganContext organ

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

    # Handle simple mode (ke, ke_title)
    ke_id = safe_get(entry, "ke")
    if ke_id and aop_node:
        ke_node = create_node(SOURCE_NAMESPACES["aopwiki"] + "aop.events/", ke_id)
        add_type(g, ke_node, NODE_TYPES["ke"])
        add_label(g, ke_node, safe_get(entry, "ke_title", ""))
        # AOP has_key_event KE
        add_triple(g, aop_node, URIRef(PREDICATES["has_key_event"]), ke_node)
        # Add organ context for simple KE if available
        ke_organ = safe_get(entry, "ke_organ")
        if ke_organ:
            ke_organ_node = URIRef(NAMESPACE_BINDINGS["obo"] + ke_organ.replace(":", "_"))
            add_triple(g, ke_node, URIRef(PREDICATES["organ_context"]), ke_organ_node)

    # Handle pathway mode - Add Adverse Outcome (AO)
    ao_id = safe_get(entry, "ao")
    if ao_id and aop_node:
        ao_node = create_node(SOURCE_NAMESPACES["aopwiki"] + "aop.events/", ao_id)
        add_type(g, ao_node, NODE_TYPES["ao"])
        add_label(g, ao_node, safe_get(entry, "ao_title", ""))
        # AOP has_adverse_outcome AO
        add_triple(g, aop_node, URIRef(PREDICATES["has_adverse_outcome"]), ao_node)

    # Handle pathway mode - Add Molecular Initiating Event (MIE)
    mie_id = safe_get(entry, "MIE")
    if mie_id and aop_node:
        mie_node = create_node(SOURCE_NAMESPACES["aopwiki"] + "aop.events/", mie_id)
        add_type(g, mie_node, NODE_TYPES["mie"])
        add_label(g, mie_node, safe_get(entry, "MIE_title", ""))
        # AOP has_molecular_initiating_event MIE
        add_triple(g, aop_node, URIRef(PREDICATES["has_molecular_initiating_event"]), mie_node)

    # Handle pathway mode - Add Key Event Relationship (KER)
    ker_id = safe_get(entry, "KER")
    ke_upstream_id = safe_get(entry, "KE_upstream")
    ke_downstream_id = safe_get(entry, "KE_downstream")

    if ker_id and aop_node:
        ker_node = create_node(SOURCE_NAMESPACES["aopwiki"] + "aop.relationships/", ker_id)
        add_type(g, ker_node, NODE_TYPES["ker"])

        # AOP has_key_event_relationship KER
        add_triple(g, aop_node, URIRef(PREDICATES["has_key_event_relationship"]), ker_node)

        # Create upstream key event
        if ke_upstream_id:
            ke_upstream_node = create_node(
                SOURCE_NAMESPACES["aopwiki"] + "aop.events/", ke_upstream_id
            )
            add_type(g, ke_upstream_node, NODE_TYPES["ke"])
            add_label(g, ke_upstream_node, safe_get(entry, "KE_upstream_title", ""))

            # KER has_upstream_key_event KE_upstream
            add_triple(g, ker_node, URIRef(PREDICATES["has_upstream_key_event"]), ke_upstream_node)

            # Add organ context if available
            ke_upstream_organ = safe_get(entry, "KE_upstream_organ")
            if ke_upstream_organ:
                ke_upstream_organ_node = URIRef(
                    NAMESPACE_BINDINGS["obo"] + ke_upstream_organ.replace(":", "_")
                )
                add_triple(
                    g, ke_upstream_node, URIRef(PREDICATES["organ_context"]), ke_upstream_organ_node
                )

        # Create downstream key event
        if ke_downstream_id:
            ke_downstream_node = create_node(
                SOURCE_NAMESPACES["aopwiki"] + "aop.events/", ke_downstream_id
            )
            add_type(g, ke_downstream_node, NODE_TYPES["ke"])
            add_label(g, ke_downstream_node, safe_get(entry, "KE_downstream_title", ""))

            # KER has_downstream_key_event KE_downstream
            add_triple(
                g, ker_node, URIRef(PREDICATES["has_downstream_key_event"]), ke_downstream_node
            )

            # Add organ context if available
            ke_downstream_organ = safe_get(entry, "KE_downstream_organ")
            if ke_downstream_organ:
                ke_downstream_organ_node = URIRef(
                    NAMESPACE_BINDINGS["obo"] + ke_downstream_organ.replace(":", "_")
                )
                add_triple(
                    g,
                    ke_downstream_node,
                    URIRef(PREDICATES["organ_context"]),
                    ke_downstream_organ_node,
                )

    # Link gene and compound nodes to AOP
    if gene_node and aop_node:
        link_associated_with(g, gene_node, aop_node, bidirectional=False)
    if compound_node and aop_node:
        link_associated_with(g, compound_node, aop_node, bidirectional=False)
