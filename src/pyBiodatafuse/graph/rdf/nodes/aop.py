#!/usr/bin/env python3

"""Populate a BDF RDF g with AOP data."""


from rdflib import Graph, Literal, URIRef
from rdflib.namespace import OWL, RDF, RDFS, XSD

from pyBiodatafuse.constants import NAMESPACE_BINDINGS, NODE_TYPES, PREDICATES, SOURCE_NAMESPACES
from pyBiodatafuse.graph.rdf.nodes.compound import add_associated_compound_node


def add_aop_data(g, entry, compound_node=None, gene_node=None):
    """
    Add Adverse Outcome Pathway (AOP) data to the RDF graph.

    This function integrates AOP data into the provided RDF graph by linking
    compound and gene nodes with the specified entry data. It also handles the creation
    of new URIs based on the base URI and updates the g accordingly.

    :param g: The RDF graph to which the AOP data will be added.
    :param compound_node: The RDF node representing the compound, defaults to None.
    :param gene_node: The RDF node representing the gene, defaults to None.
    :param entry: The data entry containing information about the AOP.
    """
    # Create AOP node
    aop_node = (
        URIRef(SOURCE_NAMESPACES["aopwiki"] + "aop/" + entry["aop"]) if entry.get("aop") else None
    )
    if aop_node:
        g.add((aop_node, RDF.type, URIRef(NODE_TYPES["aop"])))
        g.add((aop_node, RDFS.label, Literal(entry.get("aop_title", ""))))
    # Add key events
    if entry.get("KE_upstream"):
        ke_upstream_node = URIRef(
            SOURCE_NAMESPACES["aopwiki"] + "aop.events/" + entry["KE_upstream"]
        )
        g.add((ke_upstream_node, RDF.type, URIRef(NODE_TYPES["ke"])))
        g.add((ke_upstream_node, RDFS.label, Literal(entry.get("KE_upstream_title", ""))))
        if entry.get("KE_upstream_organ"):
            ke_upstream_organ = URIRef(
                NAMESPACE_BINDINGS["obo"] + entry["KE_upstream_organ"].replace(":", "/")
            )
            g.add((ke_upstream_node, URIRef(PREDICATES["occurs_in"]), ke_upstream_organ))
        if aop_node:
            g.add((ke_upstream_node, URIRef(PREDICATES["has_downstream_key_event"]), aop_node))
    if entry.get("KE_downstream"):
        ke_downstream_node = URIRef(
            SOURCE_NAMESPACES["aopwiki"] + "aop.events/" + entry["KE_downstream"]
        )
        g.add((ke_downstream_node, RDF.type, URIRef(NODE_TYPES["ke"])))
        g.add((ke_downstream_node, RDFS.label, Literal(entry.get("KE_downstream_title", ""))))
        if entry.get("KE_downstream_organ"):
            ke_downstream_organ = URIRef(
                NAMESPACE_BINDINGS["obo"] + entry["KE_downstream_organ"].replace(":", "/")
            )
            g.add((ke_downstream_node, URIRef(PREDICATES["occurs_in"]), ke_downstream_organ))
        if aop_node:
            g.add((ke_downstream_node, URIRef(PREDICATES["has_upstream_key_event"]), aop_node))
    if gene_node:
        g.add((gene_node, URIRef(PREDICATES["sio_is_associated_with"]), aop_node))
    if compound_node:
        g.add((compound_node, URIRef(PREDICATES["sio_is_associated_with"]), aop_node))
