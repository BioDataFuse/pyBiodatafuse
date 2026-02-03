# coding: utf-8

"""
Base utilities for RDF node generation.

This module provides common utilities and base classes for creating RDF nodes
in the BioDataFuse knowledge graph. It establishes consistent patterns for
node creation, triple generation, and data extraction.

Usage
-----
Import the utilities you need::

    from pyBiodatafuse.graph.rdf.nodes.base import (
        create_node,
        add_label,
        add_type,
        add_triple,
        get_uri,
    )

Extending
---------
To add a new node type:

1. Create a new module in the nodes/ directory
2. Import utilities from this base module
3. Define functions that use create_node() and add_triple()
4. Register any new node types in constants.py

Example::

    from pyBiodatafuse.graph.rdf.nodes.base import create_node, add_label, add_type
    import pyBiodatafuse.constants as Cons

    def add_my_node(g, data, parent_node):
        node_id = data.get("id")
        if not node_id:
            return None

        node = create_node(Cons.NODE_URI_PREFIXES["my_source"], node_id)
        add_type(g, node, Cons.NODE_TYPES["my_node_type"])
        add_label(g, node, data.get("name", node_id))
        add_triple(g, parent_node, Cons.PREDICATES["has_part"], node)

        return node
"""

from typing import Any, Dict, List, Optional, Tuple, Union

from rdflib import Graph, Literal, URIRef
from rdflib.namespace import OWL, RDF, RDFS, XSD

import pyBiodatafuse.constants as Cons

# =============================================================================
# Node Creation
# =============================================================================


def create_node(prefix: str, identifier: str) -> URIRef:
    """
    Create a URIRef node from a prefix and identifier.

    :param prefix: URI prefix for the node.
    :param identifier: Local identifier for the node.
    :return: URIRef for the node.
    """
    return URIRef(f"{prefix}{identifier}")


def get_uri(source: str, identifier: str) -> Optional[URIRef]:
    """
    Get a URIRef for a source and identifier using configured prefixes.

    :param source: Source name (must be in NODE_URI_PREFIXES or SOURCE_NAMESPACES).
    :param identifier: Local identifier.
    :return: URIRef or None if source not found.
    """
    # Normalize sources (patch for bdb)
    if source == "PubChem-compound":
        source = "PubChem Compound"
    prefix = Cons.NODE_URI_PREFIXES.get(source) or Cons.SOURCE_NAMESPACES.get(source)
    if prefix and identifier:
        return URIRef(f"{prefix}{identifier}")
    return None


# =============================================================================
# Triple Addition
# =============================================================================


def add_triple(
    g: Graph,
    subject: Optional[URIRef],
    predicate: Union[str, URIRef],
    obj: Union[URIRef, Literal, str],
) -> bool:
    """
    Add a triple to the graph with predicate lookup.

    :param g: RDF graph.
    :param subject: Subject URIRef (can be None, will be ignored).
    :param predicate: Predicate URI string (looked up in PREDICATES) or URIRef.
    :param obj: Object (URIRef, Literal, or string to convert to Literal).
    :return: True if triple was added.
    """
    if subject is None:
        return False

    # Resolve predicate
    if isinstance(predicate, str):
        pred_uri = Cons.PREDICATES.get(predicate)
        if pred_uri:
            predicate = URIRef(pred_uri)
        else:
            predicate = URIRef(predicate)

    # Convert string objects to Literals (but preserve URIRef objects)
    if isinstance(obj, str) and not isinstance(obj, URIRef):
        obj = Literal(obj, datatype=XSD.string)

    g.add((subject, predicate, obj))
    return True


def add_type(g: Graph, node: Optional[URIRef], node_type: str) -> None:
    """
    Add an rdf:type triple for a node.

    :param g: RDF graph.
    :param node: Node URIRef (can be None, will be ignored).
    :param node_type: Type URI string.
    """
    if node and node_type:
        g.add((node, RDF.type, URIRef(node_type)))


def add_label(g: Graph, node: Optional[URIRef], label: Optional[str]) -> None:
    """
    Add an rdfs:label triple for a node.

    :param g: RDF graph.
    :param node: Node URIRef (can be None, will be ignored).
    :param label: Label string.
    """
    if node and label:
        g.add((node, RDFS.label, Literal(str(label), datatype=XSD.string)))


def add_same_as(
    g: Graph, node1: Optional[URIRef], node2: Optional[URIRef], bidirectional: bool = True
) -> None:
    """
    Add owl:sameAs relationship between two nodes.

    :param g: RDF graph.
    :param node1: First node (can be None, will be ignored).
    :param node2: Second node (can be None, will be ignored).
    :param bidirectional: If True, add both directions.
    """
    if node1 and node2:
        g.add((node1, OWL.sameAs, node2))
        if bidirectional:
            g.add((node2, OWL.sameAs, node1))


def add_value(g: Graph, node: Optional[URIRef], value: Any, datatype: URIRef = XSD.string) -> None:
    """
    Add a sio:has-value triple for a node.

    :param g: RDF graph.
    :param node: Node URIRef (can be None, will be ignored).
    :param value: Value to add.
    :param datatype: XSD datatype for the literal.
    """
    if node and value is not None:
        g.add((node, URIRef(Cons.PREDICATES["sio_has_value"]), Literal(value, datatype=datatype)))


# =============================================================================
# Relationship Helpers
# =============================================================================


def link_has_part(
    g: Graph, parent: Optional[URIRef], child: Optional[URIRef], bidirectional: bool = True
) -> None:
    """
    Add has-part/is-part-of relationship between nodes.

    :param g: RDF graph.
    :param parent: Parent node (can be None, will be ignored).
    :param child: Child node (can be None, will be ignored).
    :param bidirectional: If True, add both directions.
    """
    if parent and child:
        g.add((parent, URIRef(Cons.PREDICATES["sio_has_part"]), child))
        if bidirectional:
            g.add((child, URIRef(Cons.PREDICATES["sio_is_part_of"]), parent))


def link_refers_to(g: Graph, subject: Optional[URIRef], target: Optional[URIRef]) -> None:
    """
    Add sio:refers-to relationship.

    :param g: RDF graph.
    :param subject: Subject node (can be None, will be ignored).
    :param target: Target node (can be None, will be ignored).
    """
    if subject and target:
        g.add((subject, URIRef(Cons.PREDICATES["sio_refers_to"]), target))


def link_has_source(g: Graph, node: Optional[URIRef], source: Optional[URIRef]) -> None:
    """
    Add sio:has-source relationship.

    :param g: RDF graph.
    :param node: Node to link (can be None, will be ignored).
    :param source: Source node (can be None, will be ignored).
    """
    if node and source:
        g.add((node, URIRef(Cons.PREDICATES["sio_has_source"]), source))


def link_associated_with(
    g: Graph, node1: Optional[URIRef], node2: Optional[URIRef], bidirectional: bool = True
) -> None:
    """
    Add sio:is-associated-with relationship.

    :param g: RDF graph.
    :param node1: First node (can be None, will be ignored).
    :param node2: Second node (can be None, will be ignored).
    :param bidirectional: If True, add both directions.
    """
    if node1 and node2:
        g.add((node1, URIRef(Cons.PREDICATES["sio_is_associated_with"]), node2))
        if bidirectional:
            g.add((node2, URIRef(Cons.PREDICATES["sio_is_associated_with"]), node1))


# =============================================================================
# Data Extraction Helpers
# =============================================================================


def safe_get(data: Dict[str, Any], key: str, default: Any = None) -> Any:
    """
    Safely get a value from a dictionary.

    :param data: Dictionary to get value from.
    :param key: Key to look up.
    :param default: Default value if key not found or value is None/NaN.
    :return: Value or default.
    """
    import pandas as pd

    value = data.get(key, default)
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return default
    return value


def split_curie(curie: str) -> Tuple[Optional[str], Optional[str]]:
    """
    Split a CURIE into prefix and local part.

    :param curie: CURIE string (e.g., "CHEBI:12345").
    :return: Tuple of (prefix, local) or (None, None) if invalid.
    """
    if not curie or ":" not in curie:
        return (None, None)
    parts = curie.split(":", 1)
    return (parts[0], parts[1]) if len(parts) == 2 else (None, None)


def extract_id(value: str, prefix: Optional[str] = None) -> str:
    """
    Extract an ID from a value, optionally removing a prefix.

    :param value: Value string.
    :param prefix: Prefix to remove (e.g., "CHEBI:").
    :return: Extracted ID.
    """
    if not value:
        return ""
    if prefix and value.startswith(prefix):
        return value[len(prefix) :]
    if ":" in value:
        return value.split(":", 1)[1]
    return value


# =============================================================================
# Batch Operations
# =============================================================================


def add_triples(g: Graph, triples: List[Tuple[URIRef, URIRef, Any]]) -> int:
    """
    Add multiple triples to the graph.

    :param g: RDF graph.
    :param triples: List of (subject, predicate, object) tuples.
    :return: Number of triples added.
    """
    count = 0
    for triple in triples:
        if all(t is not None for t in triple):
            g.add(triple)
            count += 1
    return count
