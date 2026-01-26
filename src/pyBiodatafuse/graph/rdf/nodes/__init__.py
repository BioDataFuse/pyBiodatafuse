# coding: utf-8

"""
RDF Node Generation Modules for BioDataFuse.

This package contains modules for creating RDF nodes and relationships
in the BioDataFuse knowledge graph. Each module handles a specific
type of biological entity or relationship.

Module Structure
----------------
- base.py: Common utilities for node creation (create_node, add_label, etc.)
- gene.py: Gene nodes from Ensembl identifiers
- compound.py: Compound nodes from PubChem, ChEMBL, etc.
- protein_protein.py: Protein-protein interaction nodes from StringDB
- gene_disease.py: Gene-disease association nodes from DisGeNET/OpenTargets
- gene_expression.py: Gene expression nodes from Bgee
- pathway.py: Pathway nodes from WikiPathways, MINERVA, Reactome
- go_terms.py: Gene Ontology term nodes from OpenTargets
- literature.py: Literature-based association nodes from PubChem
- aop.py: Adverse Outcome Pathway nodes from AOP-Wiki
- dataset_provenance.py: Dataset metadata and provenance nodes

Usage Pattern
-------------
Each module follows a consistent pattern::

    from pyBiodatafuse.graph.rdf.nodes.base import create_node, add_label, add_type
    import pyBiodatafuse.constants as Cons

    def add_entity_node(g, data, parent_node=None):
        # 1. Extract and validate required data
        entity_id = data.get("id")
        if not entity_id:
            return None

        # 2. Create the node
        node = create_node(Cons.NODE_URI_PREFIXES["source"], entity_id)

        # 3. Add type and label
        add_type(g, node, Cons.NODE_TYPES["entity_type"])
        add_label(g, node, data.get("name", entity_id))

        # 4. Add relationships
        if parent_node:
            link_has_part(g, parent_node, node)

        # 5. Return the node
        return node

Adding New Node Types
---------------------
1. Create a new module in this package (e.g., new_entity.py)
2. Import utilities from base.py
3. Add any new node types to constants.NODE_TYPES
4. Add any new predicates to constants.PREDICATES
5. Implement the add_*_node function following the pattern above
6. Import and use in rdf.py process methods

Constants Required
------------------
When adding new node types, ensure these are defined in constants.py:

- NODE_URI_PREFIXES: URI prefix for the entity source
- NODE_TYPES: RDF type URI for the entity
- PREDICATES: Any new predicates needed
- DATA column name for extraction

Testing
-------
Each module should have corresponding tests in tests/graph/test_*.py
that verify:

1. Node creation with valid data
2. Graceful handling of missing/invalid data
3. Correct RDF types and labels
4. Proper relationship creation
"""

# Re-export commonly used utilities
from pyBiodatafuse.graph.rdf.nodes.base import (
    add_label,
    add_same_as,
    add_triple,
    add_triples,
    add_type,
    add_value,
    create_node,
    extract_id,
    get_uri,
    link_associated_with,
    link_has_part,
    link_has_source,
    link_refers_to,
    safe_get,
    split_curie,
)

__all__ = [
    # Base utilities
    "create_node",
    "get_uri",
    "add_triple",
    "add_triples",
    "add_type",
    "add_label",
    "add_same_as",
    "add_value",
    "link_associated_with",
    "link_has_part",
    "link_refers_to",
    "link_has_source",
    "safe_get",
    "split_curie",
    "extract_id",
]
