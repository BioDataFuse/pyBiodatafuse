# experimental_process.py

"""Populate a BDF RDF graph with PubChem assay/experimental process data."""

from typing import Optional

from rdflib import Graph, URIRef

import pyBiodatafuse.constants as Cons
from pyBiodatafuse.graph.rdf.nodes.base import (
    add_label,
    add_triple,
    add_type,
    create_node,
    extract_id,
    safe_get,
)
from pyBiodatafuse.graph.rdf.nodes.compound import add_tested_compound_node


def add_pubchem_assay_node(
    g: Graph,
    data: dict,
    gene_node: URIRef,
) -> Optional[URIRef]:
    """Create and add a PubChem bioassay node to the RDF graph.

    This function creates an experimental process node representing a PubChem
    bioassay result. It links the gene to the assay and the assay to the
    tested compound.

    :param g: RDF graph to which the assay node will be added.
    :param data: Dictionary containing assay information with keys:
        - pubchem_assay_id: PubChem assay identifier (e.g., 'AID:1340965')
        - assay_type: Type of assay (e.g., 'IC50')
        - outcome: Assay outcome (e.g., 'active', 'inactive')
        - compound_cid: PubChem compound ID
        - compound_name: Name of the compound
        - smiles: SMILES representation
        - inchi: InChI representation
    :param gene_node: URIRef of the gene node to link to this assay.
    :return: URIRef for the created assay node or None if data is invalid.
    """
    # Get and validate assay ID
    pubchem_assay_id = safe_get(data, Cons.PUBCHEM_ASSAY_ID)
    if not pubchem_assay_id:
        return None

    # Extract numeric ID from 'AID:1340965' format
    assay_id = extract_id(pubchem_assay_id, "AID:")
    if not assay_id:
        assay_id = extract_id(pubchem_assay_id, "AID")
    if not assay_id:
        assay_id = str(pubchem_assay_id)

    # Create assay node
    assay_node = create_node(Cons.NODE_URI_PREFIXES["pubchem_assay"], assay_id)
    add_type(g, assay_node, Cons.NODE_TYPES["experimental_process_node"])

    # Add assay label (use assay type if available)
    assay_type = safe_get(data, "assay_type")
    if assay_type:
        add_label(g, assay_node, assay_type)
    else:
        add_label(g, assay_node, f"PubChem Assay {assay_id}")

    # Add identifier
    add_triple(g, assay_node, "dc_identifier", pubchem_assay_id)

    # Link gene to assay (gene is subject of assay)
    add_triple(g, gene_node, "sio_is_subject_of", assay_node)

    # Add outcome
    outcome = safe_get(data, "outcome")
    if outcome:
        add_triple(g, assay_node, "sio_has_output", outcome)

    # Create and link tested compound
    compound_cid = safe_get(data, "compound_cid")
    compound_name = safe_get(data, Cons.PUBCHEM_COMPOUND_NAME)
    smiles = safe_get(data, Cons.PUBCHEM_SMILES)
    inchi = safe_get(data, Cons.PUBCHEM_INCHI)

    if compound_cid or smiles or inchi:
        tested_compound_node = add_tested_compound_node(
            g=g,
            inchi=inchi or "",
            smiles=smiles or "",
            compound_id=compound_cid or "",
            compound_name=compound_name or "",
        )
        if tested_compound_node:
            add_triple(g, assay_node, "sio_has_input", tested_compound_node)

    # Add data source
    add_triple(g, assay_node, "sio_has_source", URIRef(Cons.DATA_SOURCES.get(Cons.PUBCHEM, "")))

    return assay_node


# Keep old function name for backwards compatibility
def add_experimental_process_node(g: Graph, data: dict) -> Optional[URIRef]:
    """Legacy function - creates assay node without gene linkage.

    For new code, use add_pubchem_assay_node() instead.

    :param g: RDF graph.
    :param data: Assay data dictionary.
    :return: URIRef for the created node or None.
    """
    # Get and validate assay ID
    pubchem_assay_id = safe_get(data, Cons.PUBCHEM_ASSAY_ID)
    if not pubchem_assay_id:
        return None

    # Extract numeric ID
    assay_id = extract_id(pubchem_assay_id, "AID:")
    if not assay_id:
        assay_id = extract_id(pubchem_assay_id, "AID")
    if not assay_id:
        assay_id = str(pubchem_assay_id)

    # Create assay node
    assay_node = create_node(Cons.NODE_URI_PREFIXES["pubchem_assay"], assay_id)
    add_type(g, assay_node, Cons.NODE_TYPES["experimental_process_node"])

    assay_type = safe_get(data, "assay_type")
    add_label(g, assay_node, assay_type or f"PubChem Assay {assay_id}")
    add_triple(g, assay_node, "dc_identifier", pubchem_assay_id)

    # Add outcome
    outcome = safe_get(data, "outcome")
    if outcome:
        add_triple(g, assay_node, "sio_has_output", outcome)

    # Create and link tested compound
    compound_cid = safe_get(data, "compound_cid")
    compound_name = safe_get(data, Cons.PUBCHEM_COMPOUND_NAME)
    smiles = safe_get(data, Cons.PUBCHEM_SMILES)
    inchi = safe_get(data, Cons.PUBCHEM_INCHI)

    if compound_cid or smiles or inchi:
        tested_compound_node = add_tested_compound_node(
            g=g,
            inchi=inchi or "",
            smiles=smiles or "",
            compound_id=compound_cid or "",
            compound_name=compound_name or "",
        )
        if tested_compound_node:
            add_triple(g, assay_node, "sio_has_input", tested_compound_node)

    # Add data source
    add_triple(g, assay_node, "sio_has_source", URIRef(Cons.DATA_SOURCES.get(Cons.PUBCHEM, "")))

    return assay_node
