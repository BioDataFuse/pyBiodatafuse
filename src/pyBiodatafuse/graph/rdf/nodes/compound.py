# compound.py

"""Populate a BDF RDF graph with compound data."""

from typing import Optional

from rdflib import Graph, Literal, URIRef
from rdflib.namespace import RDFS, XSD

import pyBiodatafuse.constants as Cons
from pyBiodatafuse.graph.rdf.nodes.base import (
    add_label,
    add_same_as,
    add_triple,
    add_type,
    create_node,
    extract_id,
    get_uri,
    link_associated_with,
    link_has_part,
    safe_get,
)


def get_compound_node(g: Graph, row) -> Optional[URIRef]:
    """Create and add a compound node and associated nodes to the RDF graph.

    :param g: RDF graph to which the compound nodes are added.
    :param row: Data row containing compound information.
    :return: URIRef for the compound node.
    """
    target = safe_get(row, "target")
    target_source = safe_get(row, "target.source")
    identifier = safe_get(row, "identifier")
    identifier_source = safe_get(row, "identifier.source")
    compound_node = None
    if identifier_source and identifier:
        compound_node = get_uri(identifier_source, identifier)
        add_type(g, compound_node, Cons.NODE_TYPES["compound_node"])
        add_label(g, compound_node, identifier)

    if target_source and target_source in Cons.SOURCE_NAMESPACES and target:
        target_node = get_uri(target_source, target)
        add_type(g, target_node, Cons.NODE_TYPES["compound_node"])
        add_label(g, target_node, target)
        if compound_node:
            add_same_as(g, compound_node, target_node, bidirectional=False)
    return compound_node


def add_associated_compound_node(g: Graph, compound: dict, gene_node: URIRef) -> Optional[URIRef]:
    """Create and add a compound node to the RDF graph.

    :param g: RDF graph to which the compound node will be added.
    :param compound: Dictionary containing compound data.
    :param gene_node: URIRef of the gene node associated with the compound.
    :return: URIRef for the created compound node or None.
    """
    chembl_id = safe_get(compound, "chembl_id")
    compoundbank_id = safe_get(compound, "compoundbank_id")
    compound_name = safe_get(compound, "compound_name")
    compound_cid = safe_get(compound, "compound_cid", "")
    is_approved = safe_get(compound, "is_approved")
    clinical_trial_phase = safe_get(compound, "clinical_trial_phase")

    if not chembl_id:
        return None

    # Extract ID from CHEMBL prefix
    chembl_id = extract_id(chembl_id, f"{Cons.CHEMBL}:")
    compound_node = create_node("https://www.ebi.ac.uk/chembl/compound_report_card/", chembl_id)
    add_label(g, compound_node, compound_name)

    # Set type based on approval status
    if is_approved:
        add_type(g, compound_node, Cons.NODE_TYPES["approved_compound"])
    else:
        add_type(g, compound_node, Cons.NODE_TYPES["tested_compound_node"])

    # Clinical trial phase
    if clinical_trial_phase:
        clinical_phase_iri = Cons.CLINICAL_PHASES.get(str(clinical_trial_phase))
        if clinical_phase_iri:
            add_triple(g, compound_node, "sio_has_value", URIRef(clinical_phase_iri))

    # PubChem identifiers
    if compound_cid:
        pubchem_id = extract_id(compound_cid)
        iri = (
            f"https://go.compoundbank.com/compounds/{compound_cid}"
            if compound_cid == compoundbank_id
            else f"https://pubchem.ncbi.nlm.nih.gov/compound/{pubchem_id}"
        )
        id_node = URIRef(iri)
        add_same_as(g, compound_node, id_node)
        add_type(g, id_node, Cons.NODE_TYPES["tested_compound_node"])

    return compound_node


def add_transporter_inhibitor_node(
    g: Graph, gene_node: URIRef, inhibitor_data: dict, base_uri: str
) -> None:
    """Add a transporter inhibitor node to the RDF graph.

    :param g: RDF graph to which the transporter inhibitor node is added.
    :param gene_node: URIRef of the gene node associated with the transporter inhibitor.
    :param inhibitor_data: Dictionary containing transporter inhibitor data.
    :param base_uri: The base URI for constructing nodes in the graph.
    """
    compound_cid = safe_get(inhibitor_data, Cons.PUBCHEM_COMPOUND_CID)
    compound_name = safe_get(inhibitor_data, Cons.PUBCHEM_COMPOUND_NAME)
    inchikey = safe_get(inhibitor_data, Cons.PUBCHEM_INCHI)
    smiles = safe_get(inhibitor_data, Cons.PUBCHEM_SMILES)
    molmedb_id = safe_get(inhibitor_data, "molmedb_id")
    uniprot_id = safe_get(inhibitor_data, Cons.MOLMEDB_UNIPROT_TREMBL_ID)

    if not compound_cid:
        return

    # Create compound node
    compound_node = get_uri(Cons.PUBCHEM, compound_cid)
    add_label(g, compound_node, compound_name)
    add_triple(g, compound_node, "chebi_inchi", Literal(inchikey, datatype=XSD.string))
    add_triple(g, compound_node, "chebi_smiles", Literal(smiles, datatype=XSD.string))

    # Link to MolMeDB identifier
    if molmedb_id:
        molmedb_node = get_uri(Cons.MOLMEDB, molmedb_id)
        add_same_as(g, compound_node, molmedb_node, bidirectional=False)

    # Create UniProt protein node
    if uniprot_id:
        uniprot_node = create_node("https://www.uniprot.org/uniprotkb/", uniprot_id)
        add_label(g, uniprot_node, uniprot_id)
        add_triple(g, gene_node, "translates_to", uniprot_node)
        add_triple(g, uniprot_node, "translation_of", gene_node)

        # Add inhibition relationship
        inhibition_node = create_node(f"{base_uri}inhibition/", f"{uniprot_id}_{compound_cid}")
        add_type(g, inhibition_node, "https://purl.obolibrary.org/GO_0032410")
        link_has_part(g, inhibition_node, uniprot_node, bidirectional=True)
        link_has_part(g, inhibition_node, compound_node, bidirectional=True)


def add_inhibitor_transporter_node(
    g: Graph, compound_node: Optional[URIRef], inhibitor_data: dict, base_uri: str
) -> None:
    """Add an inhibitor transporter node to the RDF graph.

    :param g: RDF graph to which the transporter inhibitor node is added.
    :param compound_node: URIRef of the compound node associated with the transporter inhibitor.
    :param inhibitor_data: Dictionary containing transporter inhibitor data.
    :param base_uri: The base URI for constructing nodes in the graph.
    :raises ValueError: If the "hgnc_symbol" key in `inhibitor_data` is None.
    """
    uniprot_id = safe_get(inhibitor_data, "uniprot_trembl_id")
    hgnc = safe_get(inhibitor_data, "hgnc_symbol")
    compoundbank_id = safe_get(inhibitor_data, "compoundbank_id")

    # Get compound ID from node label
    compound_cid = None
    if compound_node:
        compound_cid = g.value(compound_node, RDFS.label)

    # Link to CompoundBank
    if compoundbank_id:
        compoundbank_node = create_node("https://www.compoundbank.ca/compounds/", compoundbank_id)
        add_same_as(g, compound_node, compoundbank_node, bidirectional=False)

    # Create inhibited gene node first
    if not hgnc:
        raise ValueError("hgnc is None and cannot be used to create a URIRef.")

    inhibited_gene_node = get_uri("hgnc", hgnc)
    if inhibited_gene_node:
        add_type(g, inhibited_gene_node, Cons.NODE_TYPES["gene_node"])
        add_label(g, inhibited_gene_node, hgnc)

    # Add inhibition relationship
    if uniprot_id and compound_cid:
        inhibition_node = create_node(f"{base_uri}inhibition/", f"{uniprot_id}_{compound_cid}")
        add_type(g, inhibition_node, "https://purl.obolibrary.org/GO_0032410")

        uniprot_node = create_node("https://www.uniprot.org/uniprotkb/", uniprot_id)
        link_has_part(g, inhibition_node, uniprot_node, bidirectional=True)
        link_has_part(g, inhibition_node, compound_node, bidirectional=True)
        if inhibited_gene_node:
            link_has_part(g, inhibition_node, inhibited_gene_node, bidirectional=True)


def add_tested_compound_node(
    g: Graph,
    inchi: str,
    smiles: str,
    compound_id: str,
    compound_name: str,
) -> URIRef:
    """Create and add a tested compound node to the RDF graph.

    :param g: RDF graph to which the tested compound node will be added.
    :param inchi: inchi identifier of the compound.
    :param smiles: smiles identifier of the compound.
    :param compound_id: Compound ID.
    :param compound_name: Name of the compound.
    :return: URIRef for the created tested compound node.
    """
    clean_id = str(compound_id).strip("CID")
    tested_compound_node = get_uri(Cons.PUBCHEM, clean_id)
    if not tested_compound_node:
        tested_compound_node = create_node(Cons.NODE_URI_PREFIXES[Cons.PUBCHEM], clean_id)

    add_type(g, tested_compound_node, Cons.NODE_TYPES["tested_compound_node"])
    add_type(g, tested_compound_node, Cons.NODE_TYPES["compound_node"])
    add_label(g, tested_compound_node, compound_name)
    add_triple(g, tested_compound_node, "chebi_inchi", Literal(inchi, datatype=XSD.string))
    add_triple(g, tested_compound_node, "chebi_smiles", Literal(smiles, datatype=XSD.string))
    add_triple(
        g, tested_compound_node, "cheminf_compound_id", Literal(compound_id, datatype=XSD.string)
    )

    return tested_compound_node


def add_compoundwiki_annotations(
    g: Graph,
    target_node: URIRef,
    compoundwiki_data: list,
) -> None:
    """Add CompoundWiki annotations to a target node (gene or protein).

    :param g: RDF graph to which the CompoundWiki annotations will be added.
    :param target_node: URIRef of the target node (gene or protein from PubChem assays).
    :param compoundwiki_data: List of CompoundWiki annotation dictionaries.
    """
    if not compoundwiki_data:
        return

    for annotation_list in compoundwiki_data:
        if not isinstance(annotation_list, list):
            continue

        for annotation in annotation_list:
            if not isinstance(annotation, dict):
                continue

            # Get compound identifiers
            pubchem_cid = safe_get(annotation, Cons.COMPOUNDWIKI_PUBCHEM_ID)
            compound_label = safe_get(annotation, Cons.COMPOUNDWIKI_LABEL)
            inchi = safe_get(annotation, Cons.COMPOUNDWIKI_INCHI)
            smiles = safe_get(annotation, Cons.COMPOUNDWIKI_SMILES)

            if not pubchem_cid:
                continue

            # Create compound node
            pubchem_cid_clean = str(pubchem_cid).strip()
            compound_node = create_node(
                "https://pubchem.ncbi.nlm.nih.gov/compound/", pubchem_cid_clean
            )

            # Add compound node with basic properties
            add_type(g, compound_node, Cons.NODE_TYPES["compound_node"])
            add_label(g, compound_node, compound_label)

            # Add chemical structure identifiers
            if inchi:
                add_triple(g, compound_node, "chebi_inchi", Literal(inchi, datatype=XSD.string))
            if smiles:
                add_triple(g, compound_node, "chebi_smiles", Literal(smiles, datatype=XSD.string))

            # Add cross-references to other databases
            chebi_id = safe_get(annotation, Cons.COMPOUNDWIKI_CHEBI_ID)
            if chebi_id:
                chebi_node = create_node(
                    "http://purl.obolibrary.org/obo/CHEBI_", str(chebi_id).strip()
                )
                add_same_as(g, compound_node, chebi_node, bidirectional=False)

            chembl_id = safe_get(annotation, Cons.COMPOUNDWIKI_CHEMBL_ID)
            if chembl_id:
                chembl_node = create_node(
                    "https://www.ebi.ac.uk/chembl/compound_report_card/", str(chembl_id).strip()
                )
                add_same_as(g, compound_node, chembl_node, bidirectional=False)

            kegg_id = safe_get(annotation, Cons.COMPOUNDWIKI_KEGG_ID)
            if kegg_id:
                kegg_node = create_node("https://www.genome.jp/entry/", str(kegg_id).strip())
                add_same_as(g, compound_node, kegg_node, bidirectional=False)

            wikidata_id = safe_get(annotation, Cons.COMPOUNDWIKI_WIKIDATA_ID)
            if wikidata_id:
                wikidata_node = create_node(
                    "http://www.wikidata.org/entity/", str(wikidata_id).strip()
                )
                add_same_as(g, compound_node, wikidata_node, bidirectional=False)

            # Add chemical formula if available
            formula = safe_get(annotation, Cons.COMPOUNDWIKI_FORMULA)
            if formula:
                add_triple(
                    g, compound_node, "molecular_formula", Literal(formula, datatype=XSD.string)
                )

            # Link compound to target node (association)
            link_associated_with(g, target_node, compound_node)
