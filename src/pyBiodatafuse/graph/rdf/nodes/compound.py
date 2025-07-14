# compound.py

"""Populate a BDF RDF graph with compound data."""


from typing import Optional

from rdflib import Graph, Literal, URIRef
from rdflib.namespace import OWL, RDF, RDFS, XSD

import pyBiodatafuse.constants as Cons


def get_compound_node(g: Graph, row) -> tuple:
    """Create and add a compound node and associated nodes to the RDF graph.

    :param g: (Graph): RDF graph to which the compound nodes are added.
    :param row: (pd.Series): Data row containing gene information.

    :return: URIRef for the gene node.
    """
    target = row.get("target", None)
    target_source = row.get("target.source", None)
    identifier = row.get("identifier", None)
    identifier_source = row.get("identifier.source", None)

    if identifier_source and identifier_source == "PubChem Compound":
        if identifier:
            compound_node = URIRef(Cons.SOURCE_NAMESPACES["PubChem Compound"] + identifier)
            g.add((compound_node, RDF.type, URIRef(Cons.NODE_TYPES["compound_node"])))
            g.add((compound_node, RDFS.label, Literal(row["identifier"], datatype=XSD.string)))
    if target_source and target_source in Cons.SOURCE_NAMESPACES:
        if target:
            target_node = URIRef(Cons.SOURCE_NAMESPACES[target_source] + target)
            g.add((target_node, RDF.type, URIRef(Cons.NODE_TYPES["compound_node"])))
            g.add((target_node, RDFS.label, Literal(row["target"], datatype=XSD.string)))
            g.add((compound_node, URIRef(Cons.PREDICATES["sameAs"]), target_node))
    return compound_node


def add_associated_compound_node(g: Graph, compound: dict, gene_node: URIRef) -> URIRef:
    """Create and add a compound node to the RDF graph.

    :param g: RDF graph to which the compound node will be added.
    :param compound: Dictionary containing compound data.
    :param gene_node: URIRef of the gene node associated with the compound.
    :return: URIRef for the created compound node.
    """
    chembl_id = compound.get("chembl_id")
    compoundbank_id = compound.get("compoundbank_id")
    compound_name = compound.get("compound_name")
    compound_cid = compound.get("compound_cid", "")
    is_approved = compound.get("is_approved")
    clinical_trial_phase = compound.get("clinical_trial_phase")
    # moa = compound.get("mechanisms_of_action", None)
    if chembl_id:
        chembl_id = chembl_id.split(f"{Cons.CHEMBL}:")[1]
        compound_node = URIRef(f"https://www.ebi.ac.uk/chembl/compound_report_card/{chembl_id}")
        g.add((compound_node, RDFS.label, Literal(compound_name, datatype=XSD.string)))

        if is_approved:
            g.add((compound_node, RDF.type, URIRef(Cons.NODE_TYPES["approved_compound"])))
        else:
            g.add((compound_node, RDF.type, URIRef(Cons.NODE_TYPES["tested_compound_node"])))

        # Clinical trial phase
        if clinical_trial_phase:
            clinical_phase_iri = Cons.CLINICAL_PHASES.get(str(clinical_trial_phase))
            if clinical_phase_iri:
                g.add(
                    (
                        compound_node,
                        URIRef(Cons.PREDICATES["sio_has_value"]),
                        URIRef(clinical_phase_iri),
                    )
                )

        # compoundBank and PubChem identifiers
        for source_id in [compound_cid]:
            if source_id:
                iri = (
                    f"https://go.compoundbank.com/compounds/{source_id}"
                    if source_id == compoundbank_id
                    else (
                        f"https://pubchem.ncbi.nlm.nih.gov/compound/{compound_cid.split(':')[1]}"
                        if ":" in compound_cid
                        else f"https://pubchem.ncbi.nlm.nih.gov/compound/{compound_cid}"
                    )
                )
                id_node = URIRef(iri)
                g.add((compound_node, OWL.sameAs, id_node))
                g.add((id_node, OWL.sameAs, compound_node))
                g.add((id_node, RDF.type, URIRef(Cons.NODE_TYPES["tested_compound_node"])))
        return compound_node
    else:
        return None


def add_transporter_inhibitor_node(
    g: Graph, gene_node, inhibitor_data: dict, base_uri: str
) -> None:
    """Add a transporter inhibitor node to the RDF graph.

    :param g: RDF graph to which the transporter inhibitor node is added.
    :param gene_node: URIRef of the gene node associated with the transporter inhibitor.
    :param inhibitor_data: Dictionary containing transporter inhibitor data (e.g., "compound_cid", "inchikey").
    :param base_uri: The base URI for constructing nodes in the graph.
    """
    compound_cid = inhibitor_data.get(Cons.PUBCHEM_COMPOUND_CID)
    compound_name = inhibitor_data.get(Cons.PUBCHEM_COMPOUND_NAME)
    inchikey = inhibitor_data.get(Cons.PUBCHEM_INCHI)
    smiles = inhibitor_data.get(Cons.PUBCHEM_SMILES)
    molmedb_id = inhibitor_data.get("molmedb_id")
    uniprot_id = inhibitor_data.get(Cons.MOLMEDB_UNIPROT_TREMBL_ID)

    if compound_cid:
        # Create compound node
        compound_node = URIRef(Cons.NODE_URI_PREFIXES[Cons.PUBCHEM] + compound_cid)
        g.add((compound_node, RDFS.label, Literal(compound_name, datatype=XSD.string)))
        g.add(
            (
                compound_node,
                URIRef(Cons.PREDICATES["chebi_inchi"]),
                Literal(inchikey, datatype=XSD.string),
            )
        )
        g.add((compound_node, URIRef(Cons.PREDICATES["chebi_smiles"]), Literal(smiles)))

        # Link to other identifiers
        if molmedb_id:
            g.add(
                (
                    compound_node,
                    OWL.sameAs,
                    URIRef(Cons.NODE_URI_PREFIXES[Cons.MOLMEDB] + molmedb_id),
                )
            )

        g.add((uniprot_id, Cons.PREDICATES["translation_of"], gene_node))
        uniprot_node = URIRef(f"https://www.uniprot.org/uniprotkb/{uniprot_id}")
        g.add((uniprot_node, RDFS.label, Literal(uniprot_id, datatype=XSD.string)))
        g.add((gene_node, Cons.PREDICATES["translates_to"], uniprot_node))
        g.add((uniprot_node, Cons.PREDICATES["translation_of"], gene_node))
        # Add inhibition relationship
        inhibition_node = URIRef(f"{base_uri}inhibition/{uniprot_id}_{compound_cid}")
        g.add(
            (
                uniprot_node,
                URIRef(Cons.PREDICATES["sio_is_part_of"]),
                inhibition_node,
            )
        )
        g.add((compound_node, URIRef(Cons.PREDICATES["sio_is_part_of"]), inhibition_node))
        g.add((inhibition_node, RDF.type, URIRef("https://purl.obolibrary.org/GO_0032410")))


def add_inhibitor_transporter_node(
    g: Graph, compound_node: URIRef, inhibitor_data: dict, base_uri: str
) -> None:
    """Add a inhibitor transporter inhibitor node to the RDF graph.

    :param g: RDF graph to which the transporter inhibitor node is added.
    :param compound_node: URIRef of the compound node associated with the transporter inhibitor.
    :param inhibitor_data: Dictionary containing transporter inhibitor data (e.g., "compound_cid", "inchikey").
    :param base_uri: The base URI for constructing nodes in the graph.
    :raises ValueError: If the "hgnc_symbol" key in `inhibitor_data` is None.
    """
    uniprot_id = inhibitor_data.get("uniprot_trembl_id")
    hgnc = inhibitor_data.get("hgnc_symbol")
    compound_cid = None
    compoundbank_id = inhibitor_data.get("compoundbank_id")
    if compound_node:
        compound_cid = g.value(compound_node, RDFS.label)
    g.add(
        (
            compound_node,
            OWL.sameAs,
            URIRef(f"https://www.compoundbank.ca/compounds/{compoundbank_id}"),
        )
    )
    # Add inhibition relationship
    inhibition_node = URIRef(f"{base_uri}inhibition/{uniprot_id}_{compound_cid}")
    g.add(
        (
            URIRef(f"https://www.uniprot.org/uniprotkb/{uniprot_id}"),
            URIRef(Cons.PREDICATES["sio_is_part_of"]),
            inhibition_node,
        )
    )
    g.add((compound_node, URIRef(Cons.PREDICATES["sio_is_part_of"]), inhibition_node))
    g.add((inhibition_node, RDF.type, URIRef("https://purl.obolibrary.org/GO_0032410")))
    if hgnc:
        inhibited_gene_node = URIRef(Cons.NAMESPACE_BINDINGS["hgnc"] + hgnc)
    else:
        raise ValueError("hgnc is None and cannot be used to create a URIRef.")
    g.add((inhibited_gene_node, RDF.type, URIRef(Cons.NODE_TYPES["gene_node"])))
    g.add((inhibited_gene_node, RDFS.label, Literal(hgnc, datatype=XSD.string)))
    g.add(
        (
            inhibited_gene_node,
            URIRef(Cons.PREDICATES["sio_is_part_of"]),
            inhibition_node,
        )
    )


def add_tested_compound_node(
    g: Graph,
    inchi: str,
    smiles: str,
    compound_id: str,
    compound_name: str,
):
    """Create and add a tested compound node to the RDF graph.

    :param g: RDF graph to which the tested compound node will be added.
    :param inchi: inchi identifier of the compound.
    :param smiles: smiles identifier of the compound.
    :param compound_id: Compound ID.
    :param compound_name: Name of the compound.
    :return: URIRef for the created tested compound node.
    """
    uri_cas = Cons.NODE_URI_PREFIXES[Cons.PUBCHEM] + str(compound_id).strip("CID")
    tested_compound_node = URIRef((uri_cas))
    g.add((tested_compound_node, RDF.type, URIRef(Cons.NODE_TYPES["tested_compound_node"])))
    g.add((tested_compound_node, RDF.type, URIRef(Cons.NODE_TYPES["compound_node"])))
    g.add((tested_compound_node, RDFS.label, Literal(compound_name, datatype=XSD.string)))
    g.add(
        (
            tested_compound_node,
            URIRef(Cons.PREDICATES["chebi_inchi"]),
            Literal(inchi, datatype=XSD.string),
        )
    )
    g.add(
        (
            tested_compound_node,
            URIRef(Cons.PREDICATES["chebi_smiles"]),
            Literal(smiles, datatype=XSD.string),
        )
    )
    g.add(
        (
            tested_compound_node,
            URIRef(Cons.PREDICATES["cheminf_compound_id"]),
            Literal(compound_id, datatype=XSD.string),
        )
    )
    return tested_compound_node
