# compound.py

"""Populate a BDF RDF graph with compound data."""


from typing import Optional

from rdflib import Graph, Literal, URIRef
from rdflib.namespace import OWL, RDF, RDFS, XSD

import pyBiodatafuse.constants as Cons


def add_compound_node(g: Graph, compound: dict, gene_node: URIRef) -> Optional[URIRef]:
    """Create and add a compound node to the RDF graph.

    :param g: RDF graph to which the compound node will be added.
    :param compound: Dictionary containing compound data.
    :param gene_node: URIRef of the gene node associated with the compound.
    :return: URIRef for the created compound node.
    """
    chembl_id = compound.get(Cons.CHEMBL_ID)
    compound_name = compound.get(Cons.PUBCHEM_COMPOUND_NAME)
    compound_cid = compound.get(Cons.PUBCHEM_COMPOUND_CID)
    is_approved = compound.get("is_approved")
    clinical_trial_phase = compound.get("clincal_trial_phase")

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
                iri = Cons.NODE_URI_PREFIXES[Cons.PUBCHEM] + str(source_id.split(":")[1])
                id_node = URIRef(iri)
                g.add((compound_node, OWL.sameAs, id_node))
                g.add((id_node, OWL.sameAs, compound_node))
                g.add((id_node, RDF.type, URIRef(Cons.NODE_TYPES["tested_compound_node"])))
        return compound_node
    else:
        return None


def add_transporter_inhibitor_node(g: Graph, inhibitor_data: dict, base_uri: str) -> None:
    """Add a transporter inhibitor node to the RDF graph.

    :param g: RDF graph to which the transporter inhibitor node is added.
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
        # Add inhibition relationship
        if uniprot_id:
            inhibition_node = URIRef(f"{base_uri}inhibition/{uniprot_id}_{compound_cid}")
            g.add(
                (
                    URIRef(Cons.NODE_URI_PREFIXES["uniprot"] + uniprot_id),
                    URIRef(Cons.PREDICATES["sio_is_part_of"]),
                    inhibition_node,
                )
            )
            g.add((compound_node, URIRef(Cons.PREDICATES["sio_is_part_of"]), inhibition_node))
            g.add((inhibition_node, RDF.type, URIRef("https://purl.obolibrary.org/GO_0032410")))


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
