# assay.py

"""Populate a BDF RDF graph with assay data."""


from rdflib import Graph, Literal, URIRef
from rdflib.namespace import DC, RDF, RDFS, XSD

import pyBiodatafuse.constants as Cons
from pyBiodatafuse.graph.rdf.nodes.compound import (  # pylint: disable=no-name-in-module
    add_tested_compound_node,
)
from pyBiodatafuse.graph.rdf.utils import add_data_source_node  # pylint: disable=no-name-in-module


def add_experimental_process_node(
    g: Graph,
    data: dict,
):
    """Create and add an experimental process node to the RDF graph.

    :param g: RDF graph to which the experimental process node will be added.
    :param data: List of dictionaries containing experimental process information.
    :return: URIRef for the created experimental process node.
    """
    pubchem_assay_id = data.get(Cons.PUBCHEM_ASSAY_ID, None)
    if pubchem_assay_id:
        pubchem_assay_iri = Cons.NODE_URI_PREFIXES["pubchem_assay"] + str(pubchem_assay_id).strip(
            "AID"
        )
        assay_type = data["assay_type"]
        inchi = data[Cons.PUBCHEM_INCHI]
        outcome = data["outcome"]
        compound_cid = data["compound_cid"]
        compound_name = data[Cons.PUBCHEM_COMPOUND_NAME]
        smiles = data[Cons.PUBCHEM_SMILES]
        experimental_process_node = URIRef(pubchem_assay_iri)
        g.add(
            (
                experimental_process_node,
                RDF.type,
                URIRef(Cons.NODE_TYPES["experimental_process_node"]),
            )
        )
        g.add((experimental_process_node, RDFS.label, Literal(assay_type, datatype=XSD.string)))
        g.add(
            (
                experimental_process_node,
                DC.identifier,
                Literal(pubchem_assay_id, datatype=XSD.string),
            )
        )
        # has tested compound
        tested_compound_node = add_tested_compound_node(
            g=g, inchi=inchi, smiles=smiles, compound_id=compound_cid, compound_name=compound_name
        )
        g.add(
            (
                experimental_process_node,
                URIRef(Cons.PREDICATES["sio_has_input"]),
                tested_compound_node,
            )
        )
        # has outcome
        g.add(
            (experimental_process_node, URIRef(Cons.PREDICATES["sio_has_output"]), Literal(outcome))
        )
        data_source_node = add_data_source_node(g, Cons.PUBCHEM)
        g.add(
            (experimental_process_node, URIRef(Cons.PREDICATES["sio_has_source"]), data_source_node)
        )
        return experimental_process_node
