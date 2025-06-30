# go_terms.py

"""Populate a BDF RDF graph with gene ontology terms related to genes."""


from typing import Optional

from bioregistry import get_iri
from rdflib import Graph, Literal, URIRef
from rdflib.namespace import RDFS, XSD

import pyBiodatafuse.constants as Cons
from pyBiodatafuse.graph.rdf.utils import add_data_source_node


def add_go_cpf(g: Graph, process_data: dict) -> Optional[URIRef]:
    """Create and add a Gene Ontology (GO) node to the RDF graph.

    :param g: RDF graph to which the GO node will be added.
    :param process_data: Dictionary containing Gene Ontology information.
    :return: URIRef for the created GO node.
    """
    curie = process_data[Cons.GO_ID]
    if curie:
        iri = get_iri(curie)
        if not iri:
            return None
        label = process_data[Cons.GO_NAME]
        go_type = (
            URIRef(Cons.GO_TYPES[process_data[Cons.GO_TYPE]])
            if process_data[Cons.GO_TYPE] in Cons.GO_TYPES
            else None
        )
        go_cpf = URIRef(iri)
        g.add((go_cpf, RDFS.label, Literal(label, datatype=XSD.string)))
        if go_type:
            g.add((go_cpf, RDFS.subClassOf, go_type))
        data_source_node = add_data_source_node(g, Cons.OPENTARGETS_REACTOME_COL)
        g.add((go_cpf, URIRef(Cons.PREDICATES["sio_has_source"]), data_source_node))
        return go_cpf
    else:
        return None
