# go_terms.py

"""Populate a BDF RDF graph with gene ontology terms related to genes."""


from bioregistry import get_iri
from rdflib import Graph, Literal, URIRef
from rdflib.namespace import RDFS, XSD

from pyBiodatafuse.constants import GO_TYPES, PREDICATES
from pyBiodatafuse.graph.rdf.utils import add_data_source_node


def add_go_cpf(g: Graph, process_data: dict) -> URIRef:
    """Create and add a Gene Ontology (GO) node to the RDF graph.

    :param g: RDF graph to which the GO node will be added.
    :param process_data: Dictionary containing Gene Ontology information.
    :return: URIRef for the created GO node.
    """
    curie = process_data["go_id"]
    if curie:
        iri = get_iri(curie)
        label = process_data["go_name"]
        go_type = (
            URIRef(GO_TYPES[process_data["go_type"]])
            if process_data["go_type"] in GO_TYPES
            else None
        )
        go_cpf = URIRef(iri)
        g.add((go_cpf, RDFS.label, Literal(label, datatype=XSD.string)))
        if go_type:
            g.add((go_cpf, RDFS.subClassOf, go_type))
        data_source_node = add_data_source_node(g, "OpenTargets_reactome")
        g.add((go_cpf, URIRef(PREDICATES["sio_has_source"]), data_source_node))
        return go_cpf
    else:
        return None
