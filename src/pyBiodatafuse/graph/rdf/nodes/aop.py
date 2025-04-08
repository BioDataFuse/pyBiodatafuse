# aop.py

"""Populate a BDF RDF graph with AOP data."""


from rdflib import Graph, Literal, URIRef
from rdflib.namespace import OWL, RDF, RDFS, XSD

from pyBiodatafuse.constants import CLINICAL_PHASES, MOAS, NODE_TYPES, PREDICATES
from pyBiodatafuse.graph.rdf.nodes.compound import add_compound_gene_node


def add_aop_compound_node(args):
    return None


def add_aop_gene_node(args):
    return None
