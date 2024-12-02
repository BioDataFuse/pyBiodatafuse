# pathway.py


"""Populate a BDF RDF graph with pathway nodes."""


from rdflib import Graph, Literal, URIRef
from rdflib.namespace import RDF, RDFS, XSD

from pyBiodatafuse.constants import NODE_TYPES, PREDICATES


def add_pathway_node(g: Graph, data: dict, source: str):
    """Create and add a pathway node to the RDF graph.

    :param g: RDF graph to which the pathway node will be added.
    :param data: Dictionary containing pathway information.
    :param source: Source of the pathway information (e.g., WikiPathways, Reactome).
    :return: URIRef for the created pathway node.
    """
    pathway_label = data.get("pathway_label", None)
    pathway_id = data.get("pathway_id", None)

    if pathway_id:
        pathway_iri = ""
        iri_source = None
        if source == "WikiPathways":
            pathway_iri = "https://www.wikipathways.org/pathways/" + str(pathway_id.split(":")[1])
            iri_source = "https://www.wikipathways.org/"

        if source == "OpenTargets_reactome":
            pathway_id = pathway_id.split(":")[1]
            pathway_iri = "https://reactome.org/content/detail/" + str(pathway_id)
            iri_source = "https://reactome.org/"
            source = "Reactome"
        if source == "MINERVA":
            pathway_iri = "https://minerva-net.lcsb.uni.lu/api/" + str(pathway_id)
            iri_source = "https://minerva-net.lcsb.uni.lu/"
        pathway_node = URIRef(pathway_iri)
        g.add((pathway_node, RDF.type, URIRef(NODE_TYPES["pathway_node"])))
        g.add((pathway_node, RDFS.label, Literal(pathway_label, datatype=XSD.string)))
        if iri_source:
            g.add((pathway_node, URIRef(PREDICATES["sio_has_source"]), URIRef(iri_source)))
            g.add((URIRef(iri_source), RDFS.label, Literal(source, datatype=XSD.string)))
            # data_source_node = add_data_source_node(g, source)
            g.add((URIRef(iri_source), RDF.type, URIRef(NODE_TYPES["data_source_node"])))
        return pathway_node
