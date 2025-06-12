# pathway.py


"""Populate a BDF RDF graph with pathway nodes."""


from rdflib import Graph, Literal, URIRef
from rdflib.namespace import RDF, RDFS, XSD

import pyBiodatafuse.constants as Cons


def add_pathway_node(g: Graph, data: dict, source: str):
    """Create and add a pathway node to the RDF graph.

    :param g: RDF graph to which the pathway node will be added.
    :param data: Dictionary containing pathway information.
    :param source: Source of the pathway information (e.g., WikiPathways, Reactome).
    :return: URIRef for the created pathway node.
    """
    pathway_label = data.get(Cons.PATHWAY_LABEL, None)
    pathway_id = data.get(Cons.PATHWAY_ID, None)

    if pathway_id:
        pathway_iri = ""
        iri_source = None
        if source == Cons.WIKIPATHWAYS_PATHWAY_COL:
            pathway_iri = Cons.NODE_URI_PREFIXES[Cons.WIKIPATHWAYS] + str(pathway_id.split(":")[1])
            iri_source = Cons.SOURCE_NAMESPACES[Cons.WIKIPATHWAYS]

        if source == Cons.OPENTARGETS_REACTOME_COL:
            pathway_id = pathway_id.split(":")[1]
            pathway_iri = Cons.NODE_URI_PREFIXES[Cons.REACTOME] + str(pathway_id)
            iri_source = Cons.SOURCE_NAMESPACES[Cons.REACTOME]
            source = Cons.REACTOME
        if source == Cons.MINERVA_PATHWAY_COL:
            pathway_iri = Cons.NODE_URI_PREFIXES[Cons.MINERVA] + str(pathway_id)
            iri_source = Cons.SOURCE_NAMESPACES["Minerva"]

        pathway_node = URIRef(pathway_iri)
        g.add((pathway_node, RDF.type, URIRef(Cons.NODE_TYPES["pathway_node"])))
        g.add((pathway_node, RDFS.label, Literal(pathway_label, datatype=XSD.string)))
        if iri_source:
            g.add((pathway_node, URIRef(Cons.PREDICATES["sio_has_source"]), URIRef(iri_source)))
            g.add((URIRef(iri_source), RDFS.label, Literal(source, datatype=XSD.string)))
            # data_source_node = add_data_source_node(g, source)
            g.add((URIRef(iri_source), RDF.type, URIRef(Cons.NODE_TYPES["data_source_node"])))
        return pathway_node
