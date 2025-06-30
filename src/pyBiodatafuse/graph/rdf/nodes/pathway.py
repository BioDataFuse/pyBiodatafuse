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


def add_molecular_pathway_node(g, el, identifier, id_number):
    """
    Add a molecular pathway node and related triples to the RDF graph.

    :param g: RDF graph.
    :param el: Dictionary with pathway data.
    :param identifier: RDF node for gene or compound.
    :param id_number: Row identifier.
    """
    from pyBiodatafuse.graph.rdf.nodes.gene import get_gene_node

    pathway_id = el.get("pathway_id", None)
    pathway_label = el.get("pathway_label", None)
    target_gene = el.get("targetGene", None)
    target_protein = el.get("targetProtein", None)
    target_metabolite = el.get("targetMetabolite", None)
    mimtype = el.get("mimtype", None)
    rhea_id = el.get("rhea_id", None)
    if mimtype:
        mimtype_node = URIRef("http://vocabularies.wikipathways.org/wp#" + mimtype)
        g.add((mimtype_node, RDF.type, URIRef(Cons.NODE_TYPES["interaction"])))
        g.add((mimtype_node, RDFS.label, Literal(mimtype, datatype=XSD.string)))
        mim_node = URIRef(g.new_uris["interaction"] + mimtype + "_" + id_number)
        g.add((mim_node, RDF.type, URIRef(Cons.NODE_TYPES["interaction"])))
        g.add((mim_node, RDFS.label, Literal(mimtype + "_" + id_number, datatype=XSD.string)))
        g.add((identifier, URIRef(Cons.PREDICATES["sio_is_part_of"]), mim_node))
        if pathway_id and pathway_label:
            pathway_node = URIRef(
                "https://www.wikipathways.org/pathways/" + pathway_id.split(":")[1]
            )
            g.add((pathway_node, URIRef(Cons.PREDICATES["sio_has_part"]), mimtype_node))
            g.add((identifier, URIRef(Cons.PREDICATES["sio_is_part_of"]), pathway_node))
            g.add((pathway_node, URIRef(Cons.PREDICATES["sio_has_part"]), identifier))
            g.add((pathway_node, URIRef(Cons.PREDICATES["sio_has_part"]), mim_node))

        if target_gene:
            target_gene_node = get_gene_node(
                g,
                {
                    "target.source": "Ensembl",
                    "target": target_gene,
                    "identifier": target_gene,
                },
            )
            if target_gene_node:
                g.add(
                    (
                        identifier,
                        URIRef(g.new_uris["interaction"] + mimtype),
                        target_gene_node,
                    )
                )
                g.add(
                    (
                        target_gene_node,
                        URIRef("http://vocabularies.wikipathways.org/wp#" + mimtype),
                        identifier,
                    )
                )
                g.add((target_gene_node, URIRef(Cons.PREDICATES["sio_is_part_of"]), mim_node))

        if target_protein:
            target_protein_node = URIRef(Cons.BASE_URLS_DBS["uniprot"] + target_protein)
            g.add((target_protein_node, RDF.type, URIRef(Cons.NODE_TYPES["protein_node"])))
            g.add(
                (
                    identifier,
                    URIRef("http://vocabularies.wikipathways.org/wp#" + mimtype),
                    target_protein_node,
                )
            )
            g.add(
                (
                    target_protein_node,
                    URIRef("http://vocabularies.wikipathways.org/wp#" + mimtype),
                    identifier,
                )
            )
            g.add((target_protein_node, URIRef(Cons.PREDICATES["sio_is_part_of"]), mim_node))

        if target_metabolite:
            target_metabolite_node = URIRef(Cons.NODE_TYPES["compound_node"] + target_metabolite)
            g.add((target_metabolite_node, RDF.type, URIRef(Cons.NODE_TYPES["compound_node"])))
            g.add(
                (
                    identifier,
                    URIRef("http://vocabularies.wikipathways.org/wp#" + mimtype),
                    target_metabolite_node,
                )
            )
            g.add(
                (
                    target_metabolite_node,
                    URIRef("http://vocabularies.wikipathways.org/wp#" + mimtype),
                    identifier,
                )
            )
            g.add((identifier, URIRef(Cons.PREDICATES["sio_is_part_of"]), mim_node))
            g.add(
                (
                    target_metabolite_node,
                    URIRef(Cons.PREDICATES["sio_is_part_of"]),
                    mim_node,
                )
            )

        if rhea_id:
            g.add((mim_node, URIRef(Cons.PREDICATES["sameAs"]), URIRef(rhea_id)))
