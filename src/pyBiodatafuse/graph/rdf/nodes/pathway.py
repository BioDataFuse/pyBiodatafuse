# pathway.py

"""Populate a BDF RDF graph with pathway nodes."""

from typing import Optional

from rdflib import Graph, URIRef

import pyBiodatafuse.constants as Cons
from pyBiodatafuse.graph.rdf.nodes.base import (
    add_label,
    add_same_as,
    add_triple,
    add_type,
    create_node,
    extract_id,
    link_has_part,
    link_has_source,
    safe_get,
)


def add_pathway_node(g: Graph, data: dict, source: str) -> Optional[URIRef]:
    """Create and add a pathway node to the RDF graph.

    :param g: RDF graph to which the pathway node will be added.
    :param data: Dictionary containing pathway information.
    :param source: Source of the pathway information (e.g., WikiPathways, Reactome).
    :return: URIRef for the created pathway node.
    """
    pathway_label = safe_get(data, Cons.PATHWAY_LABEL)
    pathway_id = safe_get(data, Cons.PATHWAY_ID)

    if not pathway_id:
        return None

    pathway_node = None
    iri_source = None

    if source == Cons.WIKIPATHWAYS_PATHWAY_COL:
        clean_id = extract_id(pathway_id)
        pathway_node = create_node(Cons.NODE_URI_PREFIXES[Cons.WIKIPATHWAYS], clean_id)
        iri_source = Cons.SOURCE_NAMESPACES[Cons.WIKIPATHWAYS]

    elif source == Cons.OPENTARGETS_REACTOME_COL:
        clean_id = extract_id(pathway_id)
        pathway_node = create_node(Cons.NODE_URI_PREFIXES[Cons.REACTOME], clean_id)
        iri_source = Cons.SOURCE_NAMESPACES[Cons.REACTOME]
        source = Cons.REACTOME

    elif source == Cons.MINERVA_PATHWAY_COL:
        pathway_node = create_node(Cons.NODE_URI_PREFIXES[Cons.MINERVA], pathway_id)
        iri_source = Cons.SOURCE_NAMESPACES["Minerva"]

    if pathway_node:
        add_type(g, pathway_node, Cons.NODE_TYPES["pathway_node"])
        add_label(g, pathway_node, pathway_label)

        if iri_source:
            source_node = URIRef(iri_source)
            link_has_source(g, pathway_node, source_node)
            add_label(g, source_node, source)
            # Add both DCAT and VoID Dataset types (aligned with dataset_provenance.py)
            add_type(g, source_node, Cons.NODE_TYPES["data_source_node"])
            add_type(g, source_node, Cons.VOID_TYPES["dataset"])

    return pathway_node


def add_molecular_pathway_node(g, el, identifier, id_number):
    """Add a molecular pathway node and related triples to the RDF graph.

    :param g: RDF graph.
    :param el: Dictionary with pathway data.
    :param identifier: RDF node for gene or compound.
    :param id_number: Row identifier.
    """
    from pyBiodatafuse.graph.rdf.nodes.gene import get_gene_node

    pathway_id = safe_get(el, "pathway_id")
    pathway_label = safe_get(el, "pathway_label")
    target_gene = safe_get(el, "targetGene")
    target_protein = safe_get(el, "targetProtein")
    target_metabolite = safe_get(el, "targetMetabolite")
    mimtype = safe_get(el, "mimtype")
    rhea_id = safe_get(el, "rhea_id")

    if not mimtype:
        return

    # Create MIM type node
    mimtype_node = URIRef("http://vocabularies.wikipathways.org/wp#" + mimtype)
    add_type(g, mimtype_node, Cons.NODE_TYPES["interaction"])
    add_label(g, mimtype_node, mimtype)

    # Create interaction instance node
    mim_node = URIRef(g.new_uris["interaction"] + mimtype + "_" + id_number)
    add_type(g, mim_node, Cons.NODE_TYPES["interaction"])
    add_label(g, mim_node, mimtype + "_" + id_number)
    link_has_part(g, mim_node, identifier, bidirectional=True)

    # Link pathway
    if pathway_id and pathway_label:
        clean_pathway_id = extract_id(pathway_id)
        pathway_node = create_node("https://www.wikipathways.org/pathways/", clean_pathway_id)
        link_has_part(g, pathway_node, mimtype_node)
        link_has_part(g, pathway_node, identifier)
        link_has_part(g, pathway_node, mim_node)

    # Link target gene
    if target_gene:
        target_gene_node = get_gene_node(
            g,
            {"target.source": "Ensembl", "target": target_gene, "identifier": target_gene},
        )
        if target_gene_node:
            interaction_predicate = URIRef(g.new_uris["interaction"] + mimtype)
            wp_predicate = URIRef("http://vocabularies.wikipathways.org/wp#" + mimtype)
            add_triple(g, identifier, interaction_predicate, target_gene_node)
            add_triple(g, target_gene_node, wp_predicate, identifier)
            link_has_part(g, mim_node, target_gene_node, bidirectional=True)

    # Link target protein
    if target_protein:
        target_protein_node = create_node(Cons.BASE_URLS_DBS["uniprot"], target_protein)
        add_type(g, target_protein_node, Cons.NODE_TYPES["protein_node"])
        wp_predicate = URIRef("http://vocabularies.wikipathways.org/wp#" + mimtype)
        add_triple(g, identifier, wp_predicate, target_protein_node)
        add_triple(g, target_protein_node, wp_predicate, identifier)
        link_has_part(g, mim_node, target_protein_node, bidirectional=True)

    # Link target metabolite
    if target_metabolite:
        target_metabolite_node = URIRef(Cons.NODE_TYPES["compound_node"] + target_metabolite)
        add_type(g, target_metabolite_node, Cons.NODE_TYPES["compound_node"])
        wp_predicate = URIRef("http://vocabularies.wikipathways.org/wp#" + mimtype)
        add_triple(g, identifier, wp_predicate, target_metabolite_node)
        add_triple(g, target_metabolite_node, wp_predicate, identifier)
        link_has_part(g, mim_node, identifier, bidirectional=True)
        link_has_part(g, mim_node, target_metabolite_node, bidirectional=True)

    # Link rhea_id
    if rhea_id:
        add_same_as(g, mim_node, URIRef(rhea_id), bidirectional=False)
