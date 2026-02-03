# literature.py

"""Populate a BDF RDF graph with literature-based evidence."""

import pyBiodatafuse.constants as Cons
from pyBiodatafuse.graph.rdf.nodes.base import (
    add_label,
    add_same_as,
    add_type,
    create_node,
    link_refers_to,
    safe_get,
)
from pyBiodatafuse.graph.rdf.nodes.gene_disease import add_gene_disease_associations


def add_literature_based_data(
    g,
    entry,
    gene_node,
    id_number,
    source_idx,
    disease_data,
    new_uris,
    i,
) -> None:
    """Add literature-based data for gene associations.

    :param g: RDF graph to which the literature-based data is added.
    :param entry: Dictionary with literature-based association information.
    :param gene_node: URIRef of the gene node associated with the literature data.
    :param id_number: Unique identifier for the expression data.
    :param source_idx: Identifier for the source of the expression data.
    :param disease_data: Dictionary with disease data.
    :param new_uris: Node URIs for the graph.
    :param i: Row index.
    """
    source = safe_get(entry, "source")
    if not source or "PMID" not in source:
        return

    source_id = source.split(":")[1].strip()
    article_node = create_node(Cons.NODE_URI_PREFIXES["pubmed"], source_id)

    umls = safe_get(entry, Cons.UMLS)
    mondo = safe_get(entry, Cons.MONDO)
    disease_name = safe_get(entry, Cons.DISEASE_NAME)

    # Create the article node and add metadata
    add_type(g, article_node, Cons.NODE_TYPES["article"])
    link_refers_to(g, article_node, gene_node)

    if umls:
        # Create disease node via gene-disease association
        disease_node = add_gene_disease_associations(
            g, id_number, source_idx, gene_node, disease_data, new_uris, i
        )
        if disease_node:
            link_refers_to(g, article_node, disease_node)
            add_label(g, disease_node, disease_name)
            add_type(g, disease_node, Cons.NODE_TYPES["disease_node"])

            if mondo:
                mondo_node = create_node(Cons.NAMESPACE_BINDINGS["mondo"], mondo)
                add_same_as(g, disease_node, mondo_node, bidirectional=False)
                add_type(g, mondo_node, Cons.NODE_TYPES["disease_node"])
