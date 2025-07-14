# literature.py


"""Populate a BDF RDF graph with literature-based evidence."""


from rdflib import Literal, URIRef
from rdflib.namespace import OWL, RDF, RDFS

import pyBiodatafuse.constants as Cons
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

    :param g: (Graph): RDF graph to which the literature-based data is added.
    :param entry: (dict): Dictionary with literature-based association information.
    :param gene_node: (URIRef): URIRef of the gene node associated with the literature data.
    :param id_number: Unique identifier for the expression data.
    :param source_idx: Identifier for the source of the expression data.
    :param disease_data: Dictionary with disease data.
    :param new_uris: Node URIs for the graph.
    :param i: Row index.
    """
    source = entry.get("source", None)
    if source and "PMID" in source:
        source_id = source.split(":")[1].rstrip().lstrip()
        source_url = Cons.NODE_URI_PREFIXES["pubmed"] + source_id
        article_node = URIRef(source_url)
        umls = entry.get(Cons.UMLS, None)
        mondo = entry.get(Cons.MONDO, None)
        disease_name = entry[Cons.DISEASE_NAME]
        # Create the source node and add metadata
        g.add((article_node, RDF.type, URIRef(Cons.NODE_TYPES["article"])))
        g.add((article_node, URIRef(Cons.PREDICATES["sio_refers_to"]), gene_node))
        if umls:
            # disease_node = URIRef(f"{NAMESPACE_BINDINGS['umls']}{umls}")
            disease_node = add_gene_disease_associations(
                g, id_number, source_idx, gene_node, disease_data, new_uris, i
            )
            g.add(
                (
                    article_node,
                    URIRef(Cons.PREDICATES["sio_refers_to"]),
                    disease_node,
                )
            )
            g.add((disease_node, RDFS.label, Literal(disease_name)))
            g.add((disease_node, RDF.type, URIRef(Cons.NODE_TYPES["disease_node"])))
        if mondo:
            mondo_node = URIRef(f"{Cons.NAMESPACE_BINDINGS['mondo']}{mondo}")
            g.add(
                (
                    disease_node,
                    OWL.sameAs,
                    mondo_node,
                )
            )
            g.add((disease_node, RDFS.label, Literal(disease_name)))
            g.add((disease_node, RDF.type, URIRef(Cons.NODE_TYPES["disease_node"])))
