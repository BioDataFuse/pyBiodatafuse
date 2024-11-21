# protein_protein.py

"""Populate a BDF RDF graph with PPI nodes."""

from rdflib import Graph, Literal, URIRef
from rdflib.namespace import OWL, RDF, RDFS, XSD

from pyBiodatafuse.constants import BASE_URLS_DBS, NODE_TYPES, PREDICATES


def add_ppi_data(
    g: Graph, protein_node: URIRef, protein_name: str, entry: dict, base_uri: str, new_uris: dict
) -> URIRef:
    """Add a protein protein interaction node.

    :param g: RDFLib graph
    :param protein_node: URIRef for the target protein
    :param protein_name: the name of the target protein
    :param entry: the ppi dictionary
    :param base_uri: the base URI for the project
    :param new_uris: dictionary with project node URIs
    :return: a ppi URIRef node
    """
    stringdb_link_to = entry.get("stringdb_link_to", None)
    ensembl = entry.get("Ensembl", None).split(":")[1]
    score = entry.get("score", None)

    if score:
        score = float(score)
        # Nodes
        ppi_node = URIRef(base_uri + f"ppi/{protein_name}_{stringdb_link_to}")
        stringdb_node = URIRef(BASE_URLS_DBS["uniprot"] + f"{stringdb_link_to}")
        ensembl_node = URIRef(BASE_URLS_DBS["ensembl"] + f"{ensembl}")
        score_node = URIRef(f"{new_uris['score_base_node']}/{protein_name}_{ensembl}")

        g.add((ppi_node, RDF.type, URIRef(NODE_TYPES["ppi_node"])))
        g.add((stringdb_node, RDF.type, URIRef(NODE_TYPES["protein_node"])))
        g.add((ppi_node, URIRef(PREDICATES["sio_has_part"]), stringdb_node))
        g.add((ppi_node, URIRef(PREDICATES["sio_has_part"]), protein_node))
        g.add((stringdb_node, URIRef(PREDICATES["sio_is_part_of"]), ppi_node))
        g.add((stringdb_node, RDF.type, URIRef(NODE_TYPES["protein_node"])))
        g.add((stringdb_node, RDFS.label, Literal(stringdb_link_to, datatype=XSD.string)))
        g.add((ensembl_node, URIRef(PREDICATES["translates_to"]), stringdb_node))
        g.add((stringdb_node, URIRef(PREDICATES["translation_of"]), ensembl_node))
        g.add((stringdb_node, OWL.sameAs, ensembl_node))
        g.add((ensembl_node, RDF.type, URIRef(NODE_TYPES["gene_node"])))
        g.add((ppi_node, URIRef(PREDICATES["sio_has_measurement_value"]), score_node))
        g.add((score_node, RDF.type, URIRef(NODE_TYPES["score_node"])))
        g.add(
            (score_node, URIRef(PREDICATES["sio_has_value"]), Literal(score, datatype=XSD.double))
        )

        return ppi_node
