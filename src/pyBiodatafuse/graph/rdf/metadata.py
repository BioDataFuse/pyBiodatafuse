# metadata.py


"""Some metadata node functions."""


from datetime import datetime

from rdflib import Graph, Literal, URIRef
from rdflib.namespace import DCTERMS, OWL, RDF, RDFS, XSD

from pyBiodatafuse.constants import DATA_SOURCES, NODE_TYPES


def add_metadata(
    g: Graph, graph_uri: str, metadata: dict, version_iri=None, author=None, orcid=None
):
    """
    Add metadata to the RDF graph, including creation date, version, author, and ORCID.

    :param g: The RDF graph to which metadata is added.
    :param graph_uri: URI identifying the RDF graph.
    :param metadata: Combined metadata for a BioDatafuse query.
    :param version_iri: Version IRI to add.
    :param author: Author's name.
    :param orcid: Author's ORCID.
    """
    # Automatically get the current date and time in ISO 8601 format
    creation_date = datetime.utcnow().isoformat() + "Z"

    # Define the URI for the graph resource
    graph_resource = URIRef(graph_uri)

    # Add the creation date to the graph resource
    g.add((graph_resource, DCTERMS.created, Literal(creation_date, datatype=XSD.dateTime)))

    if version_iri:
        g.add((graph_resource, DCTERMS.identifier, URIRef(version_iri)))

    if orcid:
        orcid_node = URIRef(orcid)
        g.add((orcid_node, OWL.sameAs, Literal(author, datatype=XSD.string)))
        g.add((graph_resource, DCTERMS.creator, orcid_node))
    # Metadata file: add version, api version, and dates to sources
    api_version = None
    version = None
    for source in [
        "DISGENET",
        "Open Targets GraphQL & REST API Beta",
        "MINERVA",
        "WikiPathways",
        "StringDB",
        "BridgeDb",
    ]:
        entries = [
            item for item in metadata if isinstance(item, dict) and item.get("datasource") == source
        ]
        for entry in entries:
            data_version = None
            year = None
            month = None
            version = None
            api_version = None
            # date = None
            url_service = None
            source_node = None
            # input_type = entry.get("query").get("input_type")
            # date = entry.get("query").get("date")
            query = entry.get("query", None)
            if query:
                url_service = query.get("url", None)
            source = entry.get("source", None)  # type: ignore
            if source == "Open Targets GraphQL & REST API Beta":
                source_node = URIRef(DATA_SOURCES["OpenTargets"])
                g.add(
                    (
                        source_node,
                        RDFS.label,
                        Literal("The Open Targets Platform", datatype=XSD.string),
                    )
                )
                metadata_f = entry.get("metadata", None)
                if metadata_f:
                    data_version = metadata_f.get("data_version", None)
                    if data_version:
                        year = data_version.get("year")
                        month = data_version.get("month")
                        version = f"{year}-{month}-01"
                source_version = metadata_f.get("source_version") if metadata_f else None
                api_version_data = source_version.get("apiVersion") if source_version else None
                if api_version_data:
                    api_version = ".".join([str(i) for i in api_version_data.values()])
                else:
                    api_version = None
            elif source == "DISGENET":
                source_node = URIRef(DATA_SOURCES[source])
                g.add(
                    (
                        source_node,
                        RDFS.label,
                        Literal("DisGeNET", datatype=XSD.string),
                    )
                )
                if metadata_f:
                    data_version = metadata_f.get("version", None)
                    version = metadata_f.get("version", None)

            elif source == "MINERVA":
                source_node = URIRef(DATA_SOURCES[source])
                # project = entry.get('query').get('MINERVA project')
                metadata_f = entry.get("metadata", None)
                if metadata_f:
                    version = metadata_f.get("source_version")
                g.add(
                    (
                        source_node,
                        RDFS.label,
                        Literal("The MINERVA Platform", datatype=XSD.string),
                    )
                )

            elif source == "WikiPathways":
                metadata_f = entry.get("metadata", None)
                if metadata_f:
                    version = metadata_f.get("source_version", None)
                source_node = URIRef(DATA_SOURCES[source])
                g.add(
                    (
                        source_node,
                        RDFS.label,
                        Literal("WikiPathways", datatype=XSD.string),
                    )
                )

            elif source == "BridgeDb":  # Several metadata fields left unRDFied
                source_node = URIRef(DATA_SOURCES[source])
                g.add(
                    (
                        source_node,
                        RDFS.label,
                        Literal("BridgeDb", datatype=XSD.string),
                    )
                )
                metadata_f = entry.get("metadata", None)
                if metadata_f:
                    version = metadata_f.get("bridgedb.version", None)
                    api_version = metadata_f.get("webservice.version", None)

            elif source == "StringDB":
                source_node = URIRef(DATA_SOURCES[source])
                g.add(
                    (
                        source_node,
                        RDFS.label,
                        Literal("STRING", datatype=XSD.string),
                    )
                )

            # Add api node for source
            api_node = None
            if url_service:
                api_node = URIRef(url_service)
                g.add((api_node, RDF.type, URIRef("https://schema.org/WebAPI")))
            if source_node:
                if api_node:
                    g.add((api_node, URIRef("https://schema.org/provider"), source_node))
                g.add((source_node, RDF.type, URIRef(NODE_TYPES["data_source_node"])))
                if api_version and api_node:
                    g.add(
                        (
                            api_node,
                            URIRef("http://www.geneontology.org/formats/oboInOwl#version"),
                            Literal(api_version),
                        )
                    )
            if version and source_node:
                g.add(
                    (
                        source_node,
                        URIRef("http://www.geneontology.org/formats/oboInOwl#version"),
                        Literal(version),
                    )
                )


def add_data_source_node(g: Graph, source: str) -> URIRef:
    """Create and add a data source node to the RDF graph.

    :param g: RDF graph to which the data source node will be added.
    :param source: String containing the name of the source of the data
    :return: URIRef for the created data source node.
    """
    data_source_name = Literal(source, datatype=XSD.string)
    data_source_url = URIRef(DATA_SOURCES[source])
    g.add((data_source_url, RDF.type, URIRef(NODE_TYPES["data_source_node"])))
    g.add((data_source_url, RDFS.label, data_source_name))
    return data_source_url
