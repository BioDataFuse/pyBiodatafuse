# metadata.py


"""Some metadata node functions."""

from datetime import datetime
from typing import Dict, List, Optional

from rdflib import Graph, Literal, Namespace, URIRef
from rdflib.namespace import DCTERMS, RDF, RDFS, XSD

from pyBiodatafuse.constants import (
    DATA_SOURCES,
    FOAF_NAMESPACE,
    FOAF_TYPES,
    NODE_TYPES,
    SCHEMA_NAMESPACE,
    SCHEMA_TYPES,
)

# Define namespace objects
FOAF = Namespace(FOAF_NAMESPACE)
SCHEMA = Namespace(SCHEMA_NAMESPACE)


def add_creator_node(
    g: Graph,
    graph_resource: URIRef,
    name: str,
    orcid: Optional[str] = None,
    url: Optional[str] = None,
) -> URIRef:
    """
    Create and add a properly modeled creator node to the graph.

    Models the creator as both foaf:Person and schema:Person with foaf:name.
    If an ORCID is provided, uses it as the node URI. Otherwise uses a provided URL
    or generates a blank node.

    :param g: The RDF graph to which the creator node is added.
    :param graph_resource: The graph resource URI to link the creator to.
    :param name: The name of the creator.
    :param orcid: The ORCID URL of the creator (e.g., 'https://orcid.org/0000-0001-2345-6789').
    :param url: An alternative URL to identify the creator (if no ORCID provided).
    :return: The URIRef of the created person node.
    """
    # Determine the URI for the person node
    if orcid:
        person_node = URIRef(orcid)
    elif url:
        person_node = URIRef(url)
    else:
        # Generate a URI based on the name (slugified)
        slug = name.lower().replace(" ", "_").replace(".", "")
        person_node = URIRef(f"{graph_resource}/creator/{slug}")

    # Add types: foaf:Person and schema:Person
    g.add((person_node, RDF.type, URIRef(FOAF_TYPES["person"])))
    g.add((person_node, RDF.type, URIRef(SCHEMA_TYPES["person"])))

    # Add the name using foaf:name
    g.add((person_node, FOAF.name, Literal(name, datatype=XSD.string)))

    # Link the creator to the graph resource
    g.add((graph_resource, DCTERMS.creator, person_node))

    return person_node


def add_metadata(
    g: Graph,
    graph_uri: str,
    metadata: dict,
    version_iri: Optional[str] = None,
    title: Optional[str] = None,
    description: Optional[str] = None,
    author: Optional[str] = None,
    orcid: Optional[str] = None,
    creators: Optional[List[Dict[str, str]]] = None,
):
    """
    Add metadata to the RDF graph, including creation date, version, title, and creators.

    :param g: The RDF graph to which metadata is added.
    :param graph_uri: URI identifying the RDF graph.
    :param metadata: Combined metadata for a BioDatafuse query.
    :param version_iri: Version IRI to add (optional).
    :param title: Title of the graph (optional).
    :param description: Description of the graph (optional).
    :param author: Author's name (optional, deprecated - use creators instead).
    :param orcid: Author's ORCID (optional, deprecated - use creators instead).
    :param creators: List of creator dictionaries with 'name' (required),
                    'orcid' (optional), and 'url' (optional) keys.
    """
    # Automatically get the current date and time in ISO 8601 format
    creation_date = datetime.utcnow().isoformat() + "Z"

    # Define the URI for the graph resource
    graph_resource = URIRef(graph_uri)

    # Add the creation date to the graph resource
    g.add((graph_resource, DCTERMS.created, Literal(creation_date, datatype=XSD.dateTime)))

    # Add title if provided
    if title:
        g.add((graph_resource, DCTERMS.title, Literal(title, datatype=XSD.string)))

    # Add description if provided
    if description:
        g.add((graph_resource, DCTERMS.description, Literal(description, datatype=XSD.string)))

    if version_iri:
        g.add((graph_resource, DCTERMS.identifier, URIRef(version_iri)))

    # Handle creators list (preferred method)
    if creators:
        for creator in creators:
            name = creator.get("name")
            if name:
                add_creator_node(
                    g,
                    graph_resource,
                    name=name,
                    orcid=creator.get("orcid"),
                    url=creator.get("url"),
                )

    # Handle legacy author/orcid parameters (for backwards compatibility)
    elif author and orcid:
        add_creator_node(g, graph_resource, name=author, orcid=orcid)
    elif author:
        add_creator_node(g, graph_resource, name=author)

    # Metadata file: add version, api version, and dates to sources
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
