# dataset_provenance.py

"""
Dataset Provenance and Versioning for BioDataFuse RDF Graphs.

This module provides functions for adding dataset provenance information
to RDF graphs using:
- PROV-O (Provenance Ontology)
- DCAT (Data Catalog Vocabulary)
- Dublin Core Terms
- PAV (Provenance, Authoring and Versioning)
- VoID (Vocabulary of Interlinked Datasets)
- Schema.org

Each data source queried during graph construction is represented as a
dcat:Dataset with proper version information, access dates, and provenance.
"""

import logging
from datetime import datetime
from typing import Any, Dict, List, Optional, Set

from rdflib import Graph, Literal, URIRef
from rdflib.namespace import RDF, RDFS, XSD

import pyBiodatafuse.constants as Cons

# Set up logger
logger = logging.getLogger(__name__)

# Mapping from metadata datasource names to normalized constant names
# This ensures consistency between metadata-recorded sources and process-recorded sources
DATASOURCE_NAME_MAP = {
    # Open Targets variations
    "Open Targets GraphQL & REST API Beta": Cons.OPENTARGETS,
    "OpenTargets GraphQL & REST API Beta": Cons.OPENTARGETS,
    "Open Targets": Cons.OPENTARGETS,
    # StringDB variations
    "STRING": Cons.STRING,
    "String": Cons.STRING,
    "string": Cons.STRING,
    # DisGeNET variations
    "DisGeNET": Cons.DISGENET,
    "Disgenet": Cons.DISGENET,
    # AOP-Wiki variations
    "AOP-Wiki": Cons.AOPWIKIRDF,
    "AOPWiki": Cons.AOPWIKIRDF,
    "AOP Wiki": Cons.AOPWIKIRDF,
    # MINERVA variations
    "Minerva": Cons.MINERVA,
    # BridgeDB variations
    "BridgeDb": "BridgeDB",
    "bridgedb": "BridgeDB",
}


def normalize_datasource_name(name: str) -> str:
    """
    Normalize a datasource name to a consistent format.

    :param name: Raw datasource name from metadata or processing.
    :return: Normalized datasource name.
    """
    return DATASOURCE_NAME_MAP.get(name, name)


class DatasetProvenanceTracker:
    """Track and record dataset provenance for all data sources used in RDF generation."""

    def __init__(self, base_uri: str):
        """
        Initialize the provenance tracker.

        :param base_uri: Base URI for the RDF graph.
        """
        self.base_uri = base_uri
        self.used_datasets: Dict[str, Dict[str, Any]] = {}
        self.query_timestamps: Dict[str, datetime] = {}
        # Track nodes created per datasource for linking
        self.nodes_by_datasource: Dict[str, Set[URIRef]] = {}
        # Track interaction types added
        self.interaction_types: Dict[str, int] = {}

    def record_dataset_usage(
        self,
        datasource: str,
        version: Optional[str] = None,
        query_time: Optional[datetime] = None,
        endpoint_url: Optional[str] = None,
    ) -> None:
        """
        Record that a dataset was used during RDF generation.

        :param datasource: Name of the data source (e.g., "Bgee", "StringDB").
        :param version: Version string of the data source.
        :param query_time: Timestamp when the data source was queried.
        :param endpoint_url: URL of the endpoint that was queried.
        """
        # Normalize datasource name for consistency
        datasource = normalize_datasource_name(datasource)

        if datasource not in self.used_datasets:
            self.used_datasets[datasource] = {
                "version": version,
                "query_time": query_time or datetime.utcnow(),
                "endpoint_url": endpoint_url,
                "query_count": 0,
            }
            self.nodes_by_datasource[datasource] = set()
        self.used_datasets[datasource]["query_count"] += 1

        # Update version if provided and not already set
        if version and not self.used_datasets[datasource].get("version"):
            self.used_datasets[datasource]["version"] = version

    def record_node(self, datasource: str, node_uri: URIRef) -> None:
        """
        Record that a node was created from a specific data source.

        :param datasource: Name of the data source.
        :param node_uri: URIRef of the created node.
        """
        # Normalize datasource name for consistency
        datasource = normalize_datasource_name(datasource)

        if datasource not in self.nodes_by_datasource:
            self.nodes_by_datasource[datasource] = set()
        self.nodes_by_datasource[datasource].add(node_uri)

    def record_interaction_type(self, interaction_type: str) -> None:
        """
        Record that an interaction type was added to the graph.

        :param interaction_type: Type of interaction (e.g., "PPI", "gene-disease").
        """
        if interaction_type not in self.interaction_types:
            self.interaction_types[interaction_type] = 0
        self.interaction_types[interaction_type] += 1

    def get_used_datasets(self) -> Dict[str, Dict[str, Any]]:
        """
        Get all datasets that were used.

        :return: Dictionary of used datasets with their metadata.
        """
        return self.used_datasets

    def get_nodes_for_datasource(self, datasource: str) -> Set[URIRef]:
        """
        Get all nodes created from a specific data source.

        :param datasource: Name of the data source.
        :return: Set of URIRefs for nodes from this data source.
        """
        return self.nodes_by_datasource.get(datasource, set())

    def log_summary(self) -> None:
        """Log a summary of all datasets and interactions added to the graph."""
        logger.info("=" * 60)
        logger.info("RDF Graph Generation Summary")
        logger.info("=" * 60)

        logger.info("\nDatasets used:")
        for datasource, metadata in self.used_datasets.items():
            node_count = len(self.nodes_by_datasource.get(datasource, set()))
            version = metadata.get("version", "unknown")
            query_count = metadata.get("query_count", 0)
            logger.info(
                f"  - {datasource}: {node_count} nodes, version={version}, queries={query_count}"
            )

        if self.interaction_types:
            logger.info("\nInteraction types added:")
            for interaction_type, count in sorted(
                self.interaction_types.items(), key=lambda x: -x[1]
            ):
                logger.info(f"  - {interaction_type}: {count}")

        total_nodes = sum(len(nodes) for nodes in self.nodes_by_datasource.values())
        logger.info(f"\nTotal unique nodes linked to datasets: {total_nodes}")
        logger.info("=" * 60)


def create_dataset_node(
    g: Graph,
    datasource: str,
    base_uri: str,
    version: Optional[str] = None,
    query_time: Optional[datetime] = None,
    endpoint_url: Optional[str] = None,
) -> URIRef:
    """
    Create a DCAT Dataset node for a data source with full provenance.

    Uses the canonical source URI from DATA_SOURCE_IDENTIFIERS as the dataset URI
    when available, falling back to a constructed URI otherwise.

    :param g: RDF graph to add the dataset node to.
    :param datasource: Name of the data source.
    :param base_uri: Base URI for the graph.
    :param version: Version string of the data source.
    :param query_time: Timestamp when the data was retrieved.
    :param endpoint_url: URL of the endpoint that was queried.
    :return: URIRef of the created dataset node.
    """
    # Use canonical source URI from DATA_SOURCE_IDENTIFIERS if available
    if datasource in Cons.DATA_SOURCE_IDENTIFIERS:
        dataset_uri = URIRef(Cons.DATA_SOURCE_IDENTIFIERS[datasource])
    else:
        # Fallback to constructed URI for unknown sources
        dataset_uri = URIRef(f"{base_uri}dataset/{datasource.lower().replace(' ', '_')}")

    # Add dataset type (DCAT Dataset)
    g.add((dataset_uri, RDF.type, URIRef(Cons.DCAT_TYPES["dataset"])))
    g.add((dataset_uri, RDF.type, URIRef(Cons.VOID_TYPES["dataset"])))

    # Add label
    g.add((dataset_uri, RDFS.label, Literal(datasource, datatype=XSD.string)))

    # Add title using Dublin Core
    g.add(
        (
            dataset_uri,
            URIRef(Cons.DCTERMS_PREDICATES["title"]),
            Literal(f"{datasource} Dataset", datatype=XSD.string),
        )
    )

    # Add version information
    if version:
        g.add(
            (
                dataset_uri,
                URIRef(Cons.PAV_PREDICATES["version"]),
                Literal(version, datatype=XSD.string),
            )
        )
        g.add(
            (
                dataset_uri,
                URIRef(Cons.DCTERMS_PREDICATES["version"]),
                Literal(version, datatype=XSD.string),
            )
        )

    # Add retrieval timestamp
    if query_time:
        timestamp_literal = Literal(query_time.isoformat() + "Z", datatype=XSD.dateTime)
        g.add((dataset_uri, URIRef(Cons.PAV_PREDICATES["retrieved_on"]), timestamp_literal))
        g.add((dataset_uri, URIRef(Cons.PROV_PREDICATES["generated_at_time"]), timestamp_literal))

    # Add source identifier
    if datasource in Cons.DATA_SOURCE_IDENTIFIERS:
        source_uri = URIRef(Cons.DATA_SOURCE_IDENTIFIERS[datasource])
        g.add((dataset_uri, URIRef(Cons.PAV_PREDICATES["retrieved_from"]), source_uri))
        g.add((dataset_uri, URIRef(Cons.DCTERMS_PREDICATES["source"]), source_uri))

    # Add endpoint URL
    if endpoint_url:
        endpoint_uri = URIRef(endpoint_url)
        g.add((dataset_uri, URIRef(Cons.DCAT_PREDICATES["access_url"]), endpoint_uri))

        # Create a DataService node for the endpoint
        service_uri = URIRef(f"{base_uri}service/{datasource.lower().replace(' ', '_')}")
        g.add((service_uri, RDF.type, URIRef(Cons.DCAT_TYPES["data_service"])))
        g.add((service_uri, RDF.type, URIRef(Cons.SCHEMA_TYPES["web_api"])))
        g.add((service_uri, URIRef(Cons.DCAT_PREDICATES["endpoint_url"]), endpoint_uri))
        g.add((service_uri, URIRef(Cons.DCAT_PREDICATES["serves_dataset"]), dataset_uri))
        g.add((service_uri, RDFS.label, Literal(f"{datasource} API", datatype=XSD.string)))

    # Add SPARQL endpoint if applicable
    if datasource in Cons.DATA_SOURCE_SPARQL_ENDPOINTS:
        sparql_uri = URIRef(Cons.DATA_SOURCE_SPARQL_ENDPOINTS[datasource])
        g.add((dataset_uri, URIRef(Cons.VOID_PREDICATES["sparql_endpoint"]), sparql_uri))

    # Add license information
    if datasource in Cons.DATA_SOURCE_LICENSES:
        license_uri = URIRef(Cons.DATA_SOURCE_LICENSES[datasource])
        g.add((dataset_uri, URIRef(Cons.DCTERMS_PREDICATES["license"]), license_uri))

    # Add landing page from DATA_SOURCES
    if datasource in Cons.DATA_SOURCES:
        landing_page_uri = URIRef(Cons.DATA_SOURCES[datasource])
        g.add((dataset_uri, URIRef(Cons.DCAT_PREDICATES["landing_page"]), landing_page_uri))
        g.add((dataset_uri, URIRef(Cons.SCHEMA_PREDICATES["url"]), landing_page_uri))

    return dataset_uri


def add_dataset_provenance_to_graph(
    g: Graph,
    base_uri: str,
    graph_uri: str,
    tracker: DatasetProvenanceTracker,
) -> None:
    """
    Add all tracked dataset provenance to the RDF graph.

    This function:
    1. Creates a DCAT Catalog node for the BioDataFuse graph
    2. Creates DCAT Dataset nodes for each data source used
    3. Links all nodes created from each data source to their source dataset

    :param g: RDF graph to add provenance to.
    :param base_uri: Base URI for the graph.
    :param graph_uri: URI of the main graph/catalog.
    :param tracker: DatasetProvenanceTracker with recorded dataset usage.
    """
    graph_resource = URIRef(graph_uri)

    # Create a DCAT Catalog node for the BioDataFuse graph
    g.add((graph_resource, RDF.type, URIRef(Cons.DCAT_TYPES["catalog"])))
    g.add(
        (
            graph_resource,
            URIRef(Cons.DCTERMS_PREDICATES["title"]),
            Literal("BioDataFuse Knowledge Graph", datatype=XSD.string),
        )
    )

    # Store dataset URIs for linking nodes
    dataset_uris: Dict[str, URIRef] = {}

    # Add each used dataset
    for datasource, metadata in tracker.get_used_datasets().items():
        # Determine endpoint URL
        endpoint_url = metadata.get("endpoint_url")
        if not endpoint_url:
            if datasource in Cons.DATA_SOURCE_SPARQL_ENDPOINTS:
                endpoint_url = Cons.DATA_SOURCE_SPARQL_ENDPOINTS[datasource]
            elif datasource in Cons.DATA_SOURCE_API_ENDPOINTS:
                endpoint_url = Cons.DATA_SOURCE_API_ENDPOINTS[datasource]

        # Create dataset node
        dataset_uri = create_dataset_node(
            g=g,
            datasource=datasource,
            base_uri=base_uri,
            version=metadata.get("version"),
            query_time=metadata.get("query_time"),
            endpoint_url=endpoint_url,
        )
        dataset_uris[datasource] = dataset_uri

        # Link graph to dataset using PROV-O
        g.add((graph_resource, URIRef(Cons.PROV_PREDICATES["was_derived_from"]), dataset_uri))

        # Also use DCAT distribution relationship
        g.add((graph_resource, URIRef(Cons.DCTERMS_PREDICATES["source"]), dataset_uri))

    # Link all tracked nodes to their source datasets
    source_predicate = URIRef(Cons.SOURCE_DATASET_PREDICATE)
    nodes_linked = 0
    for datasource, nodes in tracker.nodes_by_datasource.items():
        if datasource in dataset_uris:
            dataset_uri = dataset_uris[datasource]
            for node in nodes:
                g.add((node, source_predicate, dataset_uri))
                nodes_linked += 1

    logger.info(f"Linked {nodes_linked} data nodes to their source datasets")

    # Log summary
    tracker.log_summary()


def record_datasource_from_metadata(
    tracker: DatasetProvenanceTracker,
    metadata_list: list,
) -> None:
    """
    Record dataset usage from BioDataFuse metadata list.

    :param tracker: DatasetProvenanceTracker instance.
    :param metadata_list: List of metadata dictionaries from BioDataFuse queries.
    """
    for entry in metadata_list:
        if not isinstance(entry, dict):
            continue

        datasource = entry.get("datasource")
        if not datasource:
            continue

        # Extract version
        version = None
        metadata_inner = entry.get("metadata", {})
        if isinstance(metadata_inner, dict):
            version = metadata_inner.get("source_version")

        # Extract query info
        query_info = entry.get("query", {})
        endpoint_url = query_info.get("url") if isinstance(query_info, dict) else None

        # Parse query date
        query_time = None
        if isinstance(query_info, dict) and "date" in query_info:
            try:
                query_time = datetime.strptime(query_info["date"], "%Y-%m-%d %H:%M:%S")
            except (ValueError, TypeError):
                query_time = datetime.utcnow()

        tracker.record_dataset_usage(
            datasource=datasource,
            version=version,
            query_time=query_time,
            endpoint_url=endpoint_url,
        )


def get_datasource_for_column(column_name: str) -> Optional[str]:
    """
    Get the data source name for a given DataFrame column.

    :param column_name: Name of the column in the BioDataFuse DataFrame.
    :return: Data source name or None if not found.
    """
    column_to_datasource = {
        Cons.BGEE_GENE_EXPRESSION_LEVELS_COL: Cons.BGEE,
        Cons.COMPOUNDWIKI_COL: Cons.COMPOUNDWIKI,
        Cons.DISGENET_DISEASE_COL: Cons.DISGENET,
        Cons.INTACT_INTERACT_COL: Cons.INTACT,
        Cons.INTACT_COMPOUND_INTERACT_COL: Cons.INTACT,
        Cons.KEGG_PATHWAY_COL: Cons.KEGG,
        Cons.KEGG_COMPOUND_COL: Cons.KEGG,
        Cons.MINERVA_PATHWAY_COL: Cons.MINERVA,
        Cons.MOLMEDB_PROTEIN_COMPOUND_COL: Cons.MOLMEDB,
        Cons.MOLMEDB_COMPOUND_PROTEIN_COL: Cons.MOLMEDB,
        Cons.OPENTARGETS_DISEASE_COL: Cons.OPENTARGETS,
        Cons.OPENTARGETS_GO_COL: Cons.OPENTARGETS,
        Cons.OPENTARGETS_GENE_COMPOUND_COL: Cons.OPENTARGETS,
        Cons.OPENTARGETS_REACTOME_COL: Cons.OPENTARGETS,
        Cons.PUBCHEM_COMPOUND_ASSAYS_COL: Cons.PUBCHEM,
        Cons.STRING_INTERACT_COL: Cons.STRING,
        Cons.WIKIDATA_CC_COL: Cons.WIKIDATA,
        Cons.WIKIPATHWAYS_MOLECULAR_COL: Cons.WIKIPATHWAYS,
        Cons.WIKIPATHWAYS_PATHWAY_COL: Cons.WIKIPATHWAYS,
        Cons.AOPWIKI_GENE_COL: Cons.AOPWIKIRDF,
        Cons.AOPWIKI_COMPOUND_COL: Cons.AOPWIKIRDF,
    }

    return column_to_datasource.get(column_name)


def get_active_datasources_from_row(row) -> Set[str]:
    """
    Determine which data sources have data in a given row.

    :param row: A pandas Series representing a row of data.
    :return: Set of data source names that have non-null data in the row.
    """
    active_sources = set()

    # Check each known column
    columns_to_check = [
        Cons.BGEE_GENE_EXPRESSION_LEVELS_COL,
        Cons.COMPOUNDWIKI_COL,
        Cons.DISGENET_DISEASE_COL,
        Cons.INTACT_INTERACT_COL,
        Cons.INTACT_COMPOUND_INTERACT_COL,
        Cons.KEGG_PATHWAY_COL,
        Cons.MOLMEDB_PROTEIN_COMPOUND_COL,
        Cons.MOLMEDB_COMPOUND_PROTEIN_COL,
        Cons.OPENTARGETS_DISEASE_COL,
        Cons.OPENTARGETS_GO_COL,
        Cons.OPENTARGETS_GENE_COMPOUND_COL,
        Cons.PUBCHEM_COMPOUND_ASSAYS_COL,
        Cons.STRING_INTERACT_COL,
        Cons.WIKIPATHWAYS_MOLECULAR_COL,
        Cons.AOPWIKI_GENE_COL,
        Cons.AOPWIKI_COMPOUND_COL,
    ]

    for col in columns_to_check:
        value = row.get(col)
        if value is not None and (not isinstance(value, list) or len(value) > 0):
            datasource = get_datasource_for_column(col)
            if datasource:
                active_sources.add(datasource)

    return active_sources
