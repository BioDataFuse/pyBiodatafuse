# coding: utf-8

"""
RDF Graph Generation for BioDataFuse.

This module defines the BDFGraph class for generating RDF knowledge graphs
from BioDataFuse annotated data. The class extends rdflib.Graph with methods
for processing biological data and generating standardized RDF output.

Architecture
------------
The module uses a layered architecture:

1. BDFGraph: Main class that orchestrates RDF generation
2. Process methods: Handle specific data types (PPI, disease, expression, etc.)
3. Node modules: Create individual RDF nodes (in nodes/ subdirectory)
4. Base utilities: Common node creation functions (in nodes/base.py)

Data Flow
---------
DataFrame row -> process_row() -> _process_gene_data() or _process_compound_data()
    -> individual process_* methods -> node module functions -> RDF triples

Error Handling
--------------
All process_* method calls are wrapped with _safe_process() which:
- Catches exceptions and logs warnings
- Continues processing remaining data
- Returns success/failure status

Extending
---------
To add a new data type:

1. Create a node module in nodes/ (see nodes/__init__.py for pattern)
2. Add a process_* method to BDFGraph
3. Call process_* from _process_gene_data or _process_compound_data
4. Add the data source to constants.py

Example::

    def process_new_data(self, data, node):
        if not data:
            return
        for entry in data:
            result = add_new_node(self, entry, node)
            if result:
                self.record_datasource(Cons.NEW_SOURCE, node=node, interaction_type="new type")
"""

import time
from typing import Any, Dict, List, Optional, Union

import pandas as pd
from bioregistry import normalize_curie
from rdflib import Graph, URIRef
from tqdm import tqdm

import pyBiodatafuse.constants as Cons
from pyBiodatafuse.graph.rdf.metadata import add_metadata
from pyBiodatafuse.graph.rdf.nodes.aop import add_aop_data
from pyBiodatafuse.graph.rdf.nodes.compound import (
    add_associated_compound_node,
    add_compoundwiki_annotations,
    add_inhibitor_transporter_node,
    add_transporter_inhibitor_node,
    get_compound_node,
)
from pyBiodatafuse.graph.rdf.nodes.dataset_provenance import (
    DatasetProvenanceTracker,
    add_dataset_provenance_to_graph,
    record_datasource_from_metadata,
)
from pyBiodatafuse.graph.rdf.nodes.experimental_process import add_pubchem_assay_node
from pyBiodatafuse.graph.rdf.nodes.gene import get_gene_node
from pyBiodatafuse.graph.rdf.nodes.gene_disease import add_gene_disease_associations
from pyBiodatafuse.graph.rdf.nodes.gene_expression import add_gene_expression_data
from pyBiodatafuse.graph.rdf.nodes.go_terms import add_go_cpf
from pyBiodatafuse.graph.rdf.nodes.literature import add_literature_based_data
from pyBiodatafuse.graph.rdf.nodes.pathway import add_molecular_pathway_node, add_pathway_node
from pyBiodatafuse.graph.rdf.nodes.protein_protein import add_ppi_data
from pyBiodatafuse.graph.rdf.utils import (
    discover_prefixes_from_graph,
    get_shacl_prefixes,
    get_shapes,
    replace_na_none,
)
from pyBiodatafuse.id_mapper import read_datasource_file
from pyBiodatafuse.logging_config import get_logger

# Set up logger
logger = get_logger(__name__)


class BDFGraph(Graph):
    """Main class for a BioDatafuse RDF Graph, superclass of rdflib.Graph."""

    def __init__(
        self,
        base_uri: str,
        version_iri: Optional[str] = None,
        title: Optional[str] = None,
        description: Optional[str] = None,
        author: Optional[str] = None,
        orcid: Optional[str] = None,
        creators: Optional[List[Dict[str, str]]] = None,
    ):
        """
        Initialize a new instance of the class with the provided metadata and URIs.

        :param base_uri: The base URI for the RDF graph.
        :param version_iri: The version IRI for the RDF graph (optional).
        :param title: The title of the BDF graph (optional).
        :param description: A description of the BDF graph (optional).
        :param author: The author of the BDF graph (optional, use creators for multiple).
        :param orcid: The ORCID identifier for the author (optional, use creators for multiple).
        :param creators: A list of creator dictionaries, each with 'name' (required),
                        'orcid' (optional), and 'url' (optional) keys. Example:
                        [{'name': 'John Doe', 'orcid': 'https://orcid.org/0000-0001-2345-6789'},
                         {'name': 'Jane Smith', 'url': 'https://example.org/jane'}]
        """
        # Initialize the rdflib.Graph superclass without passing extra arguments
        super().__init__()

        # Assign parameters to instance attributes
        self.base_uri = base_uri
        self.version_iri = version_iri
        self.title = title
        self.description = description
        self.author = author
        self.orcid = orcid
        self.creators = creators or []

        # Create and bind custom URIs and namespaces
        self.new_uris = {key: self.base_uri + value for key, value in Cons.URIS.items()}
        self._shex_path = None
        self._shacl_path = None
        self._prefixes_path = None
        self._namespaces = {}
        self.include_variants = False  # TODO: Allow user to set options that can affect graph size

        # Initialize dataset provenance tracker
        self.provenance_tracker = DatasetProvenanceTracker(base_uri)

        # Bind prefixes
        for key, new_value in self.new_uris.items():
            self.bind(key, new_value)
        for key, value in Cons.NAMESPACE_BINDINGS.items():
            self.bind(key, value)

        # Bind provenance and dataset vocabulary namespaces
        self.bind("prov", Cons.PROV_NAMESPACE)
        self.bind("dcat", Cons.DCAT_NAMESPACE)
        self.bind("dcterms", Cons.DCTERMS_NAMESPACE)
        self.bind("pav", Cons.PAV_NAMESPACE)
        self.bind("void", Cons.VOID_NAMESPACE)
        self.bind("schema", Cons.SCHEMA_NAMESPACE)
        self.bind("foaf", Cons.FOAF_NAMESPACE)

    def generate_rdf(
        self, df: pd.DataFrame, metadata: Dict[str, Any], open_only: bool = False
    ) -> None:
        """
        Generate an RDF graph from the provided DataFrame and metadata.

        :param df: The DataFrame containing the data to be converted into RDF.
        :param metadata: A dictionary containing metadata information for RDF generation.
        :param open_only: A flag indicating whether to process only open data. Defaults to False.
        :param metadata: Metadata information to be added to the RDF graph.
        """
        start_time = time.time()
        logger.info("Starting RDF graph generation...")

        # Record dataset usage from metadata for provenance tracking
        if isinstance(metadata, list):
            record_datasource_from_metadata(self.provenance_tracker, metadata)

        df = df.applymap(replace_na_none)
        datasources = read_datasource_file()
        if not self.include_variants:
            # Filter to keep only rows with valid gene identifier sources from datasources file
            valid_gene_sources = list(datasources[datasources["type"] == "gene"]["source"])
            df = df[df[Cons.TARGET_SOURCE_COL].isin(valid_gene_sources)]

        total_rows = df.shape[0]
        logger.info(f"Processing {total_rows} rows...")

        for j, (_, row) in enumerate(
            tqdm(df.iterrows(), total=total_rows, desc="Building RDF graph")
        ):
            self.process_row(row, j, datasources)

        self._add_metadata(metadata)

        # Add dataset provenance nodes to the graph
        logger.info("Adding dataset provenance nodes...")
        graph_uri = self.version_iri if self.version_iri is not None else self.base_uri
        add_dataset_provenance_to_graph(
            g=self,
            base_uri=self.base_uri,
            graph_uri=graph_uri,
            tracker=self.provenance_tracker,
        )

        elapsed = time.time() - start_time
        logger.info(f"RDF graph generation completed in {elapsed:.2f} seconds")
        logger.info(f"Total triples in graph: {len(self)}")

        # Discover additional prefixes from URIs in the graph using bioregistry
        logger.info("Discovering prefixes from graph URIs using bioregistry...")
        discovered = discover_prefixes_from_graph(self)
        if discovered:
            logger.info(f"Discovered {len(discovered)} additional prefixes")
            # Store discovered prefixes for use in shacl_prefixes()
            self._namespaces.update(discovered)
            # Also bind them to the graph for serialization
            for prefix, ns_uri in discovered.items():
                self.bind(prefix, ns_uri)

    def record_datasource(
        self, datasource: str, node: Optional[URIRef] = None, interaction_type: Optional[str] = None
    ) -> None:
        """
        Record that data from a specific data source was added to the graph.

        This method is called by process_* methods when they successfully
        add data from a data source to the RDF graph.

        :param datasource: Name of the data source (e.g., "StringDB", "Bgee").
        :param node: Optional URIRef of a node created from this data source.
        :param interaction_type: Optional type of interaction being added.
        """
        self.provenance_tracker.record_dataset_usage(datasource=datasource)
        if node is not None:
            self.provenance_tracker.record_node(datasource, node)
        if interaction_type is not None:
            self.provenance_tracker.record_interaction_type(interaction_type)

    def _safe_process(self, func, *args, error_msg: str = "Processing failed", **kwargs) -> bool:
        """
        Safely execute a processing function with error handling.

        :param func: Function to execute.
        :param args: Arguments to pass to the function.
        :param error_msg: Error message prefix for logging.
        :param kwargs: Keyword arguments to pass to the function.
        :return: True if successful, False otherwise.
        """
        try:
            func(*args, **kwargs)
            return True
        except Exception as e:
            logger.warning("%s: %s", error_msg, e)
            return False

    def _safe_iterate(
        self, items: list, func, *args, error_msg: str = "Item processing failed", **kwargs
    ) -> int:
        """
        Safely iterate over items and process each one, continuing on errors.

        :param items: List of items to iterate over.
        :param func: Function to call for each item. First argument will be the item.
        :param args: Additional arguments to pass to the function after the item.
        :param error_msg: Error message prefix for logging.
        :param kwargs: Keyword arguments to pass to the function.
        :return: Number of successfully processed items.
        """
        if not items:
            return 0

        success_count = 0
        for item in items:
            try:
                func(item, *args, **kwargs)
                success_count += 1
            except Exception as e:
                logger.warning("%s: %s", error_msg, e)
        return success_count

    def process_row(self, row: pd.Series, i: int, datasources: pd.DataFrame) -> None:
        """
        Process a single row of the DataFrame and update the RDF graph.

        :param row: A dictionary-like object representing a single row of the DataFrame.
        :param i: An integer representing the index of the row.
        :param datasources: The BDF datasource table.
        """
        # Determine entity type and create node
        gene_node, compound_node = None, None
        is_gene, is_compound = False, False

        identifier_source = row.get("identifier.source")
        gene_sources = list(datasources[datasources["type"] == "gene"]["source"]) + ["Entrez Gene"]
        compound_sources = list(datasources[datasources["type"] == "metabolite"]["source"])

        if identifier_source in gene_sources:
            gene_node = self.add_gene_node(row)
            is_gene = True
        elif identifier_source in compound_sources:
            compound_node = self.add_compound_node(row)
            is_compound = True

        # Validate row
        source_idx = row.get(Cons.IDENTIFIER_COL)
        source_namespace = row.get(Cons.IDENTIFIER_SOURCE_COL)
        target_idx = row.get(Cons.TARGET_COL)
        target_namespace = row.get(Cons.TARGET_SOURCE_COL)

        if not self.valid_indices(source_idx, source_namespace, target_idx, target_namespace):
            return

        # Map BridgeDB namespaces to bioregistry-compatible prefixes
        source_prefix = Cons.BRIDGEDB_TO_BIOREGISTRY.get(source_namespace, source_namespace)
        target_prefix = Cons.BRIDGEDB_TO_BIOREGISTRY.get(target_namespace, target_namespace)

        # Try to normalize CURIEs but don't fail if normalization doesn't work
        source_curie = normalize_curie(f"{source_prefix}:{source_idx}")
        target_curie = normalize_curie(f"{target_prefix}:{target_idx}")
        if not source_curie:
            logger.debug("Could not normalize source CURIE: %s:%s", source_prefix, source_idx)
        if not target_curie:
            logger.debug("Could not normalize target CURIE: %s:%s", target_prefix, target_idx)

        id_number = f"{i:06d}"

        # Get protein nodes for genes
        protein_nodes = []
        if is_gene and gene_node:
            try:
                protein_nodes = list(
                    self.objects(gene_node, URIRef(Cons.PREDICATES["translation_of"]))
                )
            except Exception as e:
                logger.warning("Failed to get protein nodes: %s", e)

        # Process gene data
        if is_gene and gene_node:
            self._process_gene_data(row, gene_node, protein_nodes, id_number, source_idx, i)

        # Process compound data
        if is_compound and compound_node:
            self._process_compound_data(row, compound_node, id_number)

    def _process_gene_data(
        self,
        row: pd.Series,
        gene_node: URIRef,
        protein_nodes: List[URIRef],
        id_number: str,
        source_idx: str,
        row_index: int,
    ) -> None:
        """Process all gene-related data types."""
        # PPI data
        ppi_data = row.get(Cons.STRING_INTERACT_COL)
        if ppi_data:
            self._safe_process(
                self.process_ppi_data, ppi_data, gene_node, error_msg="Failed to process PPI data"
            )

        # Disease data
        disease_data = self.collect_disease_data(row)
        if disease_data:
            self._safe_process(
                self.process_disease_data,
                disease_data,
                id_number,
                source_idx,
                gene_node,
                error_msg="Failed to process disease data",
            )

        # Expression data
        expression_data = row.get(Cons.BGEE_GENE_EXPRESSION_LEVELS_COL)
        if expression_data:
            self._safe_process(
                self.process_expression_data,
                expression_data,
                id_number,
                source_idx,
                gene_node,
                error_msg="Failed to process expression data",
            )

        # Pathways
        self._safe_process(
            self.process_pathways,
            row,
            gene_node,
            protein_nodes,
            error_msg="Failed to process pathways",
        )

        # GO terms
        processes_data = row.get(Cons.OPENTARGETS_GO_COL)
        if processes_data:
            self._safe_process(
                self.process_processes_data,
                processes_data,
                gene_node,
                error_msg="Failed to process GO terms",
            )

        # Compounds
        compound_data = row.get(Cons.OPENTARGETS_GENE_COMPOUND_COL)
        if compound_data:
            self._safe_process(
                self.process_compound_data,
                compound_data,
                gene_node,
                error_msg="Failed to process compound data",
            )

        # Literature
        literature_data = row.get(Cons.LITERATURE_DISEASE_COL)
        if literature_data:
            self._safe_process(
                self.process_literature_data,
                literature_data,
                gene_node,
                id_number,
                source_idx,
                self.new_uris,
                row_index,
                error_msg="Failed to process literature data",
            )

        # Transporter inhibitors
        ti_data = row.get(Cons.MOLMEDB_PROTEIN_COMPOUND_COL)
        if ti_data:
            self._safe_process(
                self.process_transporter_inhibitor_data,
                gene_node,
                ti_data,
                error_msg="Failed to process transporter inhibitor data",
            )

        # Protein variants (not yet implemented)
        # if self.include_variants and protein_nodes:
        #     self._safe_process(
        #         self.process_protein_variants, protein_nodes,
        #         error_msg="Failed to process protein variants"
        #     )

        # AOP data
        aop_data = row.get(Cons.AOPWIKI_GENE_COL)
        if aop_data:
            self._safe_process(
                self.process_aop_data,
                aop_data,
                gene_node,
                None,
                error_msg="Failed to process AOP data",
            )

        # Molecular pathways
        pathways_data = row.get(Cons.WIKIPATHWAYS_MOLECULAR_COL)
        if pathways_data:
            self._safe_process(
                self.process_molecular_pathway,
                pathways_data,
                gene_node,
                id_number,
                error_msg="Failed to process molecular pathway",
            )

        # PubChem assay data
        pubchem_assays = row.get(Cons.PUBCHEM_COMPOUND_ASSAYS_COL)
        if pubchem_assays:
            self._safe_process(
                self.process_pubchem_assay_data,
                pubchem_assays,
                gene_node,
                error_msg="Failed to process PubChem assay data",
            )

        # CompoundWiki
        self._safe_process(
            self.process_compoundwiki_data,
            gene_node,
            row,
            error_msg="Failed to process CompoundWiki data",
        )

    def _process_compound_data(self, row: pd.Series, compound_node: URIRef, id_number: str) -> None:
        """Process all compound-related data types."""
        # Pathways
        self._safe_process(
            self.process_pathways, row, compound_node, [], error_msg="Failed to process pathways"
        )

        # Inhibitor transporter
        it_data = row.get(Cons.MOLMEDB_COMPOUND_PROTEIN_COL)
        if it_data:
            self._safe_process(
                self.process_inhibitor_transporter_data,
                compound_node,
                it_data,
                error_msg="Failed to process inhibitor transporter data",
            )

        # AOP data
        aop_data = row.get(Cons.AOPWIKI_COMPOUND_COL)
        if aop_data:
            self._safe_process(
                self.process_aop_data,
                aop_data,
                None,
                compound_node,
                error_msg="Failed to process AOP data",
            )

        # Molecular pathways
        pathways_data = row.get(Cons.WIKIPATHWAYS_MOLECULAR_COL)
        if pathways_data:
            self._safe_process(
                self.process_molecular_pathway,
                pathways_data,
                compound_node,
                id_number,
                error_msg="Failed to process molecular pathway",
            )

        # CompoundWiki
        self._safe_process(
            self.process_compoundwiki_data,
            compound_node,
            row,
            error_msg="Failed to process CompoundWiki data",
        )

    # =========================================================================
    # Data collection methods
    # =========================================================================

    def collect_disease_data(self, row: pd.Series) -> List[Dict[str, Any]]:
        """
        Collect disease data from the row.

        :param row: A dictionary representing a row of data.
        :return: A list of collected disease data.
        """
        disease_data = []
        for source_col in [Cons.DISGENET_DISEASE_COL, Cons.OPENTARGETS_DISEASE_COL]:
            source_data = row.get(source_col, None)
            if source_data is not None:
                disease_data.extend(source_data)
        return disease_data

    def valid_indices(
        self,
        source_idx: Optional[str],
        source_namespace: Optional[str],
        target_idx: Optional[str],
        target_namespace: Optional[str],
    ) -> bool:
        """
        Check if the row is valid.

        This method verifies that none of the provided indices or namespaces are NaN (Not a Number).

        :param source_idx: The index of the source node.
        :param source_namespace: The namespace of the source node.
        :param target_idx: The index of the target node.
        :param target_namespace: The namespace of the target node.
        :return: True if all indices and namespaces are valid (not NaN), False otherwise.
        """
        return not any(
            pd.isna(val) for val in [source_idx, source_namespace, target_idx, target_namespace]
        )

    def add_gene_node(self, row: pd.Series) -> Optional[URIRef]:
        """
        Get gene node.

        Dynamically creates gene nodes based on the target source.
        The target source must have a corresponding entry in NODE_URI_PREFIXES.

        :param row: A series containing the data for a single row.
                    It must include the key "target.source".
        :return: A URIRef for the gene, else None.
        """
        target_source = row.get("target.source")
        # Check if we have a URI prefix for this source
        if target_source and target_source in Cons.NODE_URI_PREFIXES:
            return get_gene_node(self, row)
        return None

    def add_compound_node(self, row: pd.Series) -> Optional[URIRef]:
        """
        Get compound node.

        :param row: A series containing the data for a single row.
                    It must include the key "target.source".
        :return: A URIRef for the compound, else None.
        """
        return get_compound_node(self, row)

    def process_disease_data(
        self, disease_data: List[Dict[str, Any]], id_number: str, source_idx: str, gene_node: URIRef
    ) -> None:
        """
        Process disease data and add to the RDF graph.

        :param disease_data: List of disease data to be processed.
        :param id_number: Identifier number for the gene.
        :param source_idx: Source index for the data.
        :param gene_node: RDF node representing the gene.
        """
        if not disease_data:
            return

        # Track which data sources contributed disease data
        sources_used = set()
        for j, disease in enumerate(disease_data):
            try:
                add_gene_disease_associations(
                    self, id_number, source_idx, gene_node, disease, self.new_uris, j
                )
                # Track data source based on disease data structure
                if disease.get("disgenet_gene_disease_score"):
                    sources_used.add(Cons.DISGENET)
                elif disease.get("opentargets_disease_assoc_score"):
                    sources_used.add(Cons.OPENTARGETS)
            except Exception as e:
                logger.warning("Failed to process disease entry %d: %s", j, e)

        for source in sources_used:
            self.record_datasource(
                source, node=gene_node, interaction_type="gene-disease association"
            )

    def process_expression_data(
        self, expression_data, id_number: str, source_idx: str, gene_node: URIRef
    ) -> None:
        """
        Process gene expression data and add to the RDF graph.

        :param expression_data: The gene expression and experimental process data.
        :param id_number: The identifier number for the gene.
        :param source_idx: The source index for the data.
        :param gene_node: The RDF node representing the gene.
        """
        if expression_data:
            add_gene_expression_data(
                self,
                id_number,
                source_idx,
                gene_node,
                expression_data,
                expression_data,
                self.new_uris,
            )
            self.record_datasource(Cons.BGEE, node=gene_node, interaction_type="gene expression")

    def process_processes_data(
        self, processes_data: Optional[List[Dict[str, Any]]], gene_node: URIRef
    ) -> None:
        """
        Process Gene Ontology (GO) terms and add to the RDF graph.

        :param processes_data: A list of GO terms related to a gene.
        :param gene_node: The RDF node representing the gene.
        """
        if not processes_data:
            return

        data_added = False
        for process_data in processes_data:
            try:
                go_cpf = add_go_cpf(self, process_data)
                if go_cpf:
                    self.add((gene_node, URIRef(Cons.PREDICATES["sio_is_part_of"]), go_cpf))
                    self.add((go_cpf, URIRef(Cons.PREDICATES["sio_has_part"]), gene_node))
                    data_added = True
            except Exception as e:
                logger.warning("Failed to process GO term: %s", e)

        if data_added:
            self.record_datasource(
                Cons.OPENTARGETS, node=gene_node, interaction_type="GO annotation"
            )

    def process_compound_data(
        self, compound_data: Optional[List[Dict[str, Any]]], gene_node: URIRef
    ) -> None:
        """
        Process compound data and add to the RDF graph.

        :param compound_data: List of compounds to be processed.
        :param gene_node: URIRef of gene node.
        """
        if not compound_data:
            return

        data_added = False
        for compound in compound_data:
            try:
                add_associated_compound_node(self, compound, gene_node)
                data_added = True
            except Exception as e:
                logger.warning("Failed to process compound: %s", e)

        if data_added:
            self.record_datasource(
                Cons.OPENTARGETS, node=gene_node, interaction_type="gene-compound interaction"
            )

    def process_pubchem_assay_data(
        self, assay_data: List[Dict[str, Any]], gene_node: URIRef
    ) -> None:
        """
        Process PubChem assay data and add to the RDF graph.

        :param assay_data: List of PubChem assay entries to be processed.
        :param gene_node: URIRef of gene node.
        """
        if not assay_data:
            logger.warning("No assay data")
            return

        data_added = False
        for assay in assay_data:
            try:
                if isinstance(assay, dict):
                    assay_node = add_pubchem_assay_node(self, assay, gene_node)
                    if assay_node:
                        data_added = True
            except Exception as e:
                logger.warning("Failed to process PubChem assay: %s", e)

        if data_added:
            self.record_datasource(
                Cons.PUBCHEM, node=gene_node, interaction_type="gene-compound assay"
            )

    def process_compoundwiki_data(self, target_node: URIRef, row: pd.Series) -> None:
        """
        Process CompoundWiki annotation data and add to the RDF graph.

        :param target_node: URIRef of the target node (gene or protein).
        :param row: Data row containing CompoundWiki annotations.
        """
        data_added = False

        # Check if we have PubChem assays with CompoundWiki annotations
        pubchem_assays = row.get(Cons.PUBCHEM_COMPOUND_ASSAYS_COL, None)
        if pubchem_assays and isinstance(pubchem_assays, list):
            for assay in pubchem_assays:
                try:
                    if isinstance(assay, dict):
                        compoundwiki_annotations = assay.get(Cons.COMPOUNDWIKI_COL, None)
                        if compoundwiki_annotations:
                            add_compoundwiki_annotations(
                                self,
                                target_node,
                                [compoundwiki_annotations],
                            )
                            data_added = True
                except Exception as e:
                    logger.warning("Failed to process CompoundWiki from PubChem assay: %s", e)

        # Check for CompoundWiki annotations in OpenTargets compound data
        opentargets_compounds = row.get(Cons.OPENTARGETS_GENE_COMPOUND_COL, None)
        if opentargets_compounds and isinstance(opentargets_compounds, list):
            for compound in opentargets_compounds:
                try:
                    if isinstance(compound, dict):
                        compoundwiki_annotations = compound.get(Cons.COMPOUNDWIKI_COL, None)
                        if compoundwiki_annotations:
                            add_compoundwiki_annotations(
                                self,
                                target_node,
                                [compoundwiki_annotations],
                            )
                            data_added = True
                except Exception as e:
                    logger.warning("Failed to process CompoundWiki from OpenTargets: %s", e)

        # Check for CompoundWiki annotations in IntAct interaction data
        for intact_col in [Cons.INTACT_INTERACT_COL, Cons.INTACT_COMPOUND_INTERACT_COL]:
            intact_data = row.get(intact_col, None)
            if intact_data and isinstance(intact_data, list):
                for interaction in intact_data:
                    try:
                        if isinstance(interaction, dict):
                            compoundwiki_annotations = interaction.get(Cons.COMPOUNDWIKI_COL, None)
                            if compoundwiki_annotations:
                                add_compoundwiki_annotations(
                                    self,
                                    target_node,
                                    [compoundwiki_annotations],
                                )
                                data_added = True
                    except Exception as e:
                        logger.warning("Failed to process CompoundWiki from IntAct: %s", e)

        # Check for CompoundWiki annotations in KEGG pathway data
        kegg_pathways = row.get(Cons.KEGG_PATHWAY_COL, None)
        if kegg_pathways and isinstance(kegg_pathways, list):
            for pathway in kegg_pathways:
                try:
                    if isinstance(pathway, dict):
                        compoundwiki_annotations = pathway.get(Cons.COMPOUNDWIKI_COL, None)
                        if compoundwiki_annotations:
                            add_compoundwiki_annotations(
                                self,
                                target_node,
                                [compoundwiki_annotations],
                            )
                            data_added = True
                except Exception as e:
                    logger.warning("Failed to process CompoundWiki from KEGG: %s", e)

        # Check for CompoundWiki annotations in MolMeDB data
        for molmedb_col in [Cons.MOLMEDB_PROTEIN_COMPOUND_COL, Cons.MOLMEDB_COMPOUND_PROTEIN_COL]:
            molmedb_data = row.get(molmedb_col, None)
            if molmedb_data and isinstance(molmedb_data, list):
                for compound in molmedb_data:
                    try:
                        if isinstance(compound, dict):
                            compoundwiki_annotations = compound.get(Cons.COMPOUNDWIKI_COL, None)
                            if compoundwiki_annotations:
                                add_compoundwiki_annotations(
                                    self,
                                    target_node,
                                    [compoundwiki_annotations],
                                )
                                data_added = True
                    except Exception as e:
                        logger.warning("Failed to process CompoundWiki from MolMeDB: %s", e)

        if data_added:
            self.record_datasource(
                Cons.COMPOUNDWIKI, node=target_node, interaction_type="compound annotation"
            )

    def process_literature_data(
        self,
        literature_based_data: Optional[Union[Dict[str, Any], List[Dict[str, Any]]]],
        gene_node: URIRef,
        id_number: str,
        source_idx: str,
        new_uris: Dict[str, str],
        i: int,
    ) -> None:
        """
        Process literature-based data and add to the RDF graph.

        :param literature_based_data: Data derived from literature sources. Can be a single entry or a list of entries.
        :param gene_node: The gene node to which the literature-based data will be added.
        :param id_number: Unique identifier for the expression data.
        :param source_idx: Identifier for the source of the expression data.
        :param new_uris: Node URIs for the graph.
        :param i: Row index.
        """
        if literature_based_data:
            data_added = False
            entries = (
                literature_based_data
                if isinstance(literature_based_data, list)
                else [literature_based_data]
            )
            for entry in entries:
                try:
                    if entry.get(Cons.UMLS, None):
                        umls_parts = entry[Cons.UMLS].split(":")
                        umlscui = umls_parts[1] if len(umls_parts) > 1 else None
                        if not umlscui:
                            continue
                        disease_data_lit = {
                            Cons.UMLS: umlscui,
                            Cons.DISGENET_SCORE: None,
                            Cons.DISGENET_EI: None,
                            Cons.DISGENET_EL: None,
                            Cons.DISEASE_NAME: entry[Cons.DISEASE_NAME],
                        }
                        add_literature_based_data(
                            self,
                            entry,
                            gene_node,
                            id_number,
                            disease_data_lit,
                            source_idx,
                            new_uris,
                            i,
                        )
                        data_added = True
                except Exception as e:
                    logger.warning("Failed to process literature entry: %s", e)
            if data_added:
                self.record_datasource(
                    Cons.PUBCHEM, node=gene_node, interaction_type="literature-based association"
                )

    def process_transporter_inhibitor_data(
        self, gene_node, transporter_inhibitor_data: Optional[List[Dict[str, Any]]]
    ) -> None:
        """
        Process transporter inhibitor data and add to the RDF graph.

        :param gene_node: An RDF node representing the gene.
        :param transporter_inhibitor_data: A list of transporter inhibitor data entries to be processed.
        """
        if transporter_inhibitor_data:
            data_added = False
            for entry in transporter_inhibitor_data:
                try:
                    add_transporter_inhibitor_node(self, gene_node, entry, self.base_uri)
                    data_added = True
                except Exception as e:
                    logger.warning("Failed to process transporter inhibitor entry: %s", e)
            if data_added:
                self.record_datasource(
                    Cons.MOLMEDB, node=gene_node, interaction_type="transporter inhibition"
                )

    def process_inhibitor_transporter_data(
        self, compound_node, inhibitor_transporter_data: Optional[List[Dict[str, Any]]]
    ) -> None:
        """
        Process inhibitor transporter data and add to the RDF graph.

        :param compound_node: An RDF node representing the compound.
        :param inhibitor_transporter_data: A list of inhibitor transporter data entries to be processed.
        """
        if inhibitor_transporter_data:
            data_added = False
            for entry in inhibitor_transporter_data:
                try:
                    add_inhibitor_transporter_node(self, compound_node, entry, self.base_uri)
                    data_added = True
                except Exception as e:
                    logger.warning("Failed to process inhibitor transporter entry: %s", e)
            if data_added:
                self.record_datasource(
                    Cons.MOLMEDB, node=compound_node, interaction_type="transporter inhibition"
                )

    def process_pathways(
        self, row: pd.Series, identifier_node: URIRef, protein_nodes: List[URIRef]
    ) -> None:
        """
        Process pathway data and add to the RDF graph.

        This method processes pathway data from various sources and adds the relevant
        information to the RDF graph. It creates pathway nodes and establishes relationships
        between gene nodes, protein nodes, and pathway nodes.

        :param row: A dictionary containing pathway data from different sources.
        :param identifier_node: An RDF node representing the identifier.
        :param protein_nodes: A list of RDF nodes representing proteins associated with the identifier.
        """
        # Map column names to constants
        source_constants = {
            "WikiPathways": Cons.WIKIPATHWAYS,
            "MINERVA": Cons.MINERVA,
            "OpenTargets": Cons.OPENTARGETS,
        }
        for source in ["WikiPathways", "MINERVA", "OpenTargets"]:
            pathway_data_list = row.get(source)
            if pathway_data_list:
                data_added = False
                for pathway_data in pathway_data_list:
                    try:
                        if pathway_data.get("pathway_id"):
                            pathway_node = add_pathway_node(self, pathway_data, source)
                            if pathway_node:
                                self.add(
                                    (
                                        identifier_node,
                                        URIRef(Cons.PREDICATES["sio_is_part_of"]),
                                        pathway_node,
                                    )
                                )
                                self.add(
                                    (
                                        pathway_node,
                                        URIRef(Cons.PREDICATES["sio_has_part"]),
                                        identifier_node,
                                    )
                                )
                                if protein_nodes:
                                    for protein_node in protein_nodes:
                                        self.add(
                                            (
                                                protein_node,
                                                URIRef(Cons.PREDICATES["sio_is_part_of"]),
                                                pathway_node,
                                            )
                                        )
                                        self.add(
                                            (
                                                pathway_node,
                                                URIRef(Cons.PREDICATES["sio_has_part"]),
                                                protein_node,
                                            )
                                        )
                                self.add(
                                    (
                                        pathway_node,
                                        URIRef(Cons.PREDICATES["sio_has_part"]),
                                        identifier_node,
                                    )
                                )
                                self.add(
                                    (
                                        pathway_node,
                                        URIRef(Cons.PREDICATES["sio_has_source"]),
                                        URIRef(Cons.DATA_SOURCES[source]),
                                    )
                                )
                                data_added = True
                    except Exception as e:
                        logger.warning("Failed to process pathway entry from %s: %s", source, e)
                if data_added:
                    self.record_datasource(
                        source_constants[source],
                        node=identifier_node,
                        interaction_type="pathway membership",
                    )

    def process_molecular_pathway(self, molecular_data, identifier, id_number) -> None:
        """
        Process molecular pathway data and add to the RDF graph.

        :param molecular_data: A list of dicts containing pathway data.
        :param identifier: An RDF node representing the gene or compound in the row.
        :param id_number: The identifier number for the row.
        """
        if not molecular_data:
            return
        data_added = False
        for el in molecular_data:
            try:
                add_molecular_pathway_node(self, el, identifier, id_number)
                data_added = True
            except Exception as e:
                logger.warning("Failed to process molecular pathway entry: %s", e)
        if data_added:
            self.record_datasource(
                Cons.WIKIPATHWAYS, node=identifier, interaction_type="molecular pathway"
            )

    def process_ppi_data(
        self, stringdb_data: Optional[List[Dict[str, Any]]], gene_node: URIRef
    ) -> None:
        """
        Process Protein-Protein Interaction (ppi) data and add to the RDF graph.

        :param stringdb_data: List of dictionaries containing ppi data from STRING database.
        :param gene_node: The gene URIRef.
        """
        if not stringdb_data:
            return

        data_added = False
        for entry in stringdb_data:
            try:
                if entry.get("Ensembl"):
                    result = add_ppi_data(g=self, gene_node=gene_node, entry=entry)
                    if result:
                        data_added = True
            except Exception as e:
                logger.warning("Failed to process PPI entry: %s", e)

        if data_added:
            self.record_datasource(
                Cons.STRING, node=gene_node, interaction_type="protein-protein interaction"
            )

    def process_aop_data(
        self,
        aop_data: Optional[List[Dict[str, Any]]] = None,
        gene_node: Optional[URIRef] = None,
        compound_node: Optional[URIRef] = None,
    ) -> None:
        """
        Process AOP-Wiki data and add to the RDF graph.

        :param aop_data: List of dictionaries containing AOP data. Defaults to None.
        :param gene_node: The gene URIRef. Defaults to None.
        :param compound_node: The compound URIRef. Defaults to None.
        """
        if not aop_data:
            return

        data_added = False

        if gene_node:
            for entry in aop_data:
                try:
                    if entry.get("aop"):
                        add_aop_data(
                            g=self,
                            gene_node=gene_node,
                            compound_node=None,
                            entry=entry,
                        )
                        data_added = True
                except Exception as e:
                    logger.warning("Failed to process AOP entry for gene: %s (entry: %s)", e, entry)

        if compound_node:
            for entry in aop_data:
                try:
                    add_aop_data(
                        g=self,
                        gene_node=None,
                        compound_node=compound_node,
                        entry=entry,
                    )
                    data_added = True
                except Exception as e:
                    logger.warning("Failed to process AOP entry for compound: %s", e)

        if data_added:
            target_node = gene_node or compound_node
            self.record_datasource(
                Cons.AOPWIKIRDF, node=target_node, interaction_type="AOP relationship"
            )

    # Other methods not related to adding nodes begin here
    def _add_metadata(self, metadata: Dict[str, Any]) -> None:
        """Add metadata to the RDF graph.

        :param metadata: Dataframe of BDF metadata to be added.
        """
        # Use version_iri or fall back to base_uri for graph_uri
        graph_uri = self.version_iri if self.version_iri else self.base_uri

        add_metadata(
            g=self,
            graph_uri=graph_uri,
            metadata=metadata,
            version_iri=self.version_iri,
            title=self.title,
            description=self.description,
            author=self.author,
            orcid=self.orcid,
            creators=self.creators,
        )

    def shex(
        self,
        path: Optional[str] = None,
        threshold: float = 0.001,
        uml_figure_path: Optional[str] = None,
        print_string_output: bool = True,
        additional_namespaces: Optional[Dict[str, str]] = None,
    ) -> Any:
        """Get ShEx shapes with optional parameters.

        :param path: Path to save the ShEx results.
        :param threshold: Validation threshold.
        :param uml_figure_path: Path to save UML diagram for shapes.
        :param print_string_output: Whether to print the output string.
        :param additional_namespaces: Additional namespaces for shapes.
        :return: ShEx graph result.
        """
        return get_shapes(
            self,
            self.base_uri,
            path,
            threshold,
            "shex",
            uml_figure_path,
            print_string_output,
            additional_namespaces,
        )

    def shacl(
        self,
        path: Optional[str] = None,
        threshold: float = 0.001,
        uml_figure_path: Optional[str] = None,
        print_string_output: bool = True,
        additional_namespaces: Optional[Dict[str, str]] = None,
    ) -> Any:
        """Get SHACL shapes with optional parameters.

        :param path: Path to save the SHACL results.
        :param threshold: Validation threshold.
        :param uml_figure_path: Path to save UML diagram for shapes.
        :param print_string_output: Whether to print the output string.
        :param additional_namespaces: Additional namespaces for shapes.
        :return: SHACL graph result.
        """
        return get_shapes(
            self,
            self.base_uri,
            path,
            threshold,
            "shacl",
            uml_figure_path,
            print_string_output,
            additional_namespaces,
        )

    def shacl_prefixes(
        self,
        path: Optional[str] = None,
        namespaces: Optional[Dict[str, str]] = None,
        print_string_output: bool = False,
    ) -> Any:
        """Get a SHACL prefixes graph, optionally add more namespaces to bind to it.

        :param path: Path to save the SHACL prefixes.
        :param namespaces: Namespaces for the prefixes.
        :param print_string_output: bool, print or not the generated TTL as a string.
        :return: SHACL prefixes.
        """
        output_path = path if path is not None else self._prefixes_path
        current_namespaces = self._namespaces
        if namespaces:
            current_namespaces = current_namespaces or {}
            current_namespaces.update(namespaces)

        return get_shacl_prefixes(
            namespaces=current_namespaces,
            path=output_path,
            new_uris=self.new_uris,
            print_string_output=print_string_output,
        )
