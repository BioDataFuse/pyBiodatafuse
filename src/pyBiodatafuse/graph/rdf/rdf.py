"""
Module for generating and managing RDF graphs within the BioDatafuse framework.

This module defines the `BDFGraph` class, which extends the functionality of a basic
`rdflib Graph` to include specific operations for generating RDF graphs from BDF data.

The class supports adding nodes for genes, proteins, pathways,
compounds, gene-disease associations, gene expression data, and
more. It also uses shexer for SHACL and ShEx elucidation, and metadata
management for the BDF graph.

Classes:
    BDFGraph(Graph): Main class for constructing and managing BioDatafuse RDF graphs.
    - `process_row`: Processes a single row of the DataFrame and updates the RDF graph.
    - `collect_disease_data`: Collects disease data from a row.
    - `valid_indices`: Checks if the required indices and namespaces are valid.
    - `get_gene_node`: Gets gene node based on row data.
    - `process_disease_data`: Processes disease data, adds it to RDF graph.
    - `process_expression_data`: Processes gene expression data, adds it to RDF graph.
    - `process_pathways`: Processes pathway data, adds it to RDF graph.
    - `process_processes_data`: Processes Gene Ontology (GO) terms and adds them to the RDF graph.
    - `process_compound_data`: Processes compound data, adds it to RDF graph.
    - `process_literature_data`: Processes literature-based data, adds it to RDF graph.
    - `process_transporter_inhibitor_data`: Processes transporter-inhibitor data, adds it to RDF graph.
    - `process_protein_variants`: Processes protein variants and adds them to the RDF graph.
    - `process_ppi_data`: Processes Protein-Protein Interaction (ppi) data, adds it to RDF graph.
    - `process_aop_data`: Processes Protein-Protein Interaction (ppi) data, adds it to RDF graph.
    - `_add_metadata`: Attaches metadata to the RDF graph.
    - `shex`: Runs shexer on the RDF graph to obtain its ShEx shapes.
    - `shacl`: Runs shexer on the RDF graph to obtain its SHACL shapes.
    - `shacl_prefixes`: Retrieves SHACL prefixes for the graph.
"""

import logging
from typing import Any, Dict, List, Optional, Union

import pandas as pd
from bioregistry import normalize_curie
from rdflib import Graph, Literal, URIRef
from rdflib.namespace import RDF, RDFS, XSD
from tqdm import tqdm

import pyBiodatafuse.constants as Cons
from pyBiodatafuse.graph.rdf.metadata import add_metadata
from pyBiodatafuse.graph.rdf.nodes.aop import add_aop_data
from pyBiodatafuse.graph.rdf.nodes.compound import (
    add_associated_compound_node,
    add_inhibitor_transporter_node,
    add_transporter_inhibitor_node,
    get_compound_node,
)
from pyBiodatafuse.graph.rdf.nodes.gene import get_gene_node
from pyBiodatafuse.graph.rdf.nodes.gene_disease import add_gene_disease_associations
from pyBiodatafuse.graph.rdf.nodes.gene_expression import add_gene_expression_data
from pyBiodatafuse.graph.rdf.nodes.go_terms import add_go_cpf
from pyBiodatafuse.graph.rdf.nodes.literature import add_literature_based_data
from pyBiodatafuse.graph.rdf.nodes.pathway import add_molecular_pathway_node, add_pathway_node
from pyBiodatafuse.graph.rdf.nodes.protein_protein import add_ppi_data
from pyBiodatafuse.graph.rdf.utils import get_shacl_prefixes, get_shapes, replace_na_none
from pyBiodatafuse.id_mapper import read_datasource_file

# Set up logger
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

# Create console handler with a higher log level
console_handler = logging.StreamHandler()
console_handler.setLevel(logging.INFO)

# Create formatter and add it to the handler
formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
console_handler.setFormatter(formatter)

# Add the handler to the logger
logger.addHandler(console_handler)


class BDFGraph(Graph):
    """Main class for a BioDatafuse RDF Graph, superclass of rdflib.Graph."""

    def __init__(self, base_uri: str, version_iri: str, author: str, orcid: str):
        """
        Initialize a new instance of the class with the provided metadata and URIs.

        :param base_uri: The base URI for the RDF graph.
        :param version_iri: The version IRI for the RDF graph.
        :param author: The author of the BDF graph.
        :param orcid: The ORCID identifier for the author.
        """
        # Initialize the rdflib.Graph superclass without passing extra arguments
        super().__init__()

        # Assign parameters to instance attributes
        self.base_uri = base_uri
        self.version_iri = version_iri
        self.author = author
        self.orcid = orcid

        # Create and bind custom URIs and namespaces
        self.new_uris = {key: self.base_uri + value for key, value in Cons.URIS.items()}
        self._shex_path = None
        self._shacl_path = None
        self._prefixes_path = None
        self._namespaces = None
        self.include_variants = False  # TODO: Allow user to set options that can affect graph size

        # Bind prefixes
        for key, new_value in self.new_uris.items():
            self.bind(key, new_value)
        for key, value in Cons.NAMESPACE_BINDINGS.items():
            self.bind(key, value)

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
        df = df.applymap(replace_na_none)
        datasources = read_datasource_file()
        if not self.include_variants:
            df = df[df[Cons.TARGET_SOURCE_COL] == Cons.ENSEMBL]
        for j, (_, row) in enumerate(
            tqdm(df.iterrows(), total=df.shape[0], desc="Building RDF graph")
        ):
            self.process_row(row, j, datasources)
        self._add_metadata(metadata)

    def process_row(self, row: pd.Series, i: int, datasources: pd.DataFrame) -> None:
        """
        Process a single row of the DataFrame and update the RDF graph.

        :param row: A dictionary-like object representing a single row of the DataFrame.
        :param i: An integer representing the index of the row.
        :param datasources: The BDF datasource table.
        """
        # Initialize variables
        compound_node: Optional[URIRef] = None
        gene_node: Optional[URIRef] = None
        protein_nodes: List[URIRef] = []

        # Determine whether it's a gene or compound row
        gene, compound = False, False

        if row["identifier.source"] in list(
            datasources[datasources["type"] == "gene"]["source"]
        ) + [
            "Entrez Gene"
        ]:  # TODO fix datasources
            gene_node = self.add_gene_node(row)
            gene = True
        elif row["identifier.source"] in datasources[datasources["type"] == "metabolite"]["source"]:
            compound_node = self.add_compound_node(row)
            compound = True
        source_idx = row.get(Cons.IDENTIFIER_COL)
        source_namespace = row.get(Cons.IDENTIFIER_SOURCE_COL)
        target_idx = row.get(Cons.TARGET_COL)
        target_namespace = row.get(Cons.TARGET_SOURCE_COL)

        if not self.valid_indices(source_idx, source_namespace, target_idx, target_namespace):
            return

        source_curie = normalize_curie(f"{source_namespace}:{source_idx}")
        target_curie = normalize_curie(f"{target_namespace}:{target_idx}")

        if not source_curie or not target_curie:
            return

        id_number = f"{i:06d}"
        disease_data = self.collect_disease_data(row)

        # Extract relevant columns before processing
        string_ppi_data = row.get(Cons.STRING_INTERACT_COL, None)
        disease_data = self.collect_disease_data(row)
        expression_data = row.get(Cons.BGEE_GENE_EXPRESSION_LEVELS_COL, None)
        pathways_data = row.get(Cons.WIKIPATHWAYS_MOLECULAR_COL, None)
        processes_data = row.get(Cons.OPENTARGETS_GO_COL, None)
        compound_data = row.get(Cons.OPENTARGETS_GENE_COMPOUND_COL, None)
        literature_data = row.get(Cons.LITERATURE_DISEASE_COL, None)
        transporter_inhibitor_data = row.get(Cons.MOLMEDB_PROTEIN_COMPOUND_COL, None)
        inhibitor_transporter_data = row.get(Cons.MOLMEDB_COMPOUND_PROTEIN_COL, None)
        aop_data_gene = row.get(Cons.AOPWIKI_GENE_COL, None)
        aop_data_compound = row.get(Cons.AOPWIKI_COMPOUND_COL, None)

        if gene:
            try:
                self.process_ppi_data(string_ppi_data, gene_node)
            except Exception as e:
                logger.warning("Failed to process PPI data for gene node: %s", e)

            try:
                protein_nodes = list(
                    self.objects(gene_node, URIRef(Cons.PREDICATES["translation_of"]))
                )
            except Exception as e:
                logger.warning("Failed to retrieve protein nodes for gene node: %s", e)
                protein_nodes = []

            try:
                self.process_disease_data(disease_data, id_number, source_idx, gene_node)
            except Exception as e:
                logger.warning("Failed to process disease data for gene node: %s", e)

            try:
                self.process_expression_data(expression_data, id_number, source_idx, gene_node)
            except Exception as e:
                logger.warning("Failed to process expression data for gene node: %s", e)

            try:
                self.process_pathways(row, gene_node, protein_nodes)
            except Exception as e:
                logger.warning("Failed to process pathways for gene node: %s", e)

            try:
                self.process_processes_data(processes_data, gene_node)
            except Exception as e:
                logger.warning("Failed to process processes data for gene node: %s", e)

            try:
                self.process_compound_data(compound_data, gene_node)
            except Exception as e:
                logger.warning("Failed to process compound data for gene node: %s", e)

            try:
                self.process_literature_data(
                    literature_data, gene_node, id_number, source_idx, self.new_uris, i
                )
            except Exception as e:
                logger.warning("Failed to process literature data for gene node: %s", e)

            try:
                self.process_transporter_inhibitor_data(gene_node, transporter_inhibitor_data)
            except Exception as e:
                logger.warning("Failed to process transporter inhibitor data for gene node: %s", e)

            if self.include_variants:
                try:
                    self.process_protein_variants(protein_nodes)
                except Exception as e:
                    logger.warning("Failed to process protein variants for gene node: %s", e)

            try:
                self.process_aop_data(aop_data_gene, gene_node, None)
            except Exception as e:
                logger.warning("Failed to process AOP data for gene node: %s", e)

            try:
                self.process_molecular_pathway(pathways_data, gene_node, id_number)
            except Exception as e:
                logger.warning("Failed to process molecular pathway data for gene node: %s", e)

        if compound:
            try:
                self.process_pathways(row, compound_node, protein_nodes=[])
            except Exception as e:
                logger.warning("Failed to process pathways for compound node: %s", e)

            try:
                self.process_inhibitor_transporter_data(compound_node, inhibitor_transporter_data)
            except Exception as e:
                logger.warning(
                    "Failed to process inhibitor transporter data for compound node: %s", e
                )

            try:
                self.process_aop_data(aop_data_compound, None, compound_node)
            except Exception as e:
                logger.warning("Failed to process AOP data for compound node: %s", e)

            try:
                self.process_molecular_pathway(pathways_data, compound_node, id_number)
            except Exception as e:
                logger.warning("Failed to process molecular pathway data for compound node: %s", e)

    # Class methods about specific nodes begin here
    # If you add a new method, try to import most of the code from another script
    # Add new methods below
    def collect_disease_data(
        self,
        row: pd.Series,
    ) -> List[Dict[str, Any]]:
        """
        Collect disease data from the row.

        :param row: A dictionary representing a row of data.
        :return: A list of collected disease data.
        """
        disease_data = []
        for source_col in [Cons.DISGENET_DISEASE_COL, Cons.OPENTARGETS_DISEASE_COL]:
            # if open_only and source_col == DISGENET_DISEASE_COL:
            #     continue  # TODO fix open data only feature
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

        :param row: A series containing the data for a single row.
                    It must include the key "target.source".
        :return: A URIRef for the gene, else None.
        """
        if row["target.source"] == "Ensembl":
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
        if disease_data:
            for j, disease in enumerate(disease_data):
                add_gene_disease_associations(
                    self, id_number, source_idx, gene_node, disease, self.new_uris, j
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

    def process_processes_data(
        self, processes_data: Optional[List[Dict[str, Any]]], gene_node: URIRef
    ) -> None:
        """
        Process Gene Ontology (GO) terms and add to the RDF graph.

        :param processes_data: A list of GO terms related to a gene.
        :param gene_node: The RDF node representing the gene.
        """
        if processes_data:
            for process_data in processes_data:
                go_cpf = add_go_cpf(self, process_data)
                if go_cpf:
                    self.add((gene_node, URIRef(Cons.PREDICATES["sio_is_part_of"]), go_cpf))
                    self.add((go_cpf, URIRef(Cons.PREDICATES["sio_has_part"]), gene_node))

    def process_compound_data(
        self, compound_data: Optional[List[Dict[str, Any]]], gene_node: URIRef
    ) -> None:
        """
        Process compound data and add to the RDF graph.

        :param compound_data: List of compounds to be processed.
        :param gene_node: URIRef of gene node.
        """
        if compound_data is not None:
            for compound in compound_data:
                add_associated_compound_node(self, compound, gene_node)

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
            entries = (
                literature_based_data
                if isinstance(literature_based_data, list)
                else [literature_based_data]
            )
            for entry in entries:
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
                        self, entry, gene_node, id_number, disease_data_lit, source_idx, new_uris, i
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
            for entry in transporter_inhibitor_data:
                add_transporter_inhibitor_node(self, gene_node, entry, self.base_uri)

    def process_inhibitor_transporter_data(
        self, compound_node, inhibitor_transporter_data: Optional[List[Dict[str, Any]]]
    ) -> None:
        """
        Process inhibitor transporter data and add to the RDF graph.

        :param compound_node: An RDF node representing the compound.
        :param inhibitor_transporter_data: A list of inhibitor transporter data entries to be processed.
        """
        if inhibitor_transporter_data:
            for entry in inhibitor_transporter_data:
                add_inhibitor_transporter_node(self, compound_node, entry, self.base_uri)

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
        for source in ["WikiPathways", "MINERVA", "OpenTargets"]:
            pathway_data_list = row.get(source)
            if pathway_data_list:
                for pathway_data in pathway_data_list:
                    if pathway_data.get("pathway_id"):
                        pathway_node = add_pathway_node(self, pathway_data, source)
                        self.add(
                            (
                                identifier_node,
                                URIRef(Cons.PREDICATES["sio_is_part_of"]),
                                pathway_node,
                            )
                        )
                        self.add(
                            (pathway_node, URIRef(Cons.PREDICATES["sio_has_part"]), identifier_node)
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
                            (pathway_node, URIRef(Cons.PREDICATES["sio_has_part"]), identifier_node)
                        )
                        self.add(
                            (
                                pathway_node,
                                URIRef(Cons.PREDICATES["sio_has_source"]),
                                URIRef(Cons.DATA_SOURCES[source]),
                            )
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
        for el in molecular_data:
            add_molecular_pathway_node(self, el, identifier, id_number)

    def process_ppi_data(
        self, stringdb_data: Optional[List[Dict[str, Any]]], gene_node: URIRef
    ) -> None:
        """
        Process Protein-Protein Interaction (ppi) data and add to the RDF graph.

        :param stringdb_data: List of dictionaries containing ppi data from STRING database.
        :param gene_node: The gene URIRef.
        """
        if stringdb_data:
            for entry in stringdb_data:
                if entry.get("Ensembl"):
                    add_ppi_data(
                        g=self,
                        gene_node=gene_node,
                        entry=entry,
                        base_uri=self.base_uri,
                        new_uris=self.new_uris,
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
        if aop_data and gene_node:
            for entry in aop_data:
                if entry.get("Ensembl"):
                    add_aop_data(
                        g=self,
                        gene_node=gene_node,
                        compound_node=None,
                        entry=entry,
                    )
        if aop_data and compound_node:
            for entry in aop_data:
                if entry.get("Ensembl"):
                    add_aop_data(
                        g=self,
                        gene_node=None,
                        compound_node=compound_node,
                        entry=entry,
                    )

    # Other methods not related to adding nodes begin here
    def _add_metadata(self, metadata: Dict[str, Any]) -> None:
        """Add metadata to the RDF graph.

        :param metadata: Dataframe of BDF metadata to be added.
        """
        add_metadata(
            g=self,
            graph_uri=self.version_iri,  # TODO fix
            version_iri=self.version_iri,
            author=self.author,
            orcid=self.orcid,
            metadata=metadata,
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
        self, path: Optional[str] = None, namespaces: Optional[Dict[str, str]] = None
    ) -> Any:
        """Get a SHACL prefixes graph, optionally add more namespaces to bind to it.

        :param path: Path to save the SHACL prefixes.
        :param namespaces: Namespaces for the prefixes.
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
        )
