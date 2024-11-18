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

Functions:
    - `generate_rdf`: Converts a BDF DataFrame into an RDF graph with metadata.
    - `add_gene_protein_nodes`: Adds gene and protein nodes to the graph.
    - `add_gene_disease_associations`: Adds gene-disease associations to the graph.
    - `add_gene_expression_data`: Adds gene expression data to the graph.
    - `add_go_cpf`: Adds Gene Ontology (GO) terms to the graph.
    - `add_compound_node`: Adds compound data and links it to proteins.
    - `add_transporter_inhibitor_node`: Adds transporter-inhibitor data to the graph.
    - `add_pathway_node`: Adds pathway data to the RDF graph.
    - `add_literature_based_data`: Incorporates literature-based data into the graph.
    - `add_ppi_data`: Adds Protein-Protein Interaction (PPI) data to the graph.
    - `add_metadata`: Attaches metadata to the RDF graph.
    - `shex`: Runs shexer on the RDF graph to obtain its ShEx shapes.
    - `shacl`: Runs shexer on the RDF graph to obtain its SHACL shapes.
    - `shacl_prefixes`: Retrieves SHACL prefixes for the graph.
"""

import pandas as pd
from bioregistry import normalize_curie
from rdflib import Graph, URIRef

from pyBiodatafuse.constants import (
    BGEE_GENE_EXPRESSION_LEVELS_COL,
    DATA_SOURCES,
    DISGENET_DISEASE_COL,
    IDENTIFIER_COL,
    IDENTIFIER_SOURCE_COL,
    LITERATURE_DISEASE_COL,
    MOLMEDB_PROTEIN_COMPOUND_COL,
    NAMESPACE_BINDINGS,
    OPENTARGETS_DISEASE_COL,
    OPENTARGETS_GENE_COMPOUND_COL,
    OPENTARGETS_GO_COL,
    PREDICATES,
    PUBCHEM_COMPOUND_ASSAYS_COL,
    STRING_PPI_COL,
    TARGET_COL,
    TARGET_SOURCE_COL,
    URIS,
)
from pyBiodatafuse.graph.rdf.metadata import add_metadata
from pyBiodatafuse.graph.rdf.nodes.compound import add_compound_node, add_transporter_inhibitor_node
from pyBiodatafuse.graph.rdf.nodes.gene_disease import add_gene_disease_associations
from pyBiodatafuse.graph.rdf.nodes.gene_expression import add_gene_expression_data
from pyBiodatafuse.graph.rdf.nodes.gene_protein import add_gene_protein_nodes
from pyBiodatafuse.graph.rdf.nodes.go_terms import add_go_cpf
from pyBiodatafuse.graph.rdf.nodes.literature import add_literature_based_data
from pyBiodatafuse.graph.rdf.nodes.pathway import add_pathway_node
from pyBiodatafuse.graph.rdf.nodes.protein_protein import add_ppi_data
from pyBiodatafuse.graph.rdf.utils import get_shacl_prefixes, get_shapes, replace_na_none


class BDFGraph(Graph):
    """Main class for a BioDatafuse RDF Graph."""

    def __init__(self, base_uri: str, version_iri: str, author: str, orcid: str):
        """
        Initialize a new instance of the class with the provided metadata and URIs.

        :param base_uri: The base URI to be used for generating node URIs.
        :param version_iri: The IRI (Internationalized Resource Identifier) representing the version of the data.
        :param author: The name of the author of the graph or dataset.
        :param orcid: The ORCID (Open Researcher and Contributor ID) of the author.
        """
        super().__init__()
        self.base_uri = base_uri
        self.version_iri = version_iri
        self.author = author
        self.orcid = orcid
        self.new_uris = {key: self.base_uri + value for key, value in URIS.items()}
        self._shex_path = None
        self._shacl_path = None
        self._prefixes_path = None
        self._namespaces = None
        # Update self.bind only after constructing new URIs
        for key, new_value in self.new_uris.items():
            self.bind(key, new_value)

        for key, value in NAMESPACE_BINDINGS.items():
            self.bind(key, value)

    def generate_rdf(self, df: pd.DataFrame, metadata: dict, open_only: bool = False):
        """
        Generate an RDF graph from the provided DataFrame and metadata.

        :param df: The DataFrame containing biological data to be converted into an RDF graph.
        :param metadata: Metadata associated with the DataFrame, used for graph annotations.
        :param open_only: If True, only include data from open data sources (exclude non-open sources).
        """
        df = df.applymap(replace_na_none)
        for i, row in df.iterrows():
            # Unpack the relevant columns
            source_idx = row.get(IDENTIFIER_COL, None)
            source_namespace = row.get(IDENTIFIER_SOURCE_COL, None)
            target_idx = row.get(TARGET_COL, None)
            target_namespace = row.get(TARGET_SOURCE_COL, None)
            expression_data = row.get(BGEE_GENE_EXPRESSION_LEVELS_COL, None)
            experimental_process_data = row.get(PUBCHEM_COMPOUND_ASSAYS_COL, None)
            processes_data = row.get(OPENTARGETS_GO_COL, None)
            compound_data = row.get(OPENTARGETS_GENE_COMPOUND_COL, None)
            literature_based_data = row.get(LITERATURE_DISEASE_COL, None)
            transporter_inhibitor_data = row.get(MOLMEDB_PROTEIN_COMPOUND_COL, None)
            stringdb_data = row.get(STRING_PPI_COL, None)
            # molmedb_data = row.get("MolMeDB_transporter_inhibitor", None)
            # transporter_inhibited_data = row.get(MOLMEDB_COMPOUND_PROTEIN_COL, None)
            # Initialize an empty list to store disease data
            disease_data = []

            # Iterate over the specified data columns
            for source_col in [DISGENET_DISEASE_COL, OPENTARGETS_DISEASE_COL]:
                if open_only and source_col == DISGENET_DISEASE_COL:
                    continue  # TODO fix open data only feature
                # Extract the source data from the current row
                source_data = row.get(source_col, [])
                # Ensure source_data is not None and is of type array or list
                if len(source_data) > 0:
                    disease_data.extend(source_data)
            # Check for missing required indices or namespaces and skip the row if any are NaN
            if any(
                pd.isna(val) for val in [source_idx, source_namespace, target_idx, target_namespace]
            ):
                continue

            # Check for missing required indices or namespaces and skip the row if any are NaN
            if any(
                pd.isna(val) for val in [source_idx, source_namespace, target_idx, target_namespace]
            ):
                continue
            source_curie = normalize_curie(f"{source_namespace}:{source_idx}")
            target_curie = normalize_curie(f"{target_namespace}:{target_idx}")

            if not source_curie or not target_curie:
                continue

            id_number = f"{i:06d}"
            gene_node, protein_node = self.add_gene_protein_nodes(row)

            if len(disease_data) > 0:
                for j, disease in enumerate(disease_data):
                    self.add_gene_disease_associations(id_number, source_idx, gene_node, disease, j)
            if expression_data:
                self.add_gene_expression_data(
                    id_number,
                    source_idx,
                    gene_node,
                    expression_data,
                    experimental_process_data,
                )

            for source in ["WikiPathways", "MINERVA", "OpenTargets_reactome"]:
                if row.get(source, None):
                    for pathway_data in row[source]:
                        if pathway_data["pathway_id"]:
                            pathway_node = self.add_pathway_node(pathway_data, source)
                            self.add(
                                (gene_node, URIRef(PREDICATES["sio_is_part_of"]), pathway_node)
                            )
                            self.add(
                                (
                                    protein_node,
                                    URIRef(PREDICATES["sio_is_part_of"]),
                                    pathway_node,
                                )
                            )
                            self.add((pathway_node, URIRef(PREDICATES["sio_has_part"]), gene_node))
                            self.add(
                                (
                                    pathway_node,
                                    URIRef(PREDICATES["sio_has_source"]),
                                    URIRef(DATA_SOURCES[source]),
                                )
                            )

            if processes_data:
                for process_data in processes_data:
                    go_cpf = self.add_go_cpf(process_data)
                    if go_cpf:
                        self.add((gene_node, URIRef(PREDICATES["sio_is_part_of"]), go_cpf))
                        self.add((go_cpf, URIRef(PREDICATES["sio_has_part"]), gene_node))

            if compound_data:
                for compound in compound_data:
                    self.add_compound_node(compound, protein_node)

            if literature_based_data:
                if isinstance(literature_based_data, list):
                    entries = literature_based_data
                    for entry in entries:
                        self.add_literature_based_data(entry, gene_node)
                elif isinstance(literature_based_data, dict):
                    entry = literature_based_data
                    self.add_literature_based_data(entry, gene_node)

            if transporter_inhibitor_data:
                for entry in transporter_inhibitor_data:
                    self.add_transporter_inhibitor_node(entry)
            if stringdb_data and isinstance(stringdb_data, list):
                protein_name = f"{row.target}_protein"
                for entry in stringdb_data:
                    if entry.get("Ensembl", None):
                        self.add_ppi_data(protein_node, protein_name, entry)

        # Add metadata to the RDF graph
        self.add_metadata(metadata)

    def add_gene_protein_nodes(self, row):
        """Add gene and protein nodes based on the provided row data.

        :param row: Data for the gene/protein node.
        :return: Gene and protein nodes.
        """
        return add_gene_protein_nodes(self, row)

    def add_gene_disease_associations(self, id_number, source_idx, gene_node, disease, j):
        """Add gene-disease associations to the RDF graph.

        :param id_number: Unique identifier for the association.
        :param source_idx: Identifier for the source.
        :param gene_node: Node representing the gene.
        :param disease: Disease associated with the gene.
        :param j: Index of the disease in the list.
        """
        add_gene_disease_associations(
            self, id_number, source_idx, gene_node, disease, self.new_uris, j
        )

    def add_gene_expression_data(
        self, id_number, source_idx, gene_node, expression_data, experimental_process_data
    ):
        """Add gene expression data to the RDF graph.

        :param id_number: Unique identifier for the expression data.
        :param source_idx: Identifier for the source of the expression data.
        :param gene_node: Node representing the gene.
        :param expression_data: Expression data to be added.
        :param experimental_process_data: Experimental process data associated with the expression.
        """
        add_gene_expression_data(
            self,
            id_number,
            source_idx,
            gene_node,
            expression_data,
            experimental_process_data,
            self.new_uris,
        )

    def add_go_cpf(self, process_data):
        """Add Gene Ontology (GO) terms to the RDF graph.

        :param process_data: Process data related to GO and CPF.
        :return: Corresponding GO node.
        """
        return add_go_cpf(self, process_data)

    def add_compound_node(self, compound, protein_node):
        """Add compound data to the RDF graph and associate it with a protein node.

        :param compound: Compound data to be added.
        :param protein_node: Node representing the associated protein.
        """
        add_compound_node(self, compound, protein_node)

    def add_transporter_inhibitor_node(self, entry):
        """Add transporter inhibitor data to the RDF graph.

        :param entry: Data for the transporter inhibitor.
        """
        add_transporter_inhibitor_node(self, entry, self.base_uri)

    def add_pathway_node(self, data, source):
        """Add pathway data to the RDF graph.

        :param data: Data for the pathway.
        :param source: Source of the pathway data.
        :return: Corresponding pathway node.
        """
        return add_pathway_node(self, data, source)

    def add_literature_based_data(self, entry, gene_node):
        """Add literature-based data to the RDF graph.

        :param entry: Literature data to be added.
        :param gene_node: Node representing the gene.
        """
        add_literature_based_data(self, entry, gene_node)

    def add_ppi_data(self, protein_node, protein_name, entry):
        """Add Protein-Protein Interaction (PPI) data to the RDF graph.

        :param protein_node: Node representing the protein.
        :param protein_name: Name of the protein.
        :param entry: PPI data to be added.
        """
        add_ppi_data(self, protein_node, protein_name, entry, self.base_uri, self.new_uris)

    def add_metadata(self, metadata):
        """Add metadata to the RDF graph.

        :param metadata: Metadata to be added.
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
        path=None,
        threshold=0.001,
        uml_figure_path=None,
        print_string_output=True,
        additional_namespaces=None,
    ):
        """Run ShEx validation with optional parameters.

        :param path: Path to save the ShEx result.
        :param threshold: Validation threshold.
        :param uml_figure_path: Path to save UML diagram for validation.
        :param print_string_output: Whether to print the output string.
        :param additional_namespaces: Additional namespaces for validation.
        :return: ShEx validation result.
        """
        return get_shapes(
            self,
            self.base_uri,
            path,  # Only serialize if path is provided
            threshold,
            "shex",
            uml_figure_path,
            print_string_output,
            additional_namespaces,
        )

    def shacl(
        self,
        path=None,
        threshold=0.001,
        uml_figure_path=None,
        print_string_output=True,
        additional_namespaces=None,
    ):
        """Run SHACL validation with optional parameters.

        :param path: Path to save the SHACL result.
        :param threshold: Validation threshold.
        :param uml_figure_path: Path to save UML diagram for validation.
        :param print_string_output: Whether to print the output string.
        :param additional_namespaces: Additional namespaces for validation.
        :return: SHACL validation result.
        """
        return get_shapes(
            self,
            self.base_uri,
            path,  # Only serialize if path is provided
            threshold,
            "shacl",
            uml_figure_path,
            print_string_output,
            additional_namespaces,
        )

    def shacl_prefixes(self, path=None, namespaces=None):
        """Get SHACL prefixes with optional parameters.

        :param path: Path to save the SHACL prefixes.
        :param namespaces: Namespaces for the prefixes.
        :return: SHACL prefixes.
        """
        output_path = path if path is not None else self._prefixes_path
        current_namespaces = namespaces if namespaces is not None else self._namespaces

        return get_shacl_prefixes(
            namespaces=current_namespaces,
            path=output_path,
            new_uris=self.new_uris,  # Will only serialize if path is provided
        )
