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
    - `process_ppi_data`: Processes Protein-Protein Interaction (PPI) data, adds it to RDF graph.
    - `_add_metadata`: Attaches metadata to the RDF graph.
    - `shex`: Runs shexer on the RDF graph to obtain its ShEx shapes.
    - `shacl`: Runs shexer on the RDF graph to obtain its SHACL shapes.
    - `shacl_prefixes`: Retrieves SHACL prefixes for the graph.
"""

import pandas as pd
from bioregistry import normalize_curie
from rdflib import Graph, URIRef
from tqdm import tqdm

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
from pyBiodatafuse.graph.rdf.nodes.gene import add_gene_nodes
from pyBiodatafuse.graph.rdf.nodes.gene_disease import add_gene_disease_associations
from pyBiodatafuse.graph.rdf.nodes.gene_expression import add_gene_expression_data
from pyBiodatafuse.graph.rdf.nodes.go_terms import add_go_cpf
from pyBiodatafuse.graph.rdf.nodes.literature import add_literature_based_data
from pyBiodatafuse.graph.rdf.nodes.pathway import add_pathway_node
from pyBiodatafuse.graph.rdf.nodes.protein_protein import add_ppi_data
from pyBiodatafuse.graph.rdf.utils import get_shacl_prefixes, get_shapes, replace_na_none


class BDFGraph(Graph):
    """Main class for a BioDatafuse RDF Graph, superclass of `rdflib.Graph`."""

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
        self.new_uris = {key: self.base_uri + value for key, value in URIS.items()}
        self._shex_path = None
        self._shacl_path = None
        self._prefixes_path = None
        self._namespaces = None
        self.include_variants = False  # TODO: Allow user to set options that can affect graph size

        # Bind prefixes
        for key, new_value in self.new_uris.items():
            self.bind(key, new_value)
        for key, value in NAMESPACE_BINDINGS.items():
            self.bind(key, value)

    def generate_rdf(self, df: pd.DataFrame, metadata: dict, open_only: bool = False):
        """
        Generate an RDF graph from the provided DataFrame and metadata.

        :param df: The DataFrame containing the data to be converted into RDF.
        :param metadata: A dictionary containing metadata information for RDF generation.
        :param open_only: A flag indicating whether to process only open data. Defaults to False.
        :param metadata: Metadata information to be added to the RDF graph.
        """
        df = df.applymap(replace_na_none)
        if not self.include_variants:
            df = df[df["target.source"] == "Ensembl"]
        for i, row in tqdm(df.iterrows(), total=df.shape[0], desc="Building RDF graph"):
            self.process_row(row, i, open_only)
        self._add_metadata(metadata)

    def process_row(self, row, i, open_only):
        """
        Process a single row of the DataFrame and update the RDF graph.

        :param row: A dictionary-like object representing a single row of the DataFrame.
        :param i: An integer representing the index of the row.
        :param open_only: A boolean indicating whether to process only open data.
        """
        source_idx = row.get(IDENTIFIER_COL)
        source_namespace = row.get(IDENTIFIER_SOURCE_COL)
        target_idx = row.get(TARGET_COL)
        target_namespace = row.get(TARGET_SOURCE_COL)
        if not self.valid_indices(source_idx, source_namespace, target_idx, target_namespace):
            return
        source_curie = normalize_curie(f"{source_namespace}:{source_idx}")
        target_curie = normalize_curie(f"{target_namespace}:{target_idx}")
        if not source_curie or not target_curie:
            return
        id_number = f"{i:06d}"
        gene_node = self.get_gene_node(row)
        if not gene_node:
            return
        disease_data = self.collect_disease_data(row)
        # New methods (e.g., new node types) can be called here
        self.process_ppi_data(row.get(STRING_PPI_COL), gene_node)
        protein_nodes = list(self.objects(gene_node, URIRef(PREDICATES["translation_of"])))
        self.process_disease_data(disease_data, id_number, source_idx, gene_node)
        self.process_expression_data(row, id_number, source_idx, gene_node)
        self.process_pathways(row, gene_node, protein_nodes)
        self.process_processes_data(row.get(OPENTARGETS_GO_COL), gene_node)
        self.process_compound_data(row.get(OPENTARGETS_GENE_COMPOUND_COL), gene_node)
        self.process_literature_data(
            row.get(LITERATURE_DISEASE_COL), gene_node, id_number, source_idx, self.new_uris, i
        )
        self.process_transporter_inhibitor_data(row.get(MOLMEDB_PROTEIN_COMPOUND_COL))
        if self.include_variants:
            self.process_protein_variants(protein_nodes)

    # Class methods about specific nodes begin here
    # If you add a new method, try to import most of the code from another script
    # Add new methods below
    def collect_disease_data(
        self,
        row,
    ):
        """
        Collect disease data from the row.

        :param row: A dictionary representing a row of data.
        :return: A list of collected disease data.
        """
        disease_data = []
        for source_col in [DISGENET_DISEASE_COL, OPENTARGETS_DISEASE_COL]:
            # if open_only and source_col == DISGENET_DISEASE_COL:
            #     continue  # TODO fix open data only feature
            source_data = row.get(source_col, None)
            if source_data:
                disease_data.extend(source_data)
        return disease_data

    def valid_indices(self, source_idx, source_namespace, target_idx, target_namespace):
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

    def get_gene_node(self, row):
        """
        Get gene node.

        :param row: A dictionary containing the data for a single row.
                    It must include the key "target.source".
        :return: A URIRef for the gene, else None.
        """
        if row["target.source"] == "Ensembl":
            return self._add_gene_nodes(row)
        return None

    def process_disease_data(self, disease_data, id_number, source_idx, gene_node):
        """
        Process disease data and add to the RDF graph.

        :param disease_data: List of disease data to be processed.
        :param id_number: Identifier number for the gene.
        :param source_idx: Source index for the data.
        :param gene_node: RDF node representing the gene.
        """
        for j, disease in enumerate(disease_data):
            self._add_gene_disease_associations(id_number, source_idx, gene_node, disease, j)

    def process_expression_data(self, row, id_number, source_idx, gene_node):
        """
        Process gene expression data and add to the RDF graph.

        :param row: The data row containing gene expression and experimental process data.
        :param id_number: The identifier number for the gene.
        :param source_idx: The source index for the data.
        :param gene_node: The RDF node representing the gene.
        """
        expression_data = row.get(BGEE_GENE_EXPRESSION_LEVELS_COL)
        experimental_process_data = row.get(PUBCHEM_COMPOUND_ASSAYS_COL)
        if expression_data:
            self._add_gene_expression_data(
                id_number,
                source_idx,
                gene_node,
                expression_data,
                experimental_process_data,
            )

    def process_pathways(self, row, gene_node, protein_nodes):
        """
        Process pathway data and add to the RDF graph.

        This method processes pathway data from various sources and adds the relevant
        information to the RDF graph. It creates pathway nodes and establishes relationships
        between gene nodes, protein nodes, and pathway nodes.

        :param row: A dictionary containing pathway data from different sources.
        :param gene_node: An RDF node representing the gene.
        :param protein_nodes: A list of RDF nodes representing proteins associated with the gene.
        """
        for source in ["WikiPathways", "MINERVA", "OpenTargets_reactome"]:
            pathway_data_list = row.get(source)
            if pathway_data_list:
                for pathway_data in pathway_data_list:
                    if pathway_data.get("pathway_id"):
                        pathway_node = self._add_pathway_node(pathway_data, source)
                        self.add((gene_node, URIRef(PREDICATES["sio_is_part_of"]), pathway_node))
                        self.add((pathway_node, URIRef(PREDICATES["sio_has_part"]), gene_node))
                        if protein_nodes:
                            for protein_node in protein_nodes:
                                self.add(
                                    (
                                        protein_node,
                                        URIRef(PREDICATES["sio_is_part_of"]),
                                        pathway_node,
                                    )
                                )
                                self.add(
                                    (pathway_node, URIRef(PREDICATES["sio_has_part"]), protein_node)
                                )
                        self.add((pathway_node, URIRef(PREDICATES["sio_has_part"]), gene_node))
                        self.add(
                            (
                                pathway_node,
                                URIRef(PREDICATES["sio_has_source"]),
                                URIRef(DATA_SOURCES[source]),
                            )
                        )

    def process_processes_data(self, processes_data, gene_node):
        """
        Process Gene Ontology (GO) terms and add to the RDF graph.

        :param processes_data: A list of GO terms related to a gene.
        :param gene_node: The RDF node representing the gene.
        """
        if processes_data:
            for process_data in processes_data:
                go_cpf = self._add_go_cpf(process_data)
                if go_cpf:
                    self.add((gene_node, URIRef(PREDICATES["sio_is_part_of"]), go_cpf))
                    self.add((go_cpf, URIRef(PREDICATES["sio_has_part"]), gene_node))

    def process_compound_data(self, compound_data, gene_node):
        """
        Process compound data and add to the RDF graph.

        :param compound_data: List of compounds to be processed.
        :param gene_node: URIRef of gene node.
        """
        if compound_data:
            for compound in compound_data:
                self._add_compound_node(compound, gene_node)

    def process_literature_data(
        self, literature_based_data, gene_node, id_number, source_idx, new_uris, i
    ):
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
                if entry.get("UMLS", None):
                    umls_parts = entry["UMLS"].split(":")
                    umlscui = umls_parts[1] if len(umls_parts) > 1 else None
                    if not umlscui:
                        continue
                    disease_data_lit = {
                        "disease_umlscui": umlscui,
                        "UMLS": umlscui,
                        "score": None,
                        "ei": None,
                        "el": None,
                        "disease_name": entry["disease_name"],
                    }
                    self._add_literature_based_data(
                        entry, gene_node, id_number, disease_data_lit, source_idx, new_uris, i
                    )

    def process_transporter_inhibitor_data(self, transporter_inhibitor_data):
        """
        Process transporter inhibitor data and add to the RDF graph.

        :param transporter_inhibitor_data: A list of transporter inhibitor data entries to be processed.
        """
        if transporter_inhibitor_data:
            for entry in transporter_inhibitor_data:
                self._add_transporter_inhibitor_node(entry)

    def process_protein_variants(self, protein_nodes):
        """
        Process protein variants and add to the RDF graph.

        This method iterates over a list of protein nodes and creates bidirectional
        "variant_of" relationships between each pair of protein nodes in the RDF graph.

        :param protein_nodes: A list of protein nodes to be processed.
        """
        if protein_nodes:
            for i, protein_node in enumerate(protein_nodes):
                for other_protein_node in protein_nodes[i + 1 :]:
                    self.add((protein_node, URIRef(PREDICATES["variant_of"]), other_protein_node))
                    self.add((other_protein_node, URIRef(PREDICATES["variant_of"]), protein_node))

    def process_ppi_data(self, stringdb_data, gene_node):
        """
        Process Protein-Protein Interaction (PPI) data and add to the RDF graph.

        :param stringdb_data: List of dictionaries containing PPI data from STRING database.
        :param gene_node: The gene URIRef.
        """
        if stringdb_data:
            for entry in stringdb_data:
                if entry.get("Ensembl"):
                    self._add_ppi_data(gene_node=gene_node, entry=entry)

    def _add_gene_nodes(self, row):
        """Add gene and protein nodes based on the provided row data.

        :param row: Data for the gene/protein node.
        :return: Gene and protein nodes.
        """
        return add_gene_nodes(self, row)

    def _add_gene_disease_associations(self, id_number, source_idx, gene_node, disease, j):
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

    def _add_gene_expression_data(
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

    def _add_go_cpf(self, process_data):
        """Add Gene Ontology (GO) terms to the RDF graph.

        :param process_data: Process data related to GO and CPF.
        :return: Corresponding GO node.
        """
        return add_go_cpf(self, process_data)

    def _add_compound_node(self, compound, protein_node):
        """Add compound data to the RDF graph and associate it with a protein node.

        :param compound: Compound data to be added.
        :param protein_node: Node representing the associated protein.
        """
        add_compound_node(self, compound, protein_node)

    def _add_transporter_inhibitor_node(self, entry):
        """Add transporter inhibitor data to the RDF graph.

        :param entry: Data for the transporter inhibitor.
        """
        add_transporter_inhibitor_node(self, entry, self.base_uri)

    def _add_pathway_node(self, data, source):
        """Add pathway data to the RDF graph.

        :param data: Data for the pathway.
        :param source: Source of the pathway data.
        :return: Corresponding pathway node.
        """
        return add_pathway_node(self, data, source)

    def _add_literature_based_data(
        self, entry, gene_node, id_number, disease_data, source_idx, new_uris, i
    ):
        """Add literature-based data to the RDF graph.

        :param entry: Literature data to be added.
        :param gene_node: Node representing the gene.
        :param id_number: Unique identifier for the expression data.
        :param disease_data: List of disease data to be processed.
        :param source_idx: Identifier for the source of the expression data.
        :param new_uris: Node URIs for the graph.
        :param i: An integer representing the index of the row.
        """
        add_literature_based_data(
            self, entry, gene_node, id_number, source_idx, disease_data, new_uris, i
        )

    def _add_ppi_data(self, gene_node, entry):
        """Add Protein-Protein Interaction (PPI) data to the RDF graph.

        :param gene_node: Node representing the gene.
        :param entry: PPI data to be added.
        """
        add_ppi_data(
            g=self, gene_node=gene_node, entry=entry, base_uri=self.base_uri, new_uris=self.new_uris
        )

    def _add_metadata(self, metadata):
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
        path=None,
        threshold=0.001,
        uml_figure_path=None,
        print_string_output=True,
        additional_namespaces=None,
    ):
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
        path=None,
        threshold=0.001,
        uml_figure_path=None,
        print_string_output=True,
        additional_namespaces=None,
    ):
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

    def shacl_prefixes(self, path=None, namespaces=None):
        """Get a SHACL prefixes graph, optionally add more namespaces to bind to it.

        :param path: Path to save the SHACL prefixes.
        :param namespaces: Namespaces for the prefixes.
        :return: SHACL prefixes.
        """
        output_path = path if path is not None else self._prefixes_path
        current_namespaces = self._namespaces
        if namespaces is not None:
            current_namespaces.update(namespaces)
        return get_shacl_prefixes(
            namespaces=current_namespaces,
            path=output_path,
            new_uris=self.new_uris,
        )


# Define new methods here
