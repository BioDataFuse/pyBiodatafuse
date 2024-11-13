# rdf.py


"""Provide function for the RDF module."""


from importlib import resources

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
from pyBiodatafuse.graph.rdf.utils import replace_na_none, matching_triples


def generate_rdf(
    df: pd.DataFrame,
    base_uri: str,
    version_iri: str,
    author: str,
    orcid: str,
    metadata: dict,
    open_only: bool,
    load_ontology: bool,
) -> Graph:
    """Generate an RDF graph from the provided DataFrame.

    :param df: DataFrame containing the data to be converted to RDF.
    :param base_uri: Base URI to use for the nodes in the graph.
    :param version_iri: Version IRI to use as the graph URI and add metadata.
    :param author: Author's name.
    :param orcid: Author's ORCID.
    :param metadata: Combined metadata for a BioDatafuse query.
    :param open_only: Boolean to include non-open data.
    :param load_ontology: Boolean to include annotations from the ontology to the graph.
    :return: RDF graph constructed from the DataFrame.
    """
    g = Graph()
    # Define the base URI and namespace bindings
    new_uris = {key: base_uri + value for key, value in URIS.items()}

    # Update g.bind only after constructing new URIs
    for key, new_value in new_uris.items():
        g.bind(key, new_value)

    for key, value in NAMESPACE_BINDINGS.items():
        g.bind(key, value)

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
        # transporter_inhibited_data = row.get(MOLMEDB_COMPOUND_PROTEIN_COL, None)
        stringdb_data = row.get(STRING_PPI_COL, None)
        # molmedb_data = row.get("MolMeDB_transporter_inhibitor", None)
        disease_data = []
        for source in [DISGENET_DISEASE_COL, OPENTARGETS_DISEASE_COL]:
            if (
                open_only and source == DISGENET_DISEASE_COL
            ):  # TODO implement open data only feature properly
                continue
            source_el = row.get(source)
            if isinstance(source_el, list):
                disease_data += source_el

        if (
            pd.isna(source_idx)
            or pd.isna(source_namespace)
            or pd.isna(target_idx)
            or pd.isna(target_namespace)
        ):
            continue

        source_curie = normalize_curie(f"{source_namespace}:{source_idx}")
        target_curie = normalize_curie(f"{target_namespace}:{target_idx}")

        if not source_curie or not target_curie:
            continue

        id_number = f"{i:06d}"
        gene_node, protein_node = add_gene_protein_nodes(g, row)

        if len(disease_data) > 0:
            for j, disease in enumerate(disease_data):
                add_gene_disease_associations(
                    g, id_number, source_idx, gene_node, disease, new_uris, j
                )
        if expression_data:
            add_gene_expression_data(
                g,
                id_number,
                source_idx,
                gene_node,
                expression_data,
                experimental_process_data,
                new_uris,
            )

        for source in ["WikiPathways", "MINERVA", "OpenTargets_reactome"]:
            if row.get(source, None):
                for pathway_data in row[source]:
                    if pathway_data["pathway_id"]:
                        pathway_node = add_pathway_node(g=g, data=pathway_data, source=source)
                        g.add((gene_node, URIRef(PREDICATES["sio_is_part_of"]), pathway_node))
                        g.add((protein_node, URIRef(PREDICATES["sio_is_part_of"]), pathway_node))
                        g.add((pathway_node, URIRef(PREDICATES["sio_has_part"]), gene_node))
                        g.add(
                            (
                                pathway_node,
                                URIRef(PREDICATES["sio_has_source"]),
                                URIRef(DATA_SOURCES[source]),
                            )
                        )

        if processes_data:
            for process_data in processes_data:
                go_cpf = add_go_cpf(g, process_data)
                if go_cpf:
                    g.add((gene_node, URIRef(PREDICATES["sio_is_part_of"]), go_cpf))
                    g.add((go_cpf, URIRef(PREDICATES["sio_has_part"]), gene_node))

        if compound_data:
            for compound in compound_data:
                add_compound_node(g, compound, protein_node)

        if literature_based_data:
            if isinstance(literature_based_data, list):
                entries = literature_based_data
                for entry in entries:
                    add_literature_based_data(g, entry, gene_node)
            elif isinstance(literature_based_data, dict):
                entry = literature_based_data
                add_literature_based_data(g, entry, gene_node)

        if transporter_inhibitor_data:
            for entry in transporter_inhibitor_data:
                add_transporter_inhibitor_node(g, entry, base_uri)
        if stringdb_data and isinstance(stringdb_data, list):
            protein_name = f"{row.target}_protein"
            protein_name = f"{row.target}_protein"
            for entry in stringdb_data:
                if entry.get("Ensembl", None):
                    add_ppi_data(g, protein_node, protein_name, entry, base_uri, new_uris)

    # Add metadata to the RDF graph
    if metadata:
        add_metadata(
            g=g,
            version_iri=version_iri,
            author=author,
            orcid=orcid,
            metadata=metadata,
            graph_uri=version_iri,
        )
        # Load ontology

    if load_ontology:
        with resources.files("pyBiodatafuse.resources").joinpath("biodatafuse.owl") as owl_path:
            temp_graph = Graph()
            temp_graph.parse(owl_path)

            # Preload subject, predicate, and object sets from `g` for quick lookups
            subj_set = {s for s, _, _ in g}
            pred_set = {p for _, p, _ in g}
            obj_set = {o for _, _, o in g}

            # Add only matching triples to `g`
            g.addN(matching_triples(temp_graph, subj_set, pred_set, obj_set))

    return g
