# coding: utf-8

"""Python module to construct a NetworkX graph from the annotated data frame."""

import json
import logging
import os
import pickle
from logging import Logger
from typing import Any, Dict

import networkx as nx
import numpy as np
import pandas as pd
from tqdm import tqdm

import pyBiodatafuse.constants as Cons

logger = Logger(__name__)
logger.setLevel(logging.INFO)


def load_dataframe_from_pickle(pickle_path: str) -> pd.DataFrame:
    """Load a previously annotated DataFrame from a pickle file.

    :param pickle_path: the path to a previously obtained annotation DataFrame dumped as a pickle file.
    :returns: a Pandas DataFrame.
    """
    with open(pickle_path, "rb") as rin:
        df = pickle.load(rin)

    return df


def merge_node(g, node_label, node_attrs):
    """Merge the attr_dict of a newly added node to the graph on duplication, otherwise, add the new node.

    :param g: the graph to which the node will be added.
    :param node_label: node label.
    :param node_attrs: dictionary of node attributes.
    """
    if node_label not in g.nodes():
        # Ensure 'labels' is set
        if Cons.LABEL not in node_attrs:
            node_attrs[Cons.LABEL] = node_attrs.get("label", "Unknown")
        g.add_node(node_label, attr_dict=node_attrs)
    else:
        if "attr_dict" in g.nodes[node_label]:
            for k, v in node_attrs.items():
                if k in g.nodes[node_label]["attr_dict"]:
                    if g.nodes[node_label]["attr_dict"][k] is not None:
                        if isinstance(v, str):
                            v_list = g.nodes[node_label]["attr_dict"][k].split("|")
                            v_list.append(v)
                            g.nodes[node_label]["attr_dict"][k] = "|".join(list(set(v_list)))
                    else:
                        g.nodes[node_label]["attr_dict"][k] = v
                else:
                    g.nodes[node_label]["attr_dict"][k] = v
        else:
            if Cons.LABEL not in node_attrs:
                node_attrs[Cons.LABEL] = node_attrs.get("label", "Unknown")
            g.add_node(node_label, attr_dict=node_attrs)


"""Adding node and edges from annotators"""


def add_gene_bgee_subgraph(g, gene_node_label, annot_list):
    """Construct part of the graph by linking the gene to a list of anatomical entities.

    :param g: the input graph to extend with new nodes and edges.
    :param gene_node_label: the gene node to be linked to annotation entities.
    :param annot_list: list of anatomical entities from Bgee with gene expression levels.
    :returns: a NetworkX MultiDiGraph
    """
    logger.debug("Adding Bgee nodes and edges")
    for annot in annot_list:
        if pd.isna(annot[Cons.ANATOMICAL_NAME]):
            continue

        annot_node_label = annot[Cons.BGEE_ANATOMICAL_NODE_MAIN_LABEL].replace(":", "_")
        entity_attrs = Cons.BGEE_ANATOMICAL_NODE_ATTRS.copy()
        entity_attrs.update(
            {
                Cons.NAME: annot[Cons.ANATOMICAL_NAME],
                Cons.ID: annot[Cons.ANATOMICAL_ID],
                Cons.DATASOURCE: Cons.BGEE,
                Cons.UBERON: annot[Cons.ANATOMICAL_ID].split(":")[1],
            }
        )

        g.add_node(annot_node_label, attr_dict=entity_attrs)

        edge_attrs = Cons.BGEE_EDGE_ATTRS.copy()
        fields = {
            Cons.CONFIDENCE_ID,
            Cons.CONFIDENCE_LEVEL_NAME,
            Cons.EXPRESSION_LEVEL,
            Cons.DEVELOPMENTAL_ID,
            Cons.DEVELOPMENTAL_STAGE_NAME,
        }

        for field in fields:
            if pd.notna(annot[field]):
                edge_attrs[field] = annot[field]

        edge_hash = hash(frozenset(edge_attrs.items()))
        edge_attrs[Cons.EDGE_HASH] = edge_hash
        edge_data = g.get_edge_data(gene_node_label, annot_node_label)
        edge_data = {} if edge_data is None else edge_data
        node_exists = [
            x for x, y in edge_data.items() if y["attr_dict"][Cons.EDGE_HASH] == edge_hash
        ]

        if len(node_exists) == 0:
            g.add_edge(
                gene_node_label,
                annot_node_label,
                label=Cons.BGEE_GENE_ANATOMICAL_EDGE_LABEL,
                attr_dict=edge_attrs,
            )

    return g


def add_disgenet_gene_disease_subgraph(g, gene_node_label, annot_list):
    """Construct part of the graph by linking the gene to diseases.

    :param g: the input graph to extend with new nodes and edges.
    :param gene_node_label: the gene node to be linked to diseases.
    :param annot_list: list of diseases from DisGeNET.
    :returns: a NetworkX MultiDiGraph
    """
    logger.debug("Adding DisGeNET nodes and edges")
    for annot in annot_list:
        if pd.isna(annot[Cons.DISEASE_NAME]):
            continue

        annot_node_label = annot[Cons.DISEASE_NODE_MAIN_LABEL]
        annot_node_attrs = Cons.DISGENET_DISEASE_NODE_ATTRS.copy()
        annot_node_attrs.update(
            {
                Cons.NAME: annot[Cons.DISEASE_NAME],
                Cons.ID: annot[Cons.UMLS],
                Cons.DATASOURCE: Cons.DISGENET,
            }
        )

        other_ids = {
            Cons.HPO: annot[Cons.HPO],
            Cons.NCI: annot[Cons.NCI],
            Cons.OMIM: annot[Cons.OMIM],
            Cons.MONDO: annot[Cons.MONDO],
            Cons.ORDO: annot[Cons.ORDO],
            Cons.EFO: annot[Cons.EFO],
            Cons.DO: annot[Cons.DO],
            Cons.MESH: annot[Cons.MESH],
            Cons.UMLS: annot[Cons.UMLS],
            Cons.DISEASE_TYPE: annot[Cons.DISEASE_TYPE],
        }

        for key, value in other_ids.items():
            if pd.notna(value):
                annot_node_attrs[key] = value

        g.add_node(annot_node_label, attr_dict=annot_node_attrs)

        edge_attrs = Cons.DISGENET_EDGE_ATTRS.copy()
        edge_attrs[Cons.DISGENET_SCORE] = annot[Cons.DISGENET_SCORE]

        if pd.notna(annot[Cons.DISGENET_EI]):
            edge_attrs[Cons.DISGENET_EI] = annot[Cons.DISGENET_EI]
        if pd.notna(annot[Cons.DISGENET_EL]):
            edge_attrs[Cons.DISGENET_EL] = annot[Cons.DISGENET_EL]

        edge_hash = hash(frozenset(edge_attrs.items()))
        edge_attrs[Cons.EDGE_HASH] = edge_hash  # type: ignore
        edge_data = g.get_edge_data(gene_node_label, annot_node_label)
        edge_data = {} if edge_data is None else edge_data
        node_exists = [
            x for x, y in edge_data.items() if y["attr_dict"][Cons.EDGE_HASH] == edge_hash
        ]

        if len(node_exists) == 0:
            g.add_edge(
                gene_node_label,
                annot_node_label,
                label=Cons.GENE_DISEASE_EDGE_LABEL,
                attr_dict=edge_attrs,
            )

    return g


def add_intact_interactions_subgraph(g, gene_node_label, annot_list):
    """Construct part of the graph by linking the gene interactions via IntAct, including all interaction attributes.

    :param g: the input graph to extend with new nodes and edges.
    :param gene_node_label: the gene node to be linked to compounds or proteins.
    :param annot_list: list of interactions from IntAct.
    :returns: a NetworkX MultiDiGraph
    """
    logger.debug("Adding IntAct nodes and edges")
    if not hasattr(add_intact_interactions_subgraph, "seen_interaction_ids"):
        add_intact_interactions_subgraph.seen_interaction_ids = set()

    seen_interaction_ids = add_intact_interactions_subgraph.seen_interaction_ids
    edges_seen = {}

    for interaction in annot_list:
        interaction_id = interaction[Cons.INTACT_INTERACTION_ID]
        if not interaction_id or interaction_id in seen_interaction_ids:
            continue

        seen_interaction_ids.add(interaction_id)

        try:
            id_a = interaction[Cons.INTACT_ID_A]
            id_b = interaction[Cons.INTACT_ID_B]
        except KeyError:
            continue

        is_a_chebi = isinstance(id_a, str) and id_a.startswith("CHEBI:")
        is_b_chebi = isinstance(id_b, str) and id_b.startswith("CHEBI:")

        if is_a_chebi:
            partner_node_label = id_a
            partner_name = interaction[Cons.INTACT_INTERACTOR_A_NAME]
            partner_species = interaction[Cons.INTACT_INTERACTOR_A_SPECIES]
            molecule = interaction[Cons.INTACT_MOLECULE_A]
            is_compound = True
        elif is_b_chebi:
            partner_node_label = id_b
            partner_name = interaction[Cons.INTACT_INTERACTOR_B_NAME]
            partner_species = interaction[Cons.INTACT_INTERACTOR_B_SPECIES]
            molecule = interaction[Cons.INTACT_MOLECULE_B]
            is_compound = True
        else:
            partner_node_label = interaction[Cons.INTACT_PPI_EDGE_MAIN_LABEL]
            is_compound = False

        if not partner_node_label or pd.isna(partner_node_label):
            continue

        edge_key = (gene_node_label, partner_node_label)
        if edge_key in edges_seen:
            existing = edges_seen[edge_key]
            method = interaction[Cons.INTACT_DETECTION_METHOD]
            if method:
                if isinstance(existing[Cons.INTACT_DETECTION_METHOD], list):
                    existing[Cons.INTACT_DETECTION_METHOD].append(method)
                else:
                    existing[Cons.INTACT_DETECTION_METHOD] = [
                        existing[Cons.INTACT_DETECTION_METHOD],
                        method,
                    ]
            continue

        if is_compound:
            annot_node_attrs = Cons.INTACT_COMPOUND_NODE_ATTRS.copy()
            annot_node_attrs[Cons.ID] = partner_node_label
            annot_node_attrs[Cons.NAME] = partner_name
            annot_node_attrs[Cons.SPECIES] = partner_species
            annot_node_attrs[Cons.MOLECULE] = molecule
            merge_node(g, partner_node_label, annot_node_attrs)

        edge_attrs = Cons.INTACT_PPI_EDGE_ATTRS.copy()

        for key, value in interaction.items():
            if value is None or (isinstance(value, str) and not value.strip()):
                continue
            edge_attrs[key] = ",".join(map(str, value)) if isinstance(value, list) else value

        edge_attrs[Cons.EDGE_HASH] = hash(frozenset(edge_attrs.items()))

        edges_seen[edge_key] = edge_attrs

    for (source, target), edge_attrs in edges_seen.items():
        g.add_edge(
            source,
            target,
            label=edge_attrs[Cons.LABEL],
            attr_dict=edge_attrs,
        )

    return g


# TODO: test this function
def add_intact_compound_interactions_subgraph(g, compound_node_label, annot_list):
    """Construct part of the graph by linking the compound interactions via IntAct, including all interaction attributes.

    :param g: the input graph to extend with new nodes and edges.
    :param compound_node_label: the compound node label (used as source node), expected as a ChEBI ID (e.g., '15361').
    :param annot_list: list of interaction dicts from IntAct.
    :returns: a NetworkX MultiDiGraph
    """
    logger.debug("Adding IntAct compound nodes and edges")
    if not hasattr(add_intact_compound_interactions_subgraph, "seen_interaction_ids"):
        add_intact_compound_interactions_subgraph.seen_interaction_ids = set()

    seen_interaction_ids = add_intact_compound_interactions_subgraph.seen_interaction_ids
    edges_seen = {}

    compound_full_id = f"CHEBI:{compound_node_label.strip()}"

    for interaction in annot_list:
        interaction_id = interaction[Cons.INTACT_BINARY_INTERACTION_ID]
        if not interaction_id or interaction_id in seen_interaction_ids:
            continue
        seen_interaction_ids.add(interaction_id)

        id_a = interaction[Cons.INTACT_ID_A].strip()
        id_b = interaction[Cons.INTACT_ID_B].strip()

        if compound_full_id == id_a:
            partner_id = id_b
            partner_name = interaction[Cons.INTACT_INTERACTOR_B_NAME]
            partner_species = interaction[Cons.INTACT_INTERACTOR_B_SPECIES]
            partner_molecule = interaction[Cons.INTACT_MOLECULE_B]
        elif compound_full_id == id_b:
            partner_id = id_a
            partner_name = interaction[Cons.INTACT_INTERACTOR_A_NAME]
            partner_species = interaction[Cons.INTACT_INTERACTOR_A_SPECIES]
            partner_molecule = interaction[Cons.INTACT_MOLECULE_A]
        else:
            continue

        if not partner_id or pd.isna(partner_id):
            continue

        edge_key = (compound_node_label, partner_id)
        if edge_key in edges_seen:
            existing_methods = edges_seen[edge_key].get("detection_method", "")
            new_method = interaction.get("detection_method", "")
            if new_method and new_method not in existing_methods:
                existing_methods += f",{new_method}"
                edges_seen[edge_key]["detection_method"] = existing_methods
            continue

        node_label_type = (
            Cons.COMPOUND_NODE_LABELS if partner_id.startswith("CHEBI:") else Cons.GENE_NODE_LABELS
        )

        partner_node_attrs = {
            "id": partner_id,
            "label": partner_name,
            "species": partner_species,
            "molecule": partner_molecule,
            Cons.LABEL: node_label_type,
        }
        merge_node(g, partner_id, partner_node_attrs)

        edge_attrs = {
            key: ",".join(map(str, value)) if isinstance(value, list) else value
            for key, value in interaction.items()
            if value is not None and not (isinstance(value, str) and not value.strip())
        }
        edge_attrs["detection_method"] = interaction.get("detection_method", "")
        edge_attrs["interaction_type"] = interaction.get("type", "compound-ppi")
        edge_attrs["edge_hash"] = hash(frozenset(edge_attrs.items()))

        edges_seen[edge_key] = edge_attrs

    for (source, target), edge_attrs in edges_seen.items():
        g.add_edge(
            source,
            target,
            label=edge_attrs.get("interaction_type", "compound-ppi"),
            attr_dict=edge_attrs,
        )

    return g


def add_literature_gene_disease_subgraph(g, gene_node_label, annot_list):
    """Construct part of the graph by linking the gene to diseases form literature.

    :param g: the input graph to extend with new nodes and edges.
    :param gene_node_label: the gene node to be linked to diseases.
    :param annot_list: list of diseases from DisGeNET.
    :returns: a NetworkX MultiDiGraph
    """
    logger.debug("Adding literature disease nodes and edges")
    for annot in annot_list:
        if pd.isna(annot["disease_name"]):
            continue
        annot_node_label = annot[Cons.LITERATURE_NODE_MAIN_LABEL]
        annot_node_attrs = Cons.LITERATURE_DISEASE_NODE_ATTRS.copy()
        annot_node_attrs.update(
            {"datasource": annot["source"], "name": annot["disease_name"], "id": annot[Cons.UMLS]}
        )

        other_ids = {
            Cons.UMLS: annot[Cons.UMLS],
            Cons.MONDO: annot[Cons.MONDO],
        }

        for key, value in other_ids.items():
            if pd.notna(value):
                annot_node_attrs[key] = value

        g.add_node(annot_node_label, attr_dict=annot_node_attrs)

        edge_attrs = Cons.LITERATURE_DISEASE_EDGE_ATTRS.copy()
        edge_attrs["datasource"] = annot["source"]

        edge_hash = hash(frozenset(edge_attrs.items()))
        edge_attrs["edge_hash"] = edge_hash
        edge_data = g.get_edge_data(gene_node_label, annot_node_label)
        edge_data = {} if edge_data is None else edge_data
        node_exists = [x for x, y in edge_data.items() if y["attr_dict"]["edge_hash"] == edge_hash]

        if len(node_exists) == 0:
            g.add_edge(
                gene_node_label,
                annot_node_label,
                label=Cons.GENE_DISEASE_EDGE_LABEL,
                attr_dict=edge_attrs,
            )

    return g


def add_minerva_gene_pathway_subgraph(g, gene_node_label, annot_list):
    """Construct part of the graph by linking the gene to MINERVA pathways.

    :param g: the input graph to extend with new nodes and edges.
    :param gene_node_label: the gene node to be linked to MINERVA pathways.
    :param annot_list: list of MINERVA pathways from MINERVA.
    :returns: a NetworkX MultiDiGraph
    """
    logger.debug("Adding MINERVA nodes and edges")
    for annot in annot_list:
        if pd.isna(annot[Cons.PATHWAY_LABEL]):
            continue

        annot_node_label = annot[Cons.MINERVA_PATHWAY_NODE_MAIN_LABEL]
        annot_node_attrs = Cons.MINERVA_PATHWAY_NODE_ATTRS.copy()
        annot_node_attrs.update(
            {
                Cons.DATASOURCE: Cons.MINERVA,
                Cons.NAME: annot[Cons.PATHWAY_LABEL],
                Cons.ID: annot[Cons.PATHWAY_ID],
                Cons.GENE_COUNTS: annot[Cons.PATHWAY_GENE_COUNTS],
            }
        )

        g.add_node(annot_node_label, attr_dict=annot_node_attrs)

        edge_attrs = Cons.GENE_PATHWAY_EDGE_ATTRS.copy()
        edge_attrs[Cons.DATASOURCE] = Cons.MINERVA

        edge_hash = hash(frozenset(edge_attrs.items()))
        edge_attrs[Cons.EDGE_HASH] = edge_hash
        edge_data = g.get_edge_data(gene_node_label, annot_node_label)
        edge_data = {} if edge_data is None else edge_data
        node_exists = [
            x for x, y in edge_data.items() if y["attr_dict"][Cons.EDGE_HASH] == edge_hash
        ]

        if len(node_exists) == 0:
            g.add_edge(
                gene_node_label,
                annot_node_label,
                label=Cons.GENE_PATHWAY_EDGE_LABEL,
                attr_dict=edge_attrs,
            )

    return g


def add_wikipathways_gene_pathway_subgraph(g, gene_node_label, annot_list):
    """Construct part of the graph by linking the gene to pathways from WikiPathways.

    :param g: the input graph to extend with new nodes and edges.
    :param gene_node_label: the gene node to be linked to pathways from WikiPathways.
    :param annot_list: list of pathways from WikiPathways.
    :returns: a NetworkX MultiDiGraph
    """
    logger.debug("Adding WikiPathways pathway nodes and edges")
    for annot in annot_list:
        if pd.isna(annot[Cons.WIKIPATHWAYS_NODE_MAIN_LABEL]):
            continue

        annot_node_label = annot[Cons.WIKIPATHWAYS_NODE_MAIN_LABEL]
        annot_node_attrs = Cons.WIKIPATHWAYS_NODE_ATTRS.copy()
        annot_node_attrs.update(
            {
                Cons.DATASOURCE: Cons.WIKIPATHWAYS,
                Cons.NAME: annot[Cons.PATHWAY_LABEL],
                Cons.ID: annot[Cons.PATHWAY_ID],
                Cons.GENE_COUNTS: annot[Cons.PATHWAY_GENE_COUNTS],
            }
        )

        g.add_node(annot_node_label, attr_dict=annot_node_attrs)

        edge_attrs = Cons.GENE_PATHWAY_EDGE_ATTRS.copy()
        edge_attrs[Cons.DATASOURCE] = Cons.WIKIPATHWAYS

        edge_hash = hash(frozenset(edge_attrs.items()))
        edge_attrs[Cons.EDGE_HASH] = edge_hash
        edge_data = g.get_edge_data(gene_node_label, annot_node_label)
        edge_data = {} if edge_data is None else edge_data
        node_exists = [
            x for x, y in edge_data.items() if y["attr_dict"][Cons.EDGE_HASH] == edge_hash
        ]

        if len(node_exists) == 0:
            g.add_edge(
                gene_node_label,
                annot_node_label,
                label=Cons.GENE_PATHWAY_EDGE_LABEL,
                attr_dict=edge_attrs,
            )

    return g


def add_kegg_gene_pathway_subgraph(g, gene_node_label, annot_list):
    """Construct part of the graph by linking the gene to pathways from KEGG.

    :param g: the input graph to extend with new nodes and edges.
    :param gene_node_label: the gene node to be linked to pathways from KEGG.
    :param annot_list: list of pathways from KEGG.
    :returns: a NetworkX MultiDiGraph
    """
    logger.debug("Adding KEGG nodes and edges")
    for annot in annot_list:
        if pd.isna(annot[Cons.PATHWAY_LABEL]):
            continue

        annot_node_label = annot[Cons.KEGG_PATHWAY_NODE_MAIN_LABEL]
        annot_node_attrs = Cons.KEGG_PATHWAY_NODE_ATTRS.copy()
        annot_node_attrs.update(
            {
                Cons.DATASOURCE: Cons.KEGG,
                Cons.NAME: annot[Cons.PATHWAY_LABEL],
                Cons.ID: annot[Cons.PATHWAY_ID],
                Cons.GENE_COUNTS: annot[Cons.PATHWAY_GENE_COUNTS],
            }
        )

        # g.add_node(annot_node_label, attr_dict=annot_node_attrs)
        merge_node(g, annot_node_label, annot_node_attrs)

        edge_attrs = Cons.GENE_PATHWAY_EDGE_ATTRS.copy()
        edge_attrs[Cons.DATASOURCE] = Cons.KEGG

        edge_hash = hash(frozenset(edge_attrs.items()))
        edge_attrs[Cons.EDGE_HASH] = edge_hash
        edge_data = g.get_edge_data(gene_node_label, annot_node_label)
        edge_data = {} if edge_data is None else edge_data
        node_exists = [
            x for x, y in edge_data.items() if y["attr_dict"][Cons.EDGE_HASH] == edge_hash
        ]

        if len(node_exists) == 0:
            g.add_edge(
                gene_node_label,
                annot_node_label,
                label=Cons.GENE_PATHWAY_EDGE_LABEL,
                attr_dict=edge_attrs,
            )

    return g


def add_kegg_compounds_subgraph(g, pathway_node_label, compounds_list, combined_df):
    """Construct part of the graph by linking the KEGG compound to its respective pathway.

    :param g: the input graph to extend with new nodes and edges.
    :param pathway_node_label: the pathway node to be linked to compound nodes.
    :param compounds_list: list of compounds from KEGG.
    :param combined_df: the combined dataframe.
    :returns: a NetworkX MultiDiGraph
    """
    logger.debug("Adding KEGG compound nodes and edges")
    for compound in compounds_list:
        if pd.isna(compound[Cons.KEGG_COMPOUND_NAME]):
            continue

        annot_node_label = compound[Cons.KEGG_IDENTIFIER]
        annot_node_attrs = Cons.KEGG_COMPOUND_NODE_ATTRS.copy()
        annot_node_attrs.update(
            {
                Cons.ID: compound[Cons.KEGG_IDENTIFIER],
                Cons.LABEL: compound[Cons.KEGG_COMPOUND_NAME],
            }
        )

        merge_node(g, annot_node_label, annot_node_attrs)

        for _, path_row in combined_df.iterrows():
            pathways = path_row.get(Cons.PATHWAYS, [])
            if isinstance(pathways, list) and pathways:
                for pathway in pathways:
                    if pathway_node_label != pathway.get(Cons.PATHWAYS, ""):
                        continue

                    if Cons.PATHWAY_COMPOUNDS in pathway:
                        pathway_compounds = [
                            comp[Cons.KEGG_IDENTIFIER] for comp in pathway[Cons.PATHWAY_COMPOUNDS]
                        ]
                        if compound[Cons.KEGG_IDENTIFIER] in pathway_compounds:
                            edge_attrs = Cons.KEGG_COMPOUND_EDGE_ATTRS.copy()
                            edge_hash = hash(frozenset(edge_attrs.items()))
                            edge_attrs[Cons.EDGE_HASH] = edge_hash
                            edge_data = g.get_edge_data(pathway_node_label, annot_node_label)
                            edge_data = {} if edge_data is None else edge_data
                            node_exists = [
                                x
                                for x, y in edge_data.items()
                                if "attr_dict" in y
                                and y["attr_dict"].get(Cons.EDGE_HASH) == edge_hash
                            ]

            if not isinstance(pathways, list):
                continue

            for pathway in pathways:
                if pathway_node_label != pathway.get(Cons.PATHWAYS):
                    continue

                if Cons.PATHWAY_COMPOUNDS not in pathway:
                    continue

                pathway_compounds = [
                    comp[Cons.KEGG_IDENTIFIER] for comp in pathway[Cons.PATHWAY_COMPOUNDS]
                ]
                if compound[Cons.KEGG_IDENTIFIER] not in pathway_compounds:
                    continue

                edge_attrs = Cons.KEGG_COMPOUND_EDGE_ATTRS.copy()
                edge_hash = hash(frozenset(edge_attrs.items()))
                edge_attrs[Cons.EDGE_HASH] = edge_hash  # type: ignore
                edge_data = g.get_edge_data(pathway_node_label, annot_node_label)
                edge_data = {} if edge_data is None else edge_data
                node_exists = [
                    x
                    for x, y in edge_data.items()
                    if "attr_dict" in y and y["attr_dict"].get(Cons.EDGE_HASH) == edge_hash
                ]

                if len(node_exists) == 0:
                    g.add_edge(
                        pathway_node_label,
                        annot_node_label,
                        label=Cons.KEGG_COMPOUND_EDGE_LABEL,
                        attr_dict=edge_attrs,
                    )

    return g


def process_kegg_pathway_compound(g, kegg_pathway_compound, combined_df):
    """Process pathway-compound relationships from KEGG and add them to the graph.

    :param g: the input graph to extend with pathway-compound relationships.
    :param kegg_pathway_compound: DataFrame containing pathway-compound relationships.
    :param combined_df: DataFrame containing KEGG pathway data.
    """
    logger.debug("Processing KEGG pathway-compound relationships")
    for _, row in kegg_pathway_compound.iterrows():
        compound_info = row[Cons.KEGG_PATHWAY_COL]

        if isinstance(compound_info, dict):
            compounds_list = [compound_info]
        elif isinstance(compound_info, list):
            compounds_list = compound_info
        else:
            compounds_list = []

        for compound in compounds_list:
            compound_id = compound[Cons.KEGG_IDENTIFIER]

            for _, pathway_row in combined_df.iterrows():
                pathway_data = pathway_row[Cons.PATHWAY_COMPOUNDS]

                if not isinstance(pathway_data, list):
                    continue

                for pathway in pathway_data:
                    pathway_id = pathway.get(Cons.PATHWAY_ID)
                    pathway_compounds = pathway.get(Cons.PATHWAY_COMPOUNDS, [])

                    if any(c.get(Cons.KEGG_IDENTIFIER) == compound_id for c in pathway_compounds):
                        add_kegg_compounds_subgraph(g, pathway_id, compounds_list, combined_df)


def _add_grofilier_genes_node_attrs(annot):
    annot_node_attrs = Cons.GPROFILER_GENE_NODE_ATTRS.copy()
    annot_node_attrs.update(
        {
            Cons.NAME: annot[Cons.GPROFILER_NAME],
            Cons.ID: annot[Cons.GPROFILER_ID],
            Cons.GENE_COUNTS: annot[Cons.GPROFILER_TERM_SIZE],
            Cons.P_VALUE: annot[Cons.P_VALUE],
            Cons.SIGNIFICANT: annot[Cons.SIGNIFICANT],
        }
    )

    return annot_node_attrs


def add_gprofiler_gene_phenotype_subgraph(g, gene_node_label, annot_list):
    """Construct part of the graph by linking the gene to phenotypes from GProfiler.

    :param g: the input graph to extend with new nodes and edges.
    :param gene_node_label: the gene node to be linked to phenotypes from GProfiler.
    :param annot_list: list of phenotypes from GProfiler.
    :returns: a NetworkX MultiDiGraph
    """
    for annot in annot_list:
        if pd.isna(annot["name"]):
            continue

        annot_node_label = annot[Cons.GPROFILER_ID]
        annot_node_attrs = _add_grofilier_genes_node_attrs(annot)
        annot_node_attrs[Cons.LABEL] = Cons.PHENOTYPE_NODE_LABEL

        g.add_node(annot_node_label, attr_dict=annot_node_attrs)

        edge_attrs = Cons.GPROFILER_EDGE_ATTRS.copy()
        edge_attrs[Cons.DATASOURCE] = annot[Cons.DATASOURCE]
        edge_hash = hash(frozenset(edge_attrs.items()))
        edge_attrs[Cons.EDGE_HASH] = edge_hash
        edge_data = g.get_edge_data(gene_node_label, annot_node_label)
        edge_data = {} if edge_data is None else edge_data
        node_exists = [
            x for x, y in edge_data.items() if y["attr_dict"][Cons.EDGE_HASH] == edge_hash
        ]
        if len(node_exists) == 0:
            g.add_edge(
                gene_node_label,
                annot_node_label,
                label=Cons.GPROFILE_GENE_HP_EDGE_LABEL,
                attr_dict=edge_attrs,
            )

    return g


def add_gprofiler_gene_hpa_subgraph(g, gene_node_label, annot_list):
    """Construct part of the graph by linking the gene to anatomical entities from GProfiler.

    :param g: the input graph to extend with new nodes and edges.
    :param gene_node_label: the gene node to be linked to anatomical entities from GProfiler.
    :param annot_list: list of pathways from WikiPathways.
    :returns: a NetworkX MultiDiGraph
    """
    for annot in annot_list:
        if pd.isna(annot[Cons.NAME]):
            continue

        annot_node_label = annot[Cons.GPROFILER_ID]
        annot_node_attrs = _add_grofilier_genes_node_attrs(annot)
        annot_node_attrs[Cons.LABEL] = Cons.ANATOMICAL_NODE_LABEL

        g.add_node(annot_node_label, attr_dict=annot_node_attrs)

        edge_attrs = Cons.GPROFILER_EDGE_ATTRS.copy()
        edge_attrs[Cons] = annot[Cons.DATASOURCE]

        edge_hash = hash(frozenset(edge_attrs.items()))
        edge_attrs[Cons.EDGE_HASH] = edge_hash
        edge_data = g.get_edge_data(gene_node_label, annot_node_label)
        edge_data = {} if edge_data is None else edge_data
        node_exists = [
            x for x, y in edge_data.items() if y["attr_dict"][Cons.EDGE_HASH] == edge_hash
        ]
        if len(node_exists) == 0:
            g.add_edge(
                gene_node_label,
                annot_node_label,
                label=Cons.GPROFILE_GENE_HPA_EDGE_LABEL,
                attr_dict=edge_attrs,
            )

    return g


def add_gprofiler_gene_kegg_subgraph(g, gene_node_label, annot_list):
    """Construct part of the graph by linking the gene to KEGG pathways from GProfiler.

    :param g: the input graph to extend with new nodes and edges.
    :param gene_node_label: the gene node to be linked to KEGG pathways from GProfiler.
    :param annot_list: list of KEGG pathways from GProfiler.
    :returns: a NetworkX MultiDiGraph
    """
    for annot in annot_list:
        if pd.isna(annot[Cons.NAME]):
            continue

        annot_node_label = annot[Cons.GPROFILER_ID]
        annot_node_attrs = _add_grofilier_genes_node_attrs(annot)
        annot_node_attrs[Cons.DATASOURCE] = annot[Cons.DATASOURCE]
        annot_node_attrs[Cons.GPROFILER_PATHWAY_TYPE] = Cons.KEGG
        annot_node_attrs[Cons.LABEL] = Cons.ANATOMICAL_NODE_LABEL

        g.add_node(annot_node_label, attr_dict=annot_node_attrs)

        edge_attrs = Cons.GENE_PATHWAY_EDGE_ATTRS.copy()
        edge_attrs[Cons.DATASOURCE] = annot[Cons.DATASOURCE]

        edge_hash = hash(frozenset(edge_attrs.items()))
        edge_attrs[Cons.EDGE_HASH] = edge_hash
        edge_data = g.get_edge_data(gene_node_label, annot_node_label)
        edge_data = {} if edge_data is None else edge_data
        node_exists = [
            x for x, y in edge_data.items() if y["attr_dict"][Cons.EDGE_HASH] == edge_hash
        ]
        if len(node_exists) == 0:
            g.add_edge(
                gene_node_label,
                annot_node_label,
                label=Cons.GENE_PATHWAY_EDGE_LABEL,
                attr_dict=edge_attrs,
            )

    return g


def add_gprofiler_gene_mirna_subgraph(g, gene_node_label, annot_list):
    """Construct part of the graph by linking the gene to miRNAs from GProfiler.

    :param g: the input graph to extend with new nodes and edges.
    :param gene_node_label: the gene node to be linked to miRNAs from GProfiler.
    :param annot_list: list of miRNAs from GProfiler.
    :returns: a NetworkX MultiDiGraph
    """
    for annot in annot_list:
        if pd.isna(annot[Cons.NAME]):
            continue

        annot_node_label = annot[Cons.GPROFILER_ID]
        annot_node_attrs = _add_grofilier_genes_node_attrs(annot)
        annot_node_attrs[Cons.DATASOURCE] = annot[Cons.DATASOURCE]
        annot_node_attrs[Cons.LABEL] = Cons.MIRNA_NODE_LABEL

        g.add_node(annot_node_label, attr_dict=annot_node_attrs)

        edge_attrs = Cons.GPROFILER_EDGE_ATTRS.copy()
        edge_attrs[Cons.DATASOURCE] = annot[Cons.DATASOURCE]

        edge_hash = hash(frozenset(edge_attrs.items()))
        edge_attrs[Cons.EDGE_HASH] = edge_hash
        edge_data = g.get_edge_data(gene_node_label, annot_node_label)
        edge_data = {} if edge_data is None else edge_data
        node_exists = [
            x for x, y in edge_data.items() if y["attr_dict"][Cons.EDGE_HASH] == edge_hash
        ]
        if len(node_exists) == 0:
            g.add_edge(
                annot_node_label,
                gene_node_label,
                label=Cons.GPROFILER_GENE_MIRNA_EDGE_LABEL,
                attr_dict=edge_attrs,
            )

    return g


def add_gprofiler_gene_reactome_subgraph(g, gene_node_label, annot_list):
    """Construct part of the graph by linking the gene to reactome pathways from GProfiler.

    :param g: the input graph to extend with new nodes and edges.
    :param gene_node_label: the gene node to be linked to reactome pathways from GProfiler.
    :param annot_list: list of pathways from GProfiler.
    :returns: a NetworkX MultiDiGraph
    """
    for annot in annot_list:
        if pd.isna(annot[Cons.NAME]):
            continue

        annot_node_label = annot[Cons.GPROFILER_ID]
        annot_node_attrs = _add_grofilier_genes_node_attrs(annot)
        annot_node_attrs[Cons.DATASOURCE] = annot[Cons.DATASOURCE]
        annot_node_attrs[Cons.GPROFILER_PATHWAY_TYPE] = Cons.REACTOME
        annot_node_attrs[Cons.LABEL] = Cons.PATHWAY_NODE_LABEL

        g.add_node(annot_node_label, attr_dict=annot_node_attrs)

        edge_attrs = Cons.GENE_PATHWAY_EDGE_ATTRS.copy()
        edge_attrs[Cons.DATASOURCE] = annot[Cons.DATASOURCE]

        edge_hash = hash(frozenset(edge_attrs.items()))
        edge_attrs[Cons.EDGE_HASH] = edge_hash
        edge_data = g.get_edge_data(gene_node_label, annot_node_label)
        edge_data = {} if edge_data is None else edge_data
        node_exists = [
            x for x, y in edge_data.items() if y["attr_dict"][Cons.EDGE_HASH] == edge_hash
        ]
        if len(node_exists) == 0:
            g.add_edge(
                gene_node_label,
                annot_node_label,
                label=Cons.GENE_PATHWAY_EDGE_LABEL,
                attr_dict=edge_attrs,
            )

    return g


def add_gprofiler_gene_transcription_factor_subgraph(g, gene_node_label, annot_list):
    """Construct part of the graph by linking the gene to transcription factors from GProfiler.

    :param g: the input graph to extend with new nodes and edges.
    :param gene_node_label: the gene node to be linked to transcription factors from GProfiler.
    :param annot_list: list of transcription factors from GProfiler.
    :returns: a NetworkX MultiDiGraph
    """
    for annot in annot_list:
        if pd.isna(annot[Cons.NAME]):
            continue

        annot_node_label = annot[Cons.GPROFILER_ID]
        annot_node_attrs = _add_grofilier_genes_node_attrs(annot)
        annot_node_attrs[Cons.DATASOURCE] = annot[Cons.DATASOURCE]
        annot_node_attrs[Cons.LABEL] = Cons.TRANS_FACTOR_NODE_LABEL

        g.add_node(annot_node_label, attr_dict=annot_node_attrs)

        edge_attrs = Cons.GPROFILER_EDGE_ATTRS.copy()
        edge_attrs[Cons.DATASOURCE] = annot[Cons.DATASOURCE]

        edge_hash = hash(frozenset(edge_attrs.items()))
        edge_attrs[Cons.EDGE_HASH] = edge_hash
        edge_data = g.get_edge_data(gene_node_label, annot_node_label)
        edge_data = {} if edge_data is None else edge_data
        node_exists = [
            x for x, y in edge_data.items() if y["attr_dict"][Cons.EDGE_HASH] == edge_hash
        ]
        if len(node_exists) == 0:
            g.add_edge(
                gene_node_label,
                annot_node_label,
                label=Cons.GPROFILER_GENE_TF_EDGE_LABEL,
                attr_dict=edge_attrs,
            )

    return g


def add_gprofiler_gene_gomf_subgraph(g, gene_node_label, annot_list):
    """Construct part of the graph by linking the gene to molecular functions from GProfiler.

    :param g: the input graph to extend with new nodes and edges.
    :param gene_node_label: the gene node to be linked to molecular functions from GProfiler.
    :param annot_list: list of molecular functions from GProfiler.
    :returns: a NetworkX MultiDiGraph
    """
    for annot in annot_list:
        if pd.isna(annot[Cons.NAME]):
            continue

        annot_node_label = annot[Cons.GPROFILER_ID]
        annot_node_attrs = _add_grofilier_genes_node_attrs(annot)
        annot_node_attrs[Cons.DATASOURCE] = annot[Cons.DATASOURCE]
        annot_node_attrs[Cons.LABEL] = Cons.GO_MF_NODE_LABEL

        g.add_node(annot_node_label, attr_dict=annot_node_attrs)

        edge_attrs = Cons.GPROFILER_EDGE_ATTRS.copy()
        edge_attrs[Cons.DATASOURCE] = annot[Cons.DATASOURCE]

        edge_hash = hash(frozenset(edge_attrs.items()))
        edge_attrs[Cons.EDGE_HASH] = edge_hash
        edge_data = g.get_edge_data(gene_node_label, annot_node_label)
        edge_data = {} if edge_data is None else edge_data
        node_exists = [
            x for x, y in edge_data.items() if y["attr_dict"][Cons.EDGE_HASH] == edge_hash
        ]
        if len(node_exists) == 0:
            g.add_edge(
                gene_node_label,
                annot_node_label,
                label=Cons.GENE_PATHWAY_EDGE_LABEL,
                attr_dict=edge_attrs,
            )

    return g


def add_gprofiler_gene_gocc_subgraph(g, gene_node_label, annot_list):
    """Construct part of the graph by linking the gene to cellular components from GProfiler.

    :param g: the input graph to extend with new nodes and edges.
    :param gene_node_label: the gene node to be linked to cellular components from GProfiler.
    :param annot_list: list of cellular components from GProfiler.
    :returns: a NetworkX MultiDiGraph
    """
    for annot in annot_list:
        if pd.isna(annot[Cons.NAME]):
            continue

        annot_node_label = annot[Cons.GPROFILER_ID]
        annot_node_attrs = _add_grofilier_genes_node_attrs(annot)
        annot_node_attrs[Cons.DATASOURCE] = annot[Cons.DATASOURCE]
        annot_node_attrs[Cons.LABEL] = Cons.GO_CC_NODE_LABEL

        g.add_node(annot_node_label, attr_dict=annot_node_attrs)

        edge_attrs = Cons.GPROFILER_EDGE_ATTRS.copy()
        edge_attrs[Cons.DATASOURCE] = annot[Cons.DATASOURCE]

        edge_hash = hash(frozenset(edge_attrs.items()))
        edge_attrs[Cons.EDGE_HASH] = edge_hash
        edge_data = g.get_edge_data(gene_node_label, annot_node_label)
        edge_data = {} if edge_data is None else edge_data
        node_exists = [
            x for x, y in edge_data.items() if y["attr_dict"][Cons.EDGE_HASH] == edge_hash
        ]
        if len(node_exists) == 0:
            g.add_edge(
                gene_node_label,
                annot_node_label,
                label=Cons.GENE_PATHWAY_EDGE_LABEL,
                attr_dict=edge_attrs,
            )

    return g


def add_gprofiler_gene_gobp_subgraph(g, gene_node_label, annot_list):
    """Construct part of the graph by linking the gene to biological processes from GProfiler.

    :param g: the input graph to extend with new nodes and edges.
    :param gene_node_label: the gene node to be linked to biological processes from GProfiler.
    :param annot_list: list of biological processes from GProfiler.
    :returns: a NetworkX MultiDiGraph
    """
    for annot in annot_list:
        if pd.isna(annot[Cons.NAME]):
            continue

        annot_node_label = annot[Cons.GPROFILER_ID]
        annot_node_attrs = _add_grofilier_genes_node_attrs(annot)
        annot_node_attrs[Cons.DATASOURCE] = annot[Cons.DATASOURCE]
        annot_node_attrs[Cons.LABEL] = Cons.GO_BP_NODE_LABEL

        g.add_node(annot_node_label, attr_dict=annot_node_attrs)

        edge_attrs = Cons.GPROFILER_EDGE_ATTRS.copy()
        edge_attrs[Cons.DATASOURCE] = annot[Cons.DATASOURCE]

        edge_hash = hash(frozenset(edge_attrs.items()))
        edge_attrs[Cons.EDGE_HASH] = edge_hash
        edge_data = g.get_edge_data(gene_node_label, annot_node_label)
        edge_data = {} if edge_data is None else edge_data
        node_exists = [
            x for x, y in edge_data.items() if y["attr_dict"][Cons.EDGE_HASH] == edge_hash
        ]
        if len(node_exists) == 0:
            g.add_edge(
                gene_node_label,
                annot_node_label,
                label=Cons.GENE_PATHWAY_EDGE_LABEL,
                attr_dict=edge_attrs,
            )

    return g


def add_gprofiler_gene_wikipathway_subgraph(g, gene_node_label, annot_list):
    """Construct part of the graph by linking the gene to pathways from WikiPathways.

    :param g: the input graph to extend with new nodes and edges.
    :param gene_node_label: the gene node to be linked to pathways from WikiPathways.
    :param annot_list: list of wikipathways from GProfiler.
    :returns: a NetworkX MultiDiGraph
    """
    for annot in annot_list:
        if pd.isna(annot[Cons.NAME]):
            continue

        annot_node_label = annot[Cons.GPROFILER_ID]
        annot_node_attrs = _add_grofilier_genes_node_attrs(annot)
        annot_node_attrs[Cons.GPROFILER_PATHWAY_TYPE] = Cons.WIKIPATHWAYS
        annot_node_attrs[Cons.DATASOURCE] = annot[Cons.DATASOURCE]
        annot_node_attrs[Cons.LABEL] = Cons.PATHWAY_NODE_LABEL

        g.add_node(annot_node_label, attr_dict=annot_node_attrs)

        edge_attrs = Cons.GPROFILER_EDGE_ATTRS.copy()
        edge_attrs[Cons.DATASOURCE] = annot[Cons.DATASOURCE]

        edge_hash = hash(frozenset(edge_attrs.items()))
        edge_attrs[Cons.EDGE_HASH] = edge_hash
        edge_data = g.get_edge_data(gene_node_label, annot_node_label)
        edge_data = {} if edge_data is None else edge_data
        node_exists = [
            x for x, y in edge_data.items() if y["attr_dict"][Cons.EDGE_HASH] == edge_hash
        ]
        if len(node_exists) == 0:
            g.add_edge(
                gene_node_label,
                annot_node_label,
                label=Cons.GENE_PATHWAY_EDGE_LABEL,
                attr_dict=edge_attrs,
            )

    return g


def add_mitocarta_gene_mito_subgraph(g, gene_node_label, annot_list):
    """Construct part of the graph by linking the gene to mitochondrial pathways from Mitocarta.

    :param g: the input graph to extend with new nodes and edges.
    :param gene_node_label: the gene node to be linked to mitochondrial pathways from Mitocarta.
    :param annot_list: list of mitochondrial pathways from Mitocarta.
    :returns: a NetworkX MultiDiGraph
    """
    for annot in annot_list:
        if pd.isna(annot[Cons.MITOCART_NODE_MAIN_LABEL]):
            continue

        annot_node_label = annot[Cons.MITOCART_NODE_MAIN_LABEL]
        annot_node_attrs = Cons.MITOCART_NODE_ATTRS.copy()
        annot_node_attrs.update(
            {
                Cons.NAME: annot[Cons.MITOCART_NODE_MAIN_LABEL],
                Cons.EVIDENCE: annot[Cons.EVIDENCE],
            }
        )
        if not pd.isna(annot[Cons.MITOCART_HPA_LOCATION]):
            annot_node_attrs[Cons.MITOCART_HPA_LOCATION] = annot[Cons.MITOCART_HPA_LOCATION]
        if not pd.isna(annot[Cons.MITOCART_SUB_MITO_LOCALIZATION]):
            annot_node_attrs[Cons.MITOCART_SUB_MITO_LOCALIZATION] = annot[
                Cons.MITOCART_SUB_MITO_LOCALIZATION
            ]

        g.add_node(annot_node_label, attr_dict=annot_node_attrs)

        edge_attrs = Cons.GENE_PATHWAY_EDGE_ATTRS.copy()
        edge_attrs[Cons.DATASOURCE] = Cons.MITOCARTA

        edge_hash = hash(frozenset(edge_attrs.items()))
        edge_attrs[Cons.EDGE_HASH] = edge_hash
        edge_data = g.get_edge_data(gene_node_label, annot_node_label)
        edge_data = {} if edge_data is None else edge_data
        node_exists = [
            x for x, y in edge_data.items() if y["attr_dict"][Cons.EDGE_HASH] == edge_hash
        ]
        if len(node_exists) == 0:
            g.add_edge(
                gene_node_label,
                annot_node_label,
                label=Cons.MITOCART_GENE_PATHWAY_EDGE_LABEL,
                attr_dict=edge_attrs,
            )

    return g


def add_opentargets_gene_reactome_pathway_subgraph(g, gene_node_label, annot_list):
    """Construct part of the graph by linking the gene to Reactome pathways.

    :param g: the input graph to extend with new nodes and edges.
    :param gene_node_label: the gene node to be linked to Reactome pathways.
    :param annot_list: list of Reactome pathways from OpenTargets.
    :returns: a NetworkX MultiDiGraph
    """
    logger.debug("Adding OpenTargets Reactome nodes and edges")
    for annot in annot_list:
        if pd.isna(annot[Cons.PATHWAY_ID]):
            continue

        annot_node_label = annot[Cons.OPENTARGETS_REACTOME_NODE_MAIN_LABEL]
        annot_node_attrs = Cons.OPENTARGETS_REACTOME_NODE_ATTRS.copy()
        annot_node_attrs.update(
            {
                Cons.DATASOURCE: Cons.OPENTARGETS,
                Cons.NAME: annot[Cons.PATHWAY_LABEL],
                Cons.ID: annot[Cons.PATHWAY_ID],
            }
        )

        g.add_node(annot_node_label, attr_dict=annot_node_attrs)

        edge_attrs = Cons.OPENTARGETS_GENE_REACTOME_EDGE_ATTRS.copy()

        edge_hash = hash(frozenset(edge_attrs.items()))
        edge_attrs[Cons.EDGE_HASH] = edge_hash
        edge_data = g.get_edge_data(gene_node_label, annot_node_label)
        edge_data = {} if edge_data is None else edge_data
        node_exists = [
            x for x, y in edge_data.items() if y["attr_dict"][Cons.EDGE_HASH] == edge_hash
        ]

        if len(node_exists) == 0:
            g.add_edge(
                gene_node_label,
                annot_node_label,
                label=Cons.GENE_PATHWAY_EDGE_LABEL,
                attr_dict=edge_attrs,
            )

    return g


def add_opentargets_gene_go_subgraph(g, gene_node_label, annot_list):
    """Construct part of the graph by linking the gene to gene ontologies.

    :param g: the input graph to extend with new nodes and edges.
    :param gene_node_label: the gene node to be linked to gene ontologies.
    :param annot_list: list of gene ontologies from OpenTargets.
    :returns: a NetworkX MultiDiGraph
    :raises ValueError: if the GO type is invalid.
    """
    logger.debug("Adding OpenTargets GO nodes and edges")
    for annot in annot_list:
        if pd.isna(annot[Cons.OPENTARGETS_GO_ID]):
            continue

        annot_node_label = annot[Cons.OPENTARGETS_GO_NODE_MAIN_LABEL]
        annot_node_attrs = Cons.OPENTARGETS_GO_NODE_ATTRS.copy()
        annot_node_attrs.update(
            {
                Cons.NAME: annot[Cons.OPENTARGETS_GO_NAME],
                Cons.ID: annot[Cons.OPENTARGETS_GO_ID],
                Cons.DATASOURCE: Cons.OPENTARGETS,
            }
        )

        if annot[Cons.OPENTARGETS_GO_TYPE] == "P":
            annot_node_attrs[Cons.LABEL] = Cons.GO_BP_NODE_LABEL
        elif annot[Cons.OPENTARGETS_GO_TYPE] == "F":
            annot_node_attrs[Cons.LABEL] = Cons.GO_MF_NODE_LABEL
        elif annot[Cons.OPENTARGETS_GO_TYPE] == "C":
            annot_node_attrs[Cons.LABEL] = Cons.GO_CC_NODE_LABEL
        else:
            raise ValueError(f"Invalid GO type: {annot[Cons.OPENTARGETS_GO_TYPE]}")

        g.add_node(annot_node_label, attr_dict=annot_node_attrs)

        edge_attrs = Cons.OPENTARGETS_GENE_GO_EDGE_ATTRS.copy()

        edge_hash = hash(frozenset(edge_attrs.items()))
        edge_attrs[Cons.EDGE_HASH] = edge_hash
        edge_data = g.get_edge_data(gene_node_label, annot_node_label)
        edge_data = {} if edge_data is None else edge_data
        node_exists = [
            x for x, y in edge_data.items() if y["attr_dict"][Cons.EDGE_HASH] == edge_hash
        ]

        if len(node_exists) == 0:
            g.add_edge(
                gene_node_label,
                annot_node_label,
                label=Cons.GENE_PATHWAY_EDGE_LABEL,
                attr_dict=edge_attrs,
            )

    return g


def add_opentargets_compound_side_effect_subgraph(g, compound_node_label, side_effects_list):
    """Construct part of the graph by linking the compound to side effects.

    :param g: the input graph to extend with new nodes and edges.
    :param compound_node_label: the compound node to be linked to side effect nodes.
    :param side_effects_list: list of side effects from OpenTargets.
    :returns: a NetworkX MultiDiGraph
    """
    logger.debug("Adding OpenTargets side effect nodes and edges")
    if not isinstance(side_effects_list, list):
        return g

    id_counter = 100
    nodes_list = g.nodes()

    for effect in side_effects_list:
        if pd.isna(effect[Cons.COMPOUND_SIDE_EFFECT_NODE_LABEL]):
            continue

        effect_node_label = effect[Cons.COMPOUND_SIDE_EFFECT_NODE_LABEL]

        # Adding a BDF id if the side effect node is not in the graph
        if effect_node_label not in nodes_list:
            effect_node_idx = f"{Cons.BIODATAFUSE}:{id_counter}"
            id_counter += 1

            # Add the side effect node to the graph
            effect_node_attrs = Cons.COMPOUND_SIDE_EFFECT_NODE_ATTRS.copy()
            effect_node_attrs.update(
                {
                    Cons.NAME: effect_node_label,
                    Cons.DATASOURCE: Cons.OPENTARGETS,
                    Cons.ID: effect_node_idx,
                }
            )
            g.add_node(effect_node_label, attr_dict=effect_node_attrs)

        # Add the edge between the compound and the side effect node
        edge_attrs = Cons.COMPOUND_SIDE_EFFECT_EDGE_ATTRS.copy()
        edge_hash = hash(frozenset(edge_attrs.items()))
        edge_attrs[Cons.EDGE_HASH] = edge_hash
        edge_data = g.get_edge_data(compound_node_label, effect_node_label)
        edge_data = {} if edge_data is None else edge_data
        node_exists = [
            x
            for x, y in edge_data.items()
            if "attr_dict" in y and y["attr_dict"].get(Cons.EDGE_HASH) == edge_hash
        ]

        if len(node_exists) == 0:
            g.add_edge(
                compound_node_label,
                effect_node_label,
                label=Cons.COMPOUND_SIDE_EFFECT_EDGE_LABEL,
                attr_dict=edge_attrs,
            )

    return g


def add_opentargets_gene_compound_subgraph(g, gene_node_label, annot_list):
    """Construct part of the graph by linking the gene to a list of compounds.

    :param g: the input graph to extend with new nodes and edges.
    :param gene_node_label: the gene node to be linked to compounds.
    :param annot_list: list of compounds from OpenTargets.
    :returns: a NetworkX MultiDiGraph
    """
    logger.debug("Adding OpenTargets compound nodes and edges")
    for annot in annot_list:
        if pd.isna(annot[Cons.OPENTARGETS_COMPOUND_RELATION]):
            continue

        if not pd.isna(annot[Cons.COMPOUND_NODE_MAIN_LABEL]):
            annot_node_label = annot[Cons.COMPOUND_NODE_MAIN_LABEL]
        else:
            annot_node_label = annot[Cons.CHEMBL_ID]

        annot_node_attrs = Cons.OPENTARGETS_COMPOUND_NODE_ATTRS.copy()
        annot_node_attrs.update(
            {
                Cons.NAME: annot_node_label,
                Cons.ID: annot[Cons.CHEMBL_ID],
                Cons.DATASOURCE: Cons.OPENTARGETS,
            }
        )

        other_info = {
            Cons.DRUGBANK_ID,
            Cons.OPENTARGETS_COMPOUND_CID,
            Cons.OPENTARGETS_COMPOUND_CLINICAL_TRIAL_PHASE,
            Cons.OPENTARGETS_COMPOUND_IS_APPROVED,
            Cons.OPENTARGETS_ADVERSE_EFFECT_COUNT,
        }

        for key in other_info:
            if not pd.isna(annot[key]):
                annot_node_attrs[key] = annot[key]

        merge_node(g, annot_node_label, annot_node_attrs)

        edge_attrs = Cons.OPENTARGETS_GENE_COMPOUND_EDGE_ATTRS.copy()
        edge_hash = hash(frozenset(edge_attrs.items()))
        edge_attrs[Cons.EDGE_HASH] = edge_hash
        edge_data = g.get_edge_data(annot_node_label, gene_node_label)
        edge_data = {} if edge_data is None else edge_data
        node_exists = [
            x for x, y in edge_data.items() if y["attr_dict"][Cons.EDGE_HASH] == edge_hash
        ]

        if len(node_exists) == 0:
            g.add_edge(
                annot_node_label,
                gene_node_label,
                label=annot[Cons.OPENTARGETS_COMPOUND_RELATION],
                attr_dict=edge_attrs,
            )

        # Add side effects
        if annot[Cons.COMPOUND_SIDE_EFFECT_NODE_MAIN_LABEL]:
            add_opentargets_compound_side_effect_subgraph(
                g, annot_node_label, annot[Cons.COMPOUND_SIDE_EFFECT_NODE_MAIN_LABEL]
            )

    return g


def add_molmedb_gene_inhibitor_subgraph(g, gene_node_label, annot_list):
    """Construct part of the graph by linking the gene to its inhibitors.

    :param g: the input graph to extend with new nodes and edges.
    :param gene_node_label: the gene node to be linked to its inhibitors.
    :param annot_list: list of gene inhibitors from MolMeDB.
    :returns: a NetworkX MultiDiGraph
    """
    logger.debug("Adding MolMeDB gene inhibitor nodes and edges")
    for annot in annot_list:
        if pd.isna(annot[Cons.MOLMEDB_COMPOUND_NAME]):
            continue

        if not pd.isna(annot[Cons.COMPOUND_NODE_MAIN_LABEL]):
            annot_node_label = annot[Cons.COMPOUND_NODE_MAIN_LABEL]
            annot_id = annot[Cons.COMPOUND_NODE_MAIN_LABEL]
        else:
            annot_node_label = annot[Cons.MOLMEDB_ID]
            annot_id = annot[Cons.MOLMEDB_ID]

        annot_node_attrs = Cons.MOLMEDB_COMPOUND_NODE_ATTRS.copy()
        annot_node_attrs.update(
            {
                Cons.NAME: annot[Cons.MOLMEDB_COMPOUND_NAME],
                Cons.ID: annot_id,
                Cons.MOLMEDB_ID: annot[Cons.MOLMEDB_ID],
                Cons.MOLMEDB_INCHIKEY: annot[Cons.MOLMEDB_INCHIKEY],
                Cons.MOLMEDB_SMILES: annot[Cons.MOLMEDB_SMILES],
                Cons.SOURCE_PMID: annot[Cons.SOURCE_PMID],
                Cons.DATASOURCE: Cons.MOLMEDB,
            }
        )

        merge_node(g, annot_node_label, annot_node_attrs)

        edge_attrs = Cons.MOLMEDB_PROTEIN_COMPOUND_EDGE_ATTRS.copy()

        edge_hash = hash(frozenset(edge_attrs.items()))
        edge_attrs[Cons.EDGE_HASH] = edge_hash
        edge_data = g.get_edge_data(gene_node_label, annot_node_label)
        edge_data = {} if edge_data is None else edge_data
        node_exists = [
            x for x, y in edge_data.items() if y["attr_dict"][Cons.EDGE_HASH] == edge_hash
        ]

        if len(node_exists) == 0:
            g.add_edge(
                annot_node_label,
                gene_node_label,
                label=Cons.MOLMEDB_PROTEIN_COMPOUND_EDGE_LABEL,
                attr_dict=edge_attrs,
            )

    return g


# TODO: Fix this function, looks like it adds compounds nodes instead of transporter nodes
def add_molmedb_compound_gene_subgraph(g, compound_node_label, annot_list):
    """Construct part of the graph by linking the compound to inhibited genes.

    :param g: the input graph to extend with new nodes and edges.
    :param compound_node_label: the compound node to be linked to its inhibitors.
    :param annot_list: list of gene inhibitors from MolMeDB.
    :returns: a NetworkX MultiDiGraph
    """
    logger.debug("Adding MolMeDB compound gene nodes and edges")
    for annot in annot_list:
        if pd.isna(annot[Cons.MOLMEDB_COMPOUND_NAME]):
            continue

        if not pd.isna(annot[Cons.COMPOUND_NODE_MAIN_LABEL]):
            annot_node_label = annot[Cons.COMPOUND_NODE_MAIN_LABEL]
        else:
            annot_node_label = annot["molmedb_id"]

        annot_node_attrs = Cons.MOLMEDB_COMPOUND_NODE_ATTRS.copy()
        annot_node_attrs.update(
            {
                "name": annot["compound_name"],
                "id": annot["molmedb_id"],
                "datasource": Cons.MOLMEDB,
            }
        )

        if not pd.isna(annot[Cons.COMPOUND_NODE_MAIN_LABEL]):
            annot_node_attrs["id"] = annot[Cons.COMPOUND_NODE_MAIN_LABEL]
        else:
            annot_node_attrs["id"] = annot["molmedb_id"]

        other_info = {
            "inchikey": annot["inchikey"],
            "smiles": annot["smiles"],
            "compound_cid": annot["compound_cid"],
            "chebi_id": annot["chebi_id"],
            "drugbank_id": annot["drugbank_id"],
            "source_pmid": annot["source_pmid"],
            "uniprot_trembl_id": annot["uniprot_trembl_id"],
            # "pdb_ligand_id": annot["pdb_ligand_id"],
        }

        for key, value in other_info.items():
            if not pd.isna(value):
                annot_node_attrs[key] = value

        merge_node(g, annot_node_label, annot_node_attrs)

        edge_attrs = Cons.MOLMEDB_PROTEIN_COMPOUND_EDGE_ATTRS.copy()

        edge_hash = hash(frozenset(edge_attrs.items()))
        edge_attrs["edge_hash"] = edge_hash  # type: ignore
        edge_data = g.get_edge_data(compound_node_label, annot_node_label)
        edge_data = {} if edge_data is None else edge_data
        node_exists = [x for x, y in edge_data.items() if y["attr_dict"]["edge_hash"] == edge_hash]

        if len(node_exists) == 0:
            g.add_edge(
                annot_node_label,
                compound_node_label,
                label=Cons.MOLMEDB_PROTEIN_COMPOUND_EDGE_LABEL,
                attr_dict=edge_attrs,
            )

    return g


def add_pubchem_assay_subgraph(g, gene_node_label, annot_list):
    """Construct part of the graph by linking the gene to a list of compounds tested on it.

    :param g: the input graph to extend with new nodes and edges.
    :param gene_node_label: the gene node to be linked to compound tested on it.
    :param annot_list: list of compounds tested on gene from PubChem.
    :returns: a NetworkX MultiDiGraph
    """
    logger.debug("Adding PubChem assay nodes and edges")
    for annot in annot_list:
        if pd.isna(annot[Cons.PUBCHEM_ASSAY_ID]):
            continue

        annot_node_label = annot[Cons.COMPOUND_NODE_MAIN_LABEL]
        annot_node_attrs = Cons.PUBCHEM_COMPOUND_NODE_ATTRS.copy()
        annot_node_attrs.update(
            {
                Cons.NAME: annot["compound_name"],
                Cons.ID: annot["compound_cid"],
                Cons.INCHI: annot["inchi"],
            }
        )
        if not pd.isna(annot["smiles"]):
            annot_node_attrs[Cons.SMILES] = annot["smiles"]

        # g.add_node(annot_node_label, attr_dict=annot_node_attrs)
        merge_node(g, annot_node_label, annot_node_attrs)

        edge_attrs = Cons.PUBCHEM_GENE_COMPOUND_EDGE_ATTRS.copy()
        edge_attrs.update(
            {
                Cons.PUBCHEM_ASSAY_TYPE: annot[Cons.PUBCHEM_ASSAY_TYPE],
                Cons.PUBCHEM_ASSAY_ID: annot[Cons.PUBCHEM_ASSAY_ID],
            }
        )

        edge_hash = hash(frozenset(edge_attrs.items()))
        edge_attrs[Cons.EDGE_HASH] = edge_hash
        edge_data = g.get_edge_data(gene_node_label, annot_node_label)
        edge_data = {} if edge_data is None else edge_data
        node_exists = [
            x for x, y in edge_data.items() if y["attr_dict"][Cons.EDGE_HASH] == edge_hash
        ]

        if len(node_exists) == 0:
            g.add_edge(
                annot_node_label,
                gene_node_label,
                label=Cons.PUBCHEM_EDGE_LABEL_MAPPER[annot["outcome"]],
                attr_dict=edge_attrs,
            )

    return g


def add_stringdb_ppi_subgraph(g, gene_node_label, annot_list):
    """Construct part of the graph by linking the gene to genes.

    :param g: the input graph to extend with new nodes and edges.
    :param gene_node_label: the gene node to be linked to other genes entities.
    :param annot_list: list of protein-protein interactions from StringDb.
    :returns: a NetworkX MultiDiGraph
    """
    logger.debug("Adding StringDb PPI nodes and edges")
    for ppi in annot_list:
        edge_attrs = Cons.STRING_PPI_EDGE_ATTRS.copy()
        edge_attrs[Cons.STRING_PPI_SCORE] = ppi[Cons.STRING_PPI_SCORE]

        edge_hash = hash(frozenset(edge_attrs.items()))
        edge_attrs[Cons.EDGE_HASH] = edge_hash  # type: ignore
        edge_data = g.get_edge_data(gene_node_label, ppi[Cons.STRING_PPI_INTERACTS_WITH])

        edge_data = {} if edge_data is None else edge_data
        node_exists = [
            x for x, y in edge_data.items() if y["attr_dict"][Cons.EDGE_HASH] == edge_hash
        ]
        if len(node_exists) == 0 and not pd.isna(ppi[Cons.STRING_PPI_INTERACTS_WITH]):
            g.add_edge(
                gene_node_label,
                ppi[Cons.STRING_PPI_INTERACTS_WITH],
                label=Cons.STRING_PPI_EDGE_MAIN_LABEL,
                attr_dict=edge_attrs,
            )

            g.add_edge(
                ppi[Cons.STRING_PPI_INTERACTS_WITH],
                gene_node_label,
                label=Cons.STRING_PPI_EDGE_MAIN_LABEL,
                attr_dict=edge_attrs,
            )

    return g


def add_opentargets_disease_compound_subgraph(g, disease_node, annot_list):
    """Construct part of the graph by linking the disease to compounds.

    :param g: the input graph to extend with new nodes and edges.
    :param disease_node: the disease node to be linked to compounds.
    :param annot_list: list of compounds from OpenTargets.
    :returns: a NetworkX MultiDiGraph
    """
    logger.debug("Adding OpenTargets disease compound nodes and edges")
    for annot in annot_list:
        if pd.isna(annot[Cons.OPENTARGETS_COMPOUND_RELATION]):
            continue

        if pd.isna(annot[Cons.COMPOUND_NODE_MAIN_LABEL]):
            annot_node_label = annot[Cons.CHEMBL_ID]
        else:
            annot_node_label = annot[Cons.COMPOUND_NODE_MAIN_LABEL]

        # create compound node and merge with existing node if it exists
        annot_node_attrs = Cons.OPENTARGETS_COMPOUND_NODE_ATTRS.copy()
        annot_node_attrs.update(
            {
                Cons.NAME: annot_node_label,
                Cons.ID: annot[Cons.CHEMBL_ID],
                Cons.DATASOURCE: Cons.OPENTARGETS,
            }
        )

        other_info = {
            Cons.DRUGBANK_ID,
            Cons.OPENTARGETS_COMPOUND_CID,
            Cons.OPENTARGETS_COMPOUND_CLINICAL_TRIAL_PHASE,
            Cons.OPENTARGETS_COMPOUND_IS_APPROVED,
            Cons.OPENTARGETS_ADVERSE_EFFECT_COUNT,
        }

        for key in other_info:
            if not pd.isna(annot[key]):
                annot_node_attrs[key] = annot[key]

        merge_node(g, annot_node_label, annot_node_attrs)

        edge_attrs = Cons.OPENTARGETS_DISEASE_COMPOUND_EDGE_ATTRS.copy()
        edge_attrs[Cons.LABEL] = annot[Cons.OPENTARGETS_COMPOUND_RELATION]
        edge_hash = hash(frozenset(edge_attrs.items()))
        edge_attrs[Cons.EDGE_HASH] = edge_hash

        edge_data = g.get_edge_data(annot_node_label, disease_node)
        edge_data = {} if edge_data is None else edge_data
        node_exists = [
            x for x, y in edge_data.items() if y["attr_dict"][Cons.EDGE_HASH] == edge_hash
        ]

        if len(node_exists) == 0:
            g.add_edge(
                annot_node_label,
                disease_node,
                label=annot[Cons.OPENTARGETS_COMPOUND_RELATION],
                attr_dict=edge_attrs,
            )

        # Add side effects
        if annot[Cons.OPENTARGETS_ADVERSE_EFFECT]:
            add_opentargets_compound_side_effect_subgraph(
                g, annot_node_label, annot[Cons.OPENTARGETS_ADVERSE_EFFECT]
            )

    return g


def add_wikipathways_molecular_subgraph(g, gene_node_label, annot_list):
    """Construct part of the graph by linking molecular entities from WP with MIMs.

    :param g: the input graph to extend with new nodes and edges.
    :param gene_node_label: the disease node to be linked to compounds.
    :param annot_list: result of querying WP for molecular interactions.
    :returns: a NetworkX MultiDiGraph
    """
    logger.debug("Adding WikiPathways molecular nodes and edges")
    for annot in annot_list:
        for target_key in [Cons.WIKIPATHWAYS_TARGET_GENE, Cons.WIKIPATHWAYS_TARGET_METABOLITE]:
            target = annot.get(target_key)
            target_node_label = None

            if target and target != gene_node_label:  # No interactions with self
                target_node_label = str(target)

            if target_node_label is None:
                continue

            target_node_label = target_node_label.replace(f"{Cons.WIKIPATHWAYS_TARGET_GENE}:", "")

            interaction_type = annot.get(Cons.WIKIPATHWAYS_MIM_TYPE, "Interaction")
            edge_attrs = Cons.MOLECULAR_INTERACTION_EDGE_ATTRS.copy()
            edge_attrs.update(
                {
                    Cons.WIKIPATHWAYS_INTERACTION_TYPE: interaction_type,
                    Cons.WIKIPATHWAYS_RHEA_ID: annot.get(Cons.WIKIPATHWAYS_RHEA_ID, ""),
                }
            )
            edge_attrs[Cons.EDGE_HASH] = hash(frozenset(edge_attrs.items()))  # type: ignore

            if not g.has_node(target_node_label):
                node_attrs = Cons.MOLECULAR_PATHWAY_NODE_ATTRS.copy()
                node_attrs.update(
                    {
                        Cons.ID: annot.get(Cons.PATHWAY_ID, ""),
                        Cons.NAME: annot.get(Cons.PATHWAY_LABEL, ""),
                        Cons.ID: target_node_label,
                    }
                )
                g.add_node(target_node_label, attr_dict=node_attrs)

            edge_exists = False
            if g.has_edge(gene_node_label, target_node_label):
                edge_data = g.get_edge_data(gene_node_label, target_node_label)
                for edge_key in edge_data:
                    if (
                        edge_data[edge_key].get("attr_dict", {}).get(Cons.EDGE_HASH)
                        == edge_attrs[Cons.EDGE_HASH]
                    ):
                        edge_exists = True
                        break

            if not edge_exists:
                g.add_edge(
                    gene_node_label,
                    target_node_label,
                    label=Cons.GENE_PATHWAY_EDGE_LABEL,
                    attr_dict=edge_attrs,
                )
    return g


def add_aopwiki_subgraph(g, entity_node_label, annot_list):
    """Construct part of the graph by linking the gene to AOP entities.

    :param g: the input graph to extend with new nodes and edges.
    :param entity_node_label: the gene node to be linked to AOP entities.
    :param annot_list: list of AOPWIKI Key Events.
    :returns: a NetworkX MultiDiGraph
    """
    for annot in annot_list:
        # Add AOP node
        if annot[Cons.AOP_NODE_MAIN_LABEL]:
            aop_node_label = f"{Cons.AOP_PATHWAY}:{annot.get(Cons.AOP_NODE_MAIN_LABEL)}"
            aop_node_attrs = Cons.AOPWIKI_NODE_ATTRS.copy()
            aop_node_attrs.update(
                {
                    Cons.ID: aop_node_label,
                    Cons.NAME: annot.get("aop_title", "Unknown"),
                    Cons.LABEL: Cons.AOP_NODE_LABEL,
                }
            )
            g.add_node(aop_node_label, attr_dict=aop_node_attrs)

            # Connect gene to AOP node
            edge_attrs = Cons.AOPWIKI_EDGE_ATTRS.copy()
            edge_attrs[Cons.EDGE_HASH] = hash(frozenset(edge_attrs.items()))
            if not edge_exists(g, entity_node_label, aop_node_label, edge_attrs):
                g.add_edge(
                    entity_node_label,
                    aop_node_label,
                    label=Cons.AOP_EDGE_LABEL,
                    attr_dict=edge_attrs,
                )

        # Add MIE node
        if annot[Cons.MIE_NODE_MAIN_LABEL]:
            mie_node_label = f"{Cons.MOL_INITIATING_EVENT}:{annot.get(Cons.MIE_NODE_MAIN_LABEL)}"
            mie_node_attrs = Cons.AOPWIKI_NODE_ATTRS.copy()
            mie_node_attrs.update(
                {
                    Cons.ID: mie_node_label,
                    Cons.NAME: annot.get("MIE_title", "Unknown"),
                    Cons.LABEL: Cons.MIE_NODE_LABEL,
                }
            )
            g.add_node(mie_node_label, attr_dict=mie_node_attrs)

            # Connect MIE to AOP node
            if aop_node_label:
                edge_attrs = Cons.AOPWIKI_EDGE_ATTRS.copy()
                edge_attrs[Cons.EDGE_HASH] = hash(frozenset(edge_attrs.items()))
                if not edge_exists(g, mie_node_label, aop_node_label, edge_attrs):
                    g.add_edge(
                        mie_node_label,
                        aop_node_label,
                        label=Cons.MIE_AOP_EDGE_LABEL,
                        attr_dict=edge_attrs,
                    )

        # Add KE upstream node
        if annot[Cons.KEY_EVENT_UPSTREAM_NODE_MAIN_LABEL]:
            ke_upstream_node_label = (
                f"{Cons.KEY_EVENT}:{annot.get(Cons.KEY_EVENT_UPSTREAM_NODE_MAIN_LABEL)}"
            )
            ke_upstream_node_attrs = Cons.AOPWIKI_NODE_ATTRS.copy()
            ke_upstream_node_attrs.update(
                {
                    Cons.ID: ke_upstream_node_label,
                    Cons.NAME: annot.get("KE_upstream_title", "Unknown"),
                    Cons.LABEL: Cons.KEY_EVENT_NODE_LABEL,
                    "organ": annot.get("KE_upstream_organ", ""),
                }
            )
            ke_upstream_node_attrs[Cons.LABEL] = Cons.KEY_EVENT_NODE_LABEL
            g.add_node(ke_upstream_node_label, attr_dict=ke_upstream_node_attrs)

            # Connect KE upstream to MIE node
            if mie_node_label:
                edge_attrs = Cons.AOPWIKI_EDGE_ATTRS.copy()
                edge_attrs[Cons.EDGE_HASH] = hash(frozenset(edge_attrs.items()))
                if not edge_exists(g, ke_upstream_node_label, mie_node_label, edge_attrs):
                    g.add_edge(
                        ke_upstream_node_label,
                        mie_node_label,
                        label=Cons.KE_UPSTREAM_MIE_EDGE_LABEL,
                        attr_dict=edge_attrs,
                    )

        # Add KE downstream node
        if annot[Cons.KEY_EVENT_DOWNSTREAM_NODE_MAIN_LABEL]:
            ke_downstream_node_label = (
                f"{Cons.KEY_EVENT}:{annot.get(Cons.KEY_EVENT_DOWNSTREAM_NODE_MAIN_LABEL)}"
            )
            ke_downstream_node_attrs = Cons.AOPWIKI_NODE_ATTRS.copy()
            ke_downstream_node_attrs.update(
                {
                    Cons.ID: ke_downstream_node_label,
                    Cons.NAME: annot.get("KE_downstream_title", "Unknown"),
                    Cons.LABEL: Cons.KEY_EVENT_NODE_LABEL,
                    "organ": annot.get("KE_downstream_organ", ""),
                }
            )
            g.add_node(ke_downstream_node_label, attr_dict=ke_downstream_node_attrs)

            # Connect KE downstream to KE upstream node
            if ke_upstream_node_label:
                edge_attrs = Cons.AOPWIKI_EDGE_ATTRS.copy()
                edge_attrs[Cons.EDGE_HASH] = hash(frozenset(edge_attrs.items()))
                if not edge_exists(g, ke_upstream_node_label, ke_downstream_node_label, edge_attrs):
                    g.add_edge(
                        ke_upstream_node_label,
                        ke_downstream_node_label,
                        label=Cons.KE_DOWNSTREAM_KE_EDGE_LABEL,
                        attr_dict=edge_attrs,
                    )

        # Add AO node
        if annot[Cons.AO_NODE_MAIN_LABEL]:
            ao_node_label = f"{Cons.ADVERSE_OUTCOME}:{annot.get(Cons.AO_NODE_MAIN_LABEL)}"
            ao_node_attrs = Cons.AOPWIKI_NODE_ATTRS.copy()
            ao_node_attrs.update(
                {
                    Cons.ID: ao_node_label,
                    Cons.NAME: annot.get("ao_title", "Unknown"),
                    Cons.LABEL: Cons.AO_NODE_LABEL,
                }
            )
            g.add_node(ao_node_label, attr_dict=ao_node_attrs)

            # Connect AO directly to KE upstream node
            if ke_upstream_node_label:
                edge_attrs = Cons.AOPWIKI_EDGE_ATTRS.copy()
                edge_attrs[Cons.EDGE_HASH] = hash(frozenset(edge_attrs.items()))
                if not edge_exists(g, ke_downstream_node_label, ao_node_label, edge_attrs):
                    g.add_edge(
                        ke_upstream_node_label,
                        ao_node_label,
                        label=Cons.AO_KE_EDGE_LABEL,
                        attr_dict=edge_attrs,
                    )

    return g


def edge_exists(g, source, target, edge_attrs):
    """Check if an edge with the same attributes already exists in the graph.

    :param g: the input graph.
    :param source: the source node of the edge.
    :param target: the target node of the edge.
    :param edge_attrs: the attributes of the edge to check.
    :returns: True if the edge exists, False otherwise.
    """
    if g.has_edge(source, target):
        edge_data = g.get_edge_data(source, target)
        for edge_key in edge_data:
            if edge_data[edge_key].get("attr_dict", {}).get("edge_hash") == edge_attrs["edge_hash"]:
                return True
    return False


def add_ensembl_homolog_subgraph(g, gene_node_label, annot_list):
    """Construct part of the graph by linking the gene to genes.

    :param g: the input graph to extend with new nodes and edges.
    :param gene_node_label: the gene node to be linked to other genes entities.
    :param annot_list: list of homologs from Ensembl.
    :returns: a NetworkX MultiDiGraph
    """
    logger.debug("Adding Ensembl homolog nodes and edges")
    for hl in annot_list:
        edge_attrs = Cons.ENSEMBL_HOMOLOG_EDGE_ATTRS.copy()

        edge_hash = hash(frozenset(edge_attrs.items()))
        edge_attrs["edge_hash"] = edge_hash
        edge_data = g.get_edge_data(gene_node_label, hl[Cons.ENSEMBL_HOMOLOG_MAIN_LABEL])

        edge_data = {} if edge_data is None else edge_data
        node_exists = [x for x, y in edge_data.items() if y["attr_dict"]["edge_hash"] == edge_hash]
        if len(node_exists) == 0 and not pd.isna(hl[Cons.ENSEMBL_HOMOLOG_MAIN_LABEL]):
            g.add_edge(
                gene_node_label,
                hl[Cons.ENSEMBL_HOMOLOG_MAIN_LABEL],
                label=Cons.ENSEMBL_HOMOLOG_EDGE_LABEL,
                attr_dict=edge_attrs,
            )

    return g


def add_tflink_gene_tf_subgraph(g, gene_node_label, annot_list):
    """Construct part of the graph by linking the gene to genes.

    :param g: the input graph to extend with new nodes and edges.
    :param gene_node_label: the gene node to be linked to other genes entities.
    :param annot_list: list of protein-protein interactions from StringDb.
    :returns: a NetworkX MultiDiGraph
    """
    logger.debug("Adding TFLink gene-TF nodes and edges")
    for tf in annot_list:
        edge_attrs = Cons.TFLINK_GENE_TF_EDGE_ATTRS.copy()
        edge_attrs.update(
            {
                Cons.NAME_TARGET: tf["Name.Target"],
                Cons.UNIPROTID_TARGET: tf["UniprotID.Target"],
                Cons.DETECTION_METHOD: tf["Detection.method"],
                Cons.PUBMEDID: tf["PubmedID"],
                Cons.SOURCE_DATABASE: tf["Source.database"],
                Cons.SMALL_SCALE_EVIDENCE: tf["Small-scale.evidence"],
            }
        )

        edge_hash = hash(frozenset(edge_attrs.items()))
        edge_attrs[Cons.EDGE_HASH] = edge_hash

        edge_data = g.get_edge_data(tf[Cons.ENSEMBL_GENE_ID_TARGET], gene_node_label)
        edge_data = {} if edge_data is None else edge_data
        node_exists = [
            x for x, y in edge_data.items() if y["attr_dict"][Cons.EDGE_HASH] == edge_hash
        ]
        if len(node_exists) == 0 and not pd.isna(tf[Cons.ENSEMBL_GENE_ID_TARGET]):
            g.add_edge(
                gene_node_label,
                tf[Cons.ENSEMBL_GENE_ID_TARGET],
                label=Cons.TFLINK_EDGE_LABEL,
                attr_dict=edge_attrs,
            )
    return g


def add_gene_node(g, row, dea_columns):
    """Add gene node from each row of the combined_df to the graph.

    :param g: the input graph to extend with gene nodes.
    :param row: row in the combined DataFrame.
    :param dea_columns: list of dea_columns.
    :returns: label for gene node
    """
    gene_node_label = row[Cons.IDENTIFIER_COL]
    gene_node_attrs = {
        Cons.DATASOURCE: Cons.BRIDGEDB,
        Cons.NAME: f"{row[Cons.IDENTIFIER_SOURCE_COL]}:{row[Cons.IDENTIFIER_COL]}",
        Cons.ID: f"{row[Cons.TARGET_SOURCE_COL]}:{row[Cons.TARGET_COL]}",
        Cons.LABEL: Cons.GENE_NODE_LABEL,
        Cons.LABEL: Cons.GENE_NODE_LABEL,  # Ensure label
        row[Cons.TARGET_SOURCE_COL]: row[Cons.TARGET_COL],
    }
    if "is_tf" in row:
        gene_node_attrs["is_tf"] = row["is_tf"]
    if "is_target" in row:
        gene_node_attrs["is_target"] = row["is_target"]
    for c in dea_columns:
        gene_node_attrs[c[:-4]] = row[c]

    g.add_node(gene_node_label, attr_dict=gene_node_attrs)
    return gene_node_label


def add_compound_node(g, row):
    """Add compound node from each row of the combined_df to the graph.

    :param g: the input graph to extend with compound nodes.
    :param row: row in the combined DataFrame.
    :returns: label for compound node
    """
    compound_node_label = row["identifier"]
    compound_node_attrs = {
        Cons.DATASOURCE: Cons.BRIDGEDB,
        Cons.NAME: f"{row[Cons.IDENTIFIER_SOURCE_COL]}:{row[Cons.IDENTIFIER_COL]}",
        Cons.ID: row[Cons.TARGET_COL],
        Cons.LABEL: Cons.COMPOUND_NODE_LABEL,
        Cons.LABEL: Cons.COMPOUND_NODE_LABEL,  # Ensure label
        row[Cons.TARGET_SOURCE_COL]: f"{row[Cons.TARGET_SOURCE_COL]}:{row[Cons.TARGET_COL]}",
    }

    g.add_node(compound_node_label, attr_dict=compound_node_attrs)
    return compound_node_label


"""Processing specific annotation types"""


def process_annotations(g, gene_node_label, row, func_dict):
    """Process the annotations for gene node from each row of the combined_df to the graph.

    :param g: the input graph to extend with gene nodes.
    :param gene_node_label: the gene node to be linked to annotation entities.
    :param row: row in the combined DataFrame.
    :param func_dict: dictionary of subgraph function.
    """
    for annot_key in func_dict:
        if annot_key not in row:
            continue

        annot_list = row[annot_key]

        if isinstance(annot_list, np.ndarray):
            annot_list = annot_list.tolist()
        elif not isinstance(annot_list, list):
            annot_list = []

        func_dict[annot_key](g, gene_node_label, annot_list)


def process_disease_compound(g, disease_compound, disease_nodes):
    """Process disease-compound relationships and add them to the graph.

    :param g: the input graph to extend with gene nodes.
    :param disease_compound: the input DataFrame containing disease_compound relationships.
    :param disease_nodes: the input dictionary containing disease nodes.
    """
    for _i, row in disease_compound.iterrows():
        disease_node_id = row[Cons.TARGET_COL].replace("_", ":")  # disease node label

        # Skip disease not in the graph
        if disease_node_id not in disease_nodes:
            annot_node_attrs = Cons.OPENTARGET_DISEASE_NODE_ATTRS.copy()
            annot_node_attrs.update(
                {
                    Cons.NAME: disease_node_id,
                    Cons.ID: disease_node_id,
                }
            )

            g.add_node(disease_node_id, attr_dict=annot_node_attrs)
        else:
            disease_node_id = disease_nodes[
                disease_node_id
            ]  # Convert the EFO to existing node label

        compound_annot_list = row[Cons.OPENTARGETS_DISEASE_COMPOUND_COL]

        if isinstance(compound_annot_list, float):
            compound_annot_list = []
        elif isinstance(compound_annot_list, np.ndarray):
            compound_annot_list = compound_annot_list.tolist()
        elif not isinstance(compound_annot_list, list):
            logger.warning(
                f"compound_annot_list of type {type(compound_annot_list)} and not list. Skipping..."
            )
            compound_annot_list = []

        add_opentargets_disease_compound_subgraph(g, disease_node_id, compound_annot_list)


def process_ppi(g, gene_node_label, row):
    """Process protein-protein interactions and add them to the graph.

    :param g: the input graph to extend with gene nodes.
    :param gene_node_label: the gene node to be linked to annotation entities.
    :param row: row in the combined DataFrame.
    """
    if Cons.STRING_INTERACT_COL in row and row[Cons.STRING_INTERACT_COL] is not None:
        try:
            ppi_list = json.loads(json.dumps(row[Cons.STRING_INTERACT_COL]))
        except (ValueError, TypeError):
            ppi_list = []

        if isinstance(ppi_list, list) and len(ppi_list) > 0:
            valid_ppi_list = [
                item for item in ppi_list if pd.notna(item.get(Cons.STRING_PPI_EDGE_MAIN_LABEL))
            ]
            if valid_ppi_list:
                add_stringdb_ppi_subgraph(g, gene_node_label, valid_ppi_list)

        if not isinstance(ppi_list, float):
            for item in ppi_list:
                if pd.isna(item["stringdb_link_to"]):
                    ppi_list = []
            add_stringdb_ppi_subgraph(g, gene_node_label, ppi_list)


def process_tf_target(g, gene_node_label, row):
    """Process tf-target interactions and add them to the graph.

    :param g: the input graph to extend with gene nodes.
    :param gene_node_label: the gene node to be linked to annotation entities.
    :param row: row in the combined DataFrame.
    """
    if Cons.ITS_TARGET_COL in row and row[Cons.ITS_TARGET_COL] is not None:
        its_target_list = json.loads(json.dumps(row[Cons.ITS_TARGET_COL]))

        if isinstance(its_target_list, list) and len(its_target_list) > 0:
            add_tflink_gene_tf_subgraph(g, gene_node_label, its_target_list)


def process_homologs(g, combined_df, homolog_df_list, func_dict, dea_columns):
    """Process homolog dataframes and combined df and add them to the graph.

    :param g: the input graph to extend with gene nodes.
    :param combined_df: dataframe without homolog information.
    :param homolog_df_list: list of dataframes from homolog queries.
    :param func_dict: list of functions for node generation.
    :param dea_columns: columns ending with _dea
    """
    func_dict_hl = {}

    homolog_cols = set()
    for homolog_df in homolog_df_list:
        if homolog_df is not None and not homolog_df.empty:
            homolog_cols.update(homolog_df.columns)

    func_dict_hl = {
        key: func
        for key, func in func_dict.items()
        if key in homolog_cols and key in combined_df.columns
    }

    for _i, row in tqdm(combined_df.iterrows(), total=combined_df.shape[0], desc="Building graph"):
        if pd.isna(row["identifier"]) or pd.isna(row["target"]):
            continue
        gene_node_label = add_gene_node(g, row, dea_columns)
        homolog_keys = set(func_dict_hl.keys())
        func_dict_non_hl = {key: func for key, func in func_dict.items() if key not in homolog_keys}
        process_annotations(g, gene_node_label, row, func_dict_non_hl)

    for _i, row in tqdm(combined_df.iterrows(), total=combined_df.shape[0]):
        if pd.isna(row["identifier"]) or pd.isna(row["Ensembl_homologs"]):
            continue

        homologs = row["Ensembl_homologs"]

        if isinstance(homologs, list) and homologs:
            for homolog_entry in homologs:
                homolog_node_label = homolog_entry.get("homolog")

                if pd.isna(homolog_node_label) or homolog_node_label == "nan":
                    continue

                if homolog_node_label:
                    annot_node_attrs = Cons.ENSEMBL_HOMOLOG_NODE_ATTRS.copy()
                    annot_node_attrs["id"] = homolog_node_label
                    annot_node_attrs[Cons.LABEL] = Cons.HOMOLOG_NODE_LABEL
                    g.add_node(homolog_node_label, attr_dict=annot_node_attrs)

                    homolog_row = pd.Series(
                        {
                            "identifier": homolog_node_label,
                            **{key: row.get(key) for key in func_dict_hl},
                        }
                    )

                    process_annotations(g, homolog_node_label, homolog_row, func_dict_hl)


def normalize_node_attributes(g):
    """Normalize node attributes by flattening the 'attr_dict'.

    :param g: the input graph to extend with gene nodes.
    """
    for node in g.nodes():
        if "attr_dict" in g.nodes[node]:
            for k, v in g.nodes[node]["attr_dict"].items():
                if v is not None:
                    g.nodes[node][k] = v

            del g.nodes[node]["attr_dict"]
        # Ensure 'labels' is present after flattening
        if Cons.LABEL not in g.nodes[node]:
            g.nodes[node][Cons.LABEL] = g.nodes[node].get("label", "Unknown")


def normalize_edge_attributes(g):
    """Normalize edge attributes by flattening the 'attr_dict'.

    :param g: the input graph to extend with gene nodes.
    """
    for u, v, k in g.edges(keys=True):
        if "attr_dict" in g[u][v][k]:
            for x, y in g[u][v][k]["attr_dict"].items():
                if y is not None and x != "edge_hash":
                    g[u][v][k][x] = y

            del g[u][v][k]["attr_dict"]


def _built_gene_based_graph(
    g: nx.MultiDiGraph,
    combined_df: pd.DataFrame,
    disease_compound=None,
    pathway_compound=None,
    homolog_df_list=None,
):
    """Build a gene-based graph."""
    combined_df = combined_df[(combined_df[Cons.TARGET_SOURCE_COL] == Cons.ENSEMBL)]

    dea_columns = [c for c in combined_df.columns if c.endswith("_dea")]

    func_dict = {
        Cons.BGEE_GENE_EXPRESSION_LEVELS_COL: add_gene_bgee_subgraph,
        Cons.DISGENET_DISEASE_COL: add_disgenet_gene_disease_subgraph,
        Cons.LITERATURE_DISEASE_COL: add_literature_gene_disease_subgraph,
        Cons.MINERVA_PATHWAY_COL: add_minerva_gene_pathway_subgraph,
        Cons.WIKIPATHWAYS: add_wikipathways_gene_pathway_subgraph,
        Cons.KEGG_PATHWAY_COL: add_kegg_gene_pathway_subgraph,
        Cons.OPENTARGETS_REACTOME_COL: add_opentargets_gene_reactome_pathway_subgraph,
        Cons.OPENTARGETS_GO_COL: add_opentargets_gene_go_subgraph,
        Cons.OPENTARGETS_GENE_COMPOUND_COL: add_opentargets_gene_compound_subgraph,
        Cons.MOLMEDB_PROTEIN_COMPOUND_COL: add_molmedb_gene_inhibitor_subgraph,
        Cons.PUBCHEM_COMPOUND_ASSAYS_COL: add_pubchem_assay_subgraph,
        Cons.WIKIPATHWAYS_PATHWAY_COL: add_wikipathways_gene_pathway_subgraph,
        Cons.WIKIPATHWAYS_MOLECULAR_COL: add_wikipathways_molecular_subgraph,
        Cons.ENSEMBL_HOMOLOG_COL: add_ensembl_homolog_subgraph,
        Cons.INTACT_INTERACT_COL: add_intact_interactions_subgraph,
        Cons.INTACT_COMPOUND_INTERACT_COL: add_intact_compound_interactions_subgraph,
        Cons.STRING_INTERACT_COL: add_stringdb_ppi_subgraph,
        Cons.AOPWIKI_GENE_COL: add_aopwiki_subgraph,
        # Cons.WIKIDATA_CC_COL: add_wikidata_gene_cc_subgraph,  # TODO: add this
        f"{Cons.GPROFILER}_wp": add_gprofiler_gene_wikipathway_subgraph,
        f"{Cons.GPROFILER}_hp": add_gprofiler_gene_phenotype_subgraph,
        f"{Cons.GPROFILER}_hpa": add_gprofiler_gene_hpa_subgraph,
        f"{Cons.GPROFILER}_kegg": add_gprofiler_gene_kegg_subgraph,
        f"{Cons.GPROFILER}_mirna": add_gprofiler_gene_mirna_subgraph,
        f"{Cons.GPROFILER}_reac": add_gprofiler_gene_reactome_subgraph,
        f"{Cons.GPROFILER}_tf": add_gprofiler_gene_transcription_factor_subgraph,
        f"{Cons.GPROFILER}_go:bp": add_gprofiler_gene_gobp_subgraph,
        f"{Cons.GPROFILER}_go:cc": add_gprofiler_gene_gocc_subgraph,
        f"{Cons.GPROFILER}_go:mf": add_gprofiler_gene_gomf_subgraph,
        Cons.MITOCART_PATHWAY_COL: add_mitocarta_gene_mito_subgraph,
        Cons.ITS_TARGET_COL: process_tf_target,
    }

    for _i, row in tqdm(combined_df.iterrows(), total=combined_df.shape[0], desc="Building graph"):
        if pd.isna(row[Cons.IDENTIFIER_COL]) or pd.isna(row[Cons.TARGET_COL]):
            continue

        gene_node_label = add_gene_node(g, row, dea_columns)
        if homolog_df_list is None:
            process_annotations(g, gene_node_label, row, func_dict)
        process_ppi(g, gene_node_label, row)

    if homolog_df_list is not None:
        process_homologs(g, combined_df, homolog_df_list, func_dict, dea_columns)

    # Process disease-compound relationships
    dnodes = {
        d["attr_dict"][Cons.EFO]: n
        for n, d in g.nodes(data=True)
        if d["attr_dict"][Cons.LABEL] == Cons.DISEASE_NODE_LABEL
        and d["attr_dict"][Cons.EFO] is not None
    }

    if disease_compound is not None:
        process_disease_compound(g, disease_compound, disease_nodes=dnodes)

    if pathway_compound is not None:
        process_kegg_pathway_compound(g, pathway_compound, combined_df)

    normalize_node_attributes(g)
    normalize_edge_attributes(g)

    return g


def _built_compound_based_graph(
    g: nx.MultiDiGraph,
    combined_df: pd.DataFrame,
    disease_compound=None,
    pathway_compound=None,
    homolog_df_list=None,
):
    """Build a gene-based graph."""
    combined_df = combined_df[(combined_df["target.source"] == "Ensembl")]

    dea_columns = [c for c in combined_df.columns if c.endswith("_dea")]

    func_dict = {
        Cons.MOLMEDB_COMPOUND_PROTEIN_COL: add_molmedb_compound_gene_subgraph,
    }  # type: ignore

    if homolog_df_list is not None:
        process_homologs(g, combined_df, homolog_df_list, func_dict, dea_columns)

    for _i, row in tqdm(combined_df.iterrows(), total=combined_df.shape[0], desc="Building graph"):
        if pd.isna(row["identifier"]) or pd.isna(row["target"]):
            continue
        gene_node_label = add_gene_node(g, row, dea_columns)
        process_annotations(g, gene_node_label, row, func_dict)
        process_ppi(g, gene_node_label, row)

    # Process disease-compound relationships
    dnodes = {
        d["attr_dict"][Cons.EFO]: n
        for n, d in g.nodes(data=True)
        if d["attr_dict"][Cons.LABEL] == Cons.DISEASE_NODE_LABEL
        and d["attr_dict"][Cons.EFO] is not None
    }
    if disease_compound is not None:
        process_disease_compound(g, disease_compound, disease_nodes=dnodes)

    if pathway_compound is not None:
        process_kegg_pathway_compound(g, pathway_compound, combined_df)

    normalize_node_attributes(g)
    normalize_edge_attributes(g)

    return g


def build_networkx_graph(
    combined_df: pd.DataFrame,
    disease_compound=None,
    pathway_compound=None,
    homolog_df_list=None,
) -> nx.MultiDiGraph:
    """Construct a NetWorkX graph from a Pandas DataFrame of genes and their multi-source annotations.

    :param combined_df: the input DataFrame to be converted into a graph.
    :param disease_compound: the input DataFrame containing disease-compound relationships.
    :param pathway_compound: the input DataFrame containing pathway-compound relationships from KEGG.
    :param homolog_df_list: a list of DataFrame generated by querying homologs.
    :returns: a NetworkX MultiDiGraph
    :raises ValueError: if the target type is not supported.
    """
    g = nx.MultiDiGraph()

    main_target_type = combined_df["target.source"].unique()[0]

    if main_target_type == Cons.ENSEMBL:
        return _built_gene_based_graph(
            g, combined_df, disease_compound, pathway_compound, homolog_df_list
        )
    elif main_target_type == Cons.PUBCHEM_COMPOUND:
        return _built_compound_based_graph(
            g, combined_df, disease_compound, pathway_compound, homolog_df_list
        )
    else:
        raise ValueError(f"Unsupported target type: {main_target_type}")
