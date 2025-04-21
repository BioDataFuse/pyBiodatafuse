# coding: utf-8

"""Python module to construct a NetworkX graph from the annotated data frame."""

import json
import os
import pickle
from logging import Logger
from typing import Any, Dict

import networkx as nx
import pandas as pd
from tqdm import tqdm

from pyBiodatafuse.constants import (
    BGEE,
    BGEE_ANATOMICAL_NODE_ATTRS,
    BGEE_ANATOMICAL_NODE_MAIN_LABEL,
    BGEE_EDGE_ATTRS,
    BGEE_GENE_ANATOMICAL_EDGE_LABEL,
    BGEE_GENE_EXPRESSION_LEVELS_COL,
    BRIDGEDB,
    COMPOUND_NODE_MAIN_LABEL,
    COMPOUND_SIDE_EFFECT_EDGE_ATTRS,
    COMPOUND_SIDE_EFFECT_EDGE_LABEL,
    DISEASE_NODE_MAIN_LABEL,
    DISGENET_DISEASE_COL,
    DISGENET_DISEASE_NODE_ATTRS,
    DISGENET_EDGE_ATTRS,
    ENSEMBL_HOMOLOG_COL,
    ENSEMBL_HOMOLOG_EDGE_ATTRS,
    ENSEMBL_HOMOLOG_EDGE_LABEL,
    ENSEMBL_HOMOLOG_MAIN_LABEL,
    ENSEMBL_HOMOLOG_NODE_ATTRS,
    GENE_DISEASE_EDGE_LABEL,
    GENE_GO_EDGE_ATTRS,
    GENE_GO_EDGE_LABEL,
    GENE_NODE_LABELS,
    GENE_PATHWAY_EDGE_ATTRS,
    GENE_PATHWAY_EDGE_LABEL,
    GO_BP_NODE_LABELS,
    GO_CC_NODE_LABELS,
    GO_MF_NODE_LABELS,
    GO_NODE_ATTRS,
    GO_NODE_MAIN_LABEL,
    HOMOLOG_NODE_LABELS,
    INTACT,
    INTACT_COMPOUND_INTERACT_COL,
    INTACT_COMPOUND_NODE_ATTRS,
    INTACT_COMPOUND_NODE_MAIN_LABEL,
    INTACT_COMPOUND_EDGE_LABEL,
    INTACT_COMPOUND_EDGE_ATTRS,
    INTACT_INTERACT_COL,
    INTACT_PPI_EDGE_ATTRS,
    INTACT_PPI_EDGE_LABEL,
    INTACT_PPI_EDGE_MAIN_LABEL,
    KEGG,
    KEGG_COL,
    KEGG_COMPOUND_COL,
    KEGG_COMPOUND_EDGE_ATTRS,
    KEGG_COMPOUND_EDGE_LABEL,
    KEGG_COMPOUND_NODE_ATTRS,
    KEGG_COMPOUND_NODE_MAIN_LABEL,
    LITERATURE_DISEASE_COL,
    LITERATURE_DISEASE_EDGE_ATTRS,
    LITERATURE_DISEASE_NODE_ATTRS,
    LITERATURE_NODE_MAIN_LABEL,
    MINERVA,
    MOLMEDB_COMPOUND_NODE_ATTRS,
    MOLMEDB_PROTEIN_COMPOUND_COL,
    MOLMEDB_PROTEIN_COMPOUND_EDGE_ATTRS,
    MOLMEDB_PROTEIN_COMPOUND_EDGE_LABEL,
    OPENTARGETS,
    OPENTARGETS_COMPOUND_NODE_ATTRS,
    OPENTARGETS_DISEASE_COMPOUND_COL,
    OPENTARGETS_DISEASE_COMPOUND_EDGE_ATTRS,
    OPENTARGETS_GENE_COMPOUND_COL,
    OPENTARGETS_GENE_COMPOUND_EDGE_ATTRS,
    OPENTARGETS_GO_COL,
    OPENTARGETS_REACTOME_COL,
    PATHWAY_NODE_ATTRS,
    PATHWAY_NODE_MAIN_LABEL,
    PUBCHEM_COMPOUND_ASSAYS_COL,
    PUBCHEM_COMPOUND_NODE_ATTRS,
    PUBCHEM_GENE_COMPOUND_EDGE_ATTRS,
    SIDE_EFFECT_NODE_ATTRS,
    SIDE_EFFECT_NODE_MAIN_LABEL,
    STRING_PPI_COL,
    STRING_PPI_EDGE_ATTRS,
    STRING_PPI_EDGE_LABEL,
    STRING_PPI_EDGE_MAIN_LABEL,
    WIKIPATHWAYS,
)

logger = Logger(__name__)


def load_dataframe_from_pickle(pickle_path: str) -> pd.DataFrame:
    """Load a previously annotated DataFrame from a pickle file.

    :param pickle_path: the path to a previously obtained annotation DataFrame dumped as a pickle file.
    :returns: a Pandas DataFrame.
    """
    with open(pickle_path, "rb") as rin:
        df = pickle.load(rin)
        df = df[(df["target.source"] == "Ensembl")]

    return df


def merge_node(g, node_label, node_attrs):
    """Merge the attr_dict of a newly added node to the graph on duplication, otherwise, add the new node.

    :param g: the graph to which the node will be added.
    :param node_label: node label.
    :param node_attrs: dictionary of node attributes.
    """
    if node_label not in g.nodes():
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
            g.add_node(node_label, attr_dict=node_attrs)


def add_gene_bgee_subgraph(g, gene_node_label, annot_list):
    """Construct part of the graph by linking the gene to a list of anatomical entities.

    :param g: the input graph to extend with new nodes and edges.
    :param gene_node_label: the gene node to be linked to annotation entities.
    :param annot_list: list of anatomical entities from Bgee with gene expression levels.
    :returns: a NetworkX MultiDiGraph
    """
    for annot in annot_list:
        if pd.isna(annot["anatomical_entity_name"]):
            continue
        annot_node_label = annot[BGEE_ANATOMICAL_NODE_MAIN_LABEL]
        annot_node_attrs = BGEE_ANATOMICAL_NODE_ATTRS.copy()
        annot_node_attrs["name"] = annot["anatomical_entity_name"]
        annot_node_attrs["id"] = annot["anatomical_entity_id"]

        g.add_node(annot_node_label, attr_dict=annot_node_attrs)

        edge_attrs = BGEE_EDGE_ATTRS.copy()
        if not pd.isna(annot["confidence_level_id"]):
            edge_attrs["confidence_level_name"] = annot["confidence_level_name"]
            edge_attrs["confidence_level_id"] = annot["confidence_level_id"]
        if not pd.isna(annot["expression_level"]):
            edge_attrs["expression_level"] = annot["expression_level"]
        if not pd.isna(annot["developmental_stage_name"]):
            edge_attrs["developmental_stage_name"] = annot["developmental_stage_name"]
        if not pd.isna(annot["developmental_stage_id"]):
            edge_attrs["developmental_stage_id"] = annot["developmental_stage_id"]

        edge_hash = hash(frozenset(edge_attrs.items()))
        edge_attrs["edge_hash"] = edge_hash
        edge_data = g.get_edge_data(gene_node_label, annot_node_label)
        edge_data = {} if edge_data is None else edge_data
        node_exists = [x for x, y in edge_data.items() if y["attr_dict"]["edge_hash"] == edge_hash]

        if len(node_exists) == 0:
            g.add_edge(
                gene_node_label,
                annot_node_label,
                label=BGEE_GENE_ANATOMICAL_EDGE_LABEL,
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
    for annot in annot_list:
        if not pd.isna(annot["disease_name"]):
            annot_node_label = annot[DISEASE_NODE_MAIN_LABEL]
            annot_node_attrs = DISGENET_DISEASE_NODE_ATTRS.copy()
            annot_node_attrs["name"] = annot["disease_name"]
            annot_node_attrs["id"] = annot["UMLS"]

            if not pd.isna(annot["HPO"]):
                annot_node_attrs["HPO"] = annot["HPO"]
            if not pd.isna(annot["NCI"]):
                annot_node_attrs["NCI"] = annot["NCI"]
            if not pd.isna(annot["OMIM"]):
                annot_node_attrs["OMIM"] = annot["OMIM"]
            if not pd.isna(annot["MONDO"]):
                annot_node_attrs["MONDO"] = annot["MONDO"]
            if not pd.isna(annot["ORDO"]):
                annot_node_attrs["ORDO"] = annot["ORDO"]
            if not pd.isna(annot["EFO"]):
                annot_node_attrs["EFO"] = annot["EFO"]
            if not pd.isna(annot["DO"]):
                annot_node_attrs["DO"] = annot["DO"]
            if not pd.isna(annot["MESH"]):
                annot_node_attrs["MESH"] = annot["MESH"]
            if not pd.isna(annot["UMLS"]):
                annot_node_attrs["UMLS"] = annot["UMLS"]
            if not pd.isna(annot["disease_type"]):
                annot_node_attrs["disease_type"] = annot["disease_type"]
            if not pd.isna(annot["disease_umlscui"]):
                annot_node_attrs["disease_umlscui"] = annot["disease_umlscui"]

            g.add_node(annot_node_label, attr_dict=annot_node_attrs)

            edge_attrs = DISGENET_EDGE_ATTRS.copy()
            edge_attrs["score"] = annot["score"]

            if not pd.isna(annot["ei"]):
                edge_attrs["ei"] = annot["ei"]
            if not pd.isna(annot["el"]):
                edge_attrs["el"] = annot["el"]

            edge_hash = hash(frozenset(edge_attrs.items()))
            edge_attrs["edge_hash"] = edge_hash
            edge_data = g.get_edge_data(gene_node_label, annot_node_label)
            edge_data = {} if edge_data is None else edge_data
            node_exists = [
                x for x, y in edge_data.items() if y["attr_dict"]["edge_hash"] == edge_hash
            ]

            if len(node_exists) == 0:
                g.add_edge(
                    gene_node_label,
                    annot_node_label,
                    label=GENE_DISEASE_EDGE_LABEL,
                    attr_dict=edge_attrs,
                )

    return g

def add_intact_ppi_subgraph(g, gene_node_label, annot_list):
    """Construct part of the graph by linking the gene to genes via IntAct PPI, filtering based on unique interaction IDs.

    :param g: the input graph to extend with new nodes and edges.
    :param gene_node_label: the gene node to be linked to other gene nodes.
    :param annot_list: list of protein-protein interactions from IntAct.
    :returns: a NetworkX MultiDiGraph
    """
    seen_interaction_ids = set()

    for interaction in annot_list:   
        partner_gene = interaction.get("intact_link_to")
        if not partner_gene or pd.isna(partner_gene):
            continue  

        interaction_id = interaction.get("binary_interaction_id", None)
        if not interaction_id or interaction_id in seen_interaction_ids:
            continue

        edge_attrs = INTACT_PPI_EDGE_ATTRS.copy()

        conf_vals = interaction.get("confidence_values")
        if isinstance(conf_vals, list):
            edge_attrs["confidence_values"] = ",".join(map(str, conf_vals))
        else:
            edge_attrs["confidence_values"] = str(conf_vals) if conf_vals else ""

        if "interaction_id" in interaction:
            edge_attrs["interaction_id"] = interaction["interaction_id"]
        if "interactor_id_A" in interaction:
            edge_attrs["interactor_id_A"] = interaction["interactor_id_A"]
        if "interactor_id_B" in interaction:
            edge_attrs["interactor_id_B"] = interaction["interactor_id_B"]
        if "binary_interaction_id" in interaction:
            edge_attrs["binary_interaction_id"] = interaction["binary_interaction_id"]
        if "intact_score" in interaction:
            edge_attrs["intact_score"] = interaction["intact_score"]
        if "biological_role_A" in interaction:
            edge_attrs["biological_role_A"] = interaction["biological_role_A"]
        if "biological_role_B" in interaction:
            edge_attrs["biological_role_B"] = interaction["biological_role_B"]
        if "type" in interaction:
            edge_attrs["type"] = interaction["type"]
        if "stoichiometry_A" in interaction:
            edge_attrs["stoichiometry_A"] = interaction["stoichiometry_A"]
        if "stoichiometry_B" in interaction:
            edge_attrs["stoichiometry_B"] = interaction["stoichiometry_B"]
        if "detection_method" in interaction:
            edge_attrs["detection_method"] = interaction["detection_method"]
        if "detection_method_id" in interaction:
            edge_attrs["detection_method_id"] = interaction["detection_method_id"]
        if "host_organism" in interaction:
            edge_attrs["host_organism"] = interaction["host_organism"]
        if "interactor_A_name" in interaction:
            edge_attrs["interactor_A_name"] = interaction["interactor_A_name"]
        if "interactor_B_name" in interaction:
            edge_attrs["interactor_B_name"] = interaction["interactor_B_name"]
        if "interactor_A_species" in interaction:
            edge_attrs["interactor_A_species"] = interaction["interactor_A_species"]
        if "interactor_B_species" in interaction:
            edge_attrs["interactor_B_species"] = interaction["interactor_B_species"]
        if "molecule_A" in interaction:
            edge_attrs["molecule_A"] = interaction["molecule_A"]
        if "molecule_B" in interaction:
            edge_attrs["molecule_B"] = interaction["molecule_B"]
        if "pubmed_publication_id" in interaction:
            edge_attrs["pubmed_publication_id"] = interaction["pubmed_publication_id"]
        if "id_A" in interaction:
            edge_attrs["id_A"] = interaction["id_A"]
        if "id_B" in interaction:
            edge_attrs["id_B"] = interaction["id_B"]

        for k, v in edge_attrs.items():
            if isinstance(v, list):
                edge_attrs[k] = ",".join(map(str, v))

        edge_hash = hash(frozenset(edge_attrs.items()))
        edge_attrs["edge_hash"] = edge_hash

        g.add_edge(
            gene_node_label,
            partner_gene,
            label=INTACT_PPI_EDGE_LABEL,
            attr_dict=edge_attrs,
        )

        edge_data = g.get_edge_data(gene_node_label, partner_gene)
        
        seen_interaction_ids.add(interaction_id)

    return g


def add_intact_compound_subgraph(g, gene_node_label, annot_list):
    """Construct part of the graph by linking the gene to compounds via IntAct compound interactions.

    :param g: the input graph to extend with new nodes and edges.
    :param gene_node_label: the gene node to be linked to compounds.
    :param annot_list: list of compound-gene interaction dicts from IntAct.
    :returns: a NetworkX MultiDiGraph
    """
    seen_edge_signatures = set()

    for interaction in annot_list:
        if str(interaction.get("id_A", "")).startswith("CHEBI:"):
            compound_id = interaction.get("id_A")
            compound_name = interaction.get("interactor_A_name", "")
            compound_species = interaction.get("interactor_A_species", "")
            molecule = interaction.get("molecule_A", "")
        elif str(interaction.get("id_B", "")).startswith("CHEBI:"):
            compound_id = interaction.get("id_B")
            compound_name = interaction.get("interactor_B_name", "")
            compound_species = interaction.get("interactor_B_species", "")
            molecule = interaction.get("molecule_B", "")
        else:
            continue

        if not compound_id or pd.isna(compound_id):
            continue

        annot_node_label = compound_id
        annot_node_attrs = INTACT_COMPOUND_NODE_ATTRS.copy()
        annot_node_attrs["id"] = compound_id
        annot_node_attrs["label"] = compound_name
        annot_node_attrs["species"] = compound_species
        annot_node_attrs["molecule"] = molecule

        merge_node(g, annot_node_label, annot_node_attrs)

        edge_attrs = INTACT_PPI_EDGE_ATTRS.copy()

        conf_vals = interaction.get("confidence_values")
        edge_attrs["confidence_values"] = ",".join(map(str, conf_vals)) if isinstance(conf_vals, list) else str(conf_vals or "")

        for key, value in interaction.items():
            if isinstance(value, list):
                edge_attrs[key] = ",".join(map(str, value))
            else:
                edge_attrs[key] = value

        unstable_keys = {"interaction_id", "interactor_id_A", "interactor_id_B", "binary_interaction_id"}
        signature_attrs = {k: v for k, v in edge_attrs.items() if k not in unstable_keys}
        edge_signature = frozenset(signature_attrs.items())

        if edge_signature in seen_edge_signatures:
            continue

        edge_hash = hash(edge_signature)
        edge_attrs["edge_hash"] = edge_hash

        g.add_edge(
            annot_node_label,
            gene_node_label,
            label="compound-gene",
            attr_dict=edge_attrs,
        )

        seen_edge_signatures.add(edge_signature)

    return g


def add_literature_gene_disease_subgraph(g, gene_node_label, annot_list):
    """Construct part of the graph by linking the gene to diseases form literature.

    :param g: the input graph to extend with new nodes and edges.
    :param gene_node_label: the gene node to be linked to diseases.
    :param annot_list: list of diseases from DisGeNET.
    :returns: a NetworkX MultiDiGraph
    """
    for annot in annot_list:
        if pd.isna(annot["disease_name"]):
            continue
        annot_node_label = annot[LITERATURE_NODE_MAIN_LABEL]
        annot_node_attrs = LITERATURE_DISEASE_NODE_ATTRS.copy()
        annot_node_attrs["datasource"] = annot["source"]
        annot_node_attrs["name"] = annot["disease_name"]
        annot_node_attrs["id"] = annot["UMLS"]
        annot_node_attrs["UMLS"] = annot["UMLS"]
        annot_node_attrs["MONDO"] = annot["MONDO"]

        g.add_node(annot_node_label, attr_dict=annot_node_attrs)

        edge_attrs = LITERATURE_DISEASE_EDGE_ATTRS.copy()
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
                label=GENE_DISEASE_EDGE_LABEL,
                attr_dict=edge_attrs,
            )

    return g


# TODO: The disease annotations are not curated and will be used again when the OpenTarget annotation improves.
# def add_opentargets_gene_disease_subgraph(g, gene_node_label, annot_list):  # TODO: should be updated
#     """Construct part of the graph by linking the gene to a list of annotation entities (disease, compound ..etc).

#     :param g: the input graph to extend with new nodes and edges.
#     :param gene_node_label: the gene node to be linked to annotation entities.
#     :param annot_list: list of annotations from a specific source (e.g. DisGeNET, WikiPathways ..etc).
#     :returns: a NetworkX MultiDiGraph
#     """
#     for annot in annot_list:
#         if not pd.isna(annot["disease_name"]):
#             annot_node_label = annot[DISEASE_NODE_MAIN_LABEL]
#             annot_node_attrs = OPENTARGETS_DISEASE_NODE_ATTRS.copy()
#             annot_node_attrs["name"] = annot["disease_name"]
#             annot_node_attrs["id"] = annot["disease_id"]
#             annot_node_attrs["therapeutic_areas"] = annot["therapeutic_areas"]

#             # g.add_node(annot_node_label, attr_dict=annot_node_attrs)
#             merge_node(g, annot_node_label, annot_node_attrs)

#             edge_attrs = OPENTARGETS_DISEASE_EDGE_ATTRS

#             edge_hash = hash(frozenset(edge_attrs.items()))
#             edge_attrs["edge_hash"] = edge_hash
#             edge_data = g.get_edge_data(gene_node_label, annot_node_label)
#             edge_data = {} if edge_data is None else edge_data
#             node_exists = [
#                 x for x, y in edge_data.items() if y["attr_dict"]["edge_hash"] == edge_hash
#             ]

#             if len(node_exists) == 0:
#                 g.add_edge(
#                     gene_node_label,
#                     annot_node_label,
#                     label=OPENTARGETS_GENE_DISEASE_EDGE_LABEL,
#                     attr_dict=edge_attrs,
#                 )

#     return g


def add_minerva_gene_pathway_subgraph(g, gene_node_label, annot_list):
    """Construct part of the graph by linking the gene to MINERVA pathways.

    :param g: the input graph to extend with new nodes and edges.
    :param gene_node_label: the gene node to be linked to MINERVA pathways.
    :param annot_list: list of MINERVA pathways from MINERVA.
    :returns: a NetworkX MultiDiGraph
    """
    for annot in annot_list:
        if pd.isna(annot.get("pathway_label")):
            continue

        annot_node_label = annot[PATHWAY_NODE_MAIN_LABEL]
        annot_node_attrs = PATHWAY_NODE_ATTRS.copy()
        annot_node_attrs["datasource"] = MINERVA
        annot_node_attrs["name"] = annot["pathway_label"]
        annot_node_attrs["id"] = annot["pathway_id"]
        annot_node_attrs["gene_count"] = annot["pathway_gene_count"]

        g.add_node(annot_node_label, attr_dict=annot_node_attrs)

        edge_attrs = GENE_PATHWAY_EDGE_ATTRS.copy()
        edge_attrs["datasource"] = MINERVA

        edge_hash = hash(frozenset(edge_attrs.items()))
        edge_attrs["edge_hash"] = edge_hash
        edge_data = g.get_edge_data(gene_node_label, annot_node_label)
        edge_data = {} if edge_data is None else edge_data
        node_exists = [x for x, y in edge_data.items() if y["attr_dict"]["edge_hash"] == edge_hash]

        if len(node_exists) == 0:
            g.add_edge(
                gene_node_label,
                annot_node_label,
                label=GENE_PATHWAY_EDGE_LABEL,
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
    for annot in annot_list:
        if pd.isna(annot["pathway_label"]):
            continue

        annot_node_label = annot[PATHWAY_NODE_MAIN_LABEL]
        annot_node_attrs = PATHWAY_NODE_ATTRS.copy()
        annot_node_attrs["datasource"] = WIKIPATHWAYS
        annot_node_attrs["name"] = annot["pathway_label"]
        annot_node_attrs["id"] = annot["pathway_id"]
        annot_node_attrs["gene_count"] = annot["pathway_gene_count"]

        g.add_node(annot_node_label, attr_dict=annot_node_attrs)

        edge_attrs = GENE_PATHWAY_EDGE_ATTRS.copy()
        edge_attrs["datasource"] = WIKIPATHWAYS

        edge_hash = hash(frozenset(edge_attrs.items()))
        edge_attrs["edge_hash"] = edge_hash
        edge_data = g.get_edge_data(gene_node_label, annot_node_label)
        edge_data = {} if edge_data is None else edge_data
        node_exists = [x for x, y in edge_data.items() if y["attr_dict"]["edge_hash"] == edge_hash]

        if len(node_exists) == 0:
            g.add_edge(
                gene_node_label,
                annot_node_label,
                label=GENE_PATHWAY_EDGE_LABEL,
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
    for annot in annot_list:
        if pd.isna(annot.get("pathway_label")):
            continue

        annot_node_label = annot[PATHWAY_NODE_MAIN_LABEL]
        annot_node_attrs = PATHWAY_NODE_ATTRS.copy()
        annot_node_attrs["datasource"] = KEGG
        annot_node_attrs["name"] = annot["pathway_label"]
        annot_node_attrs["id"] = annot["pathway_id"]
        annot_node_attrs["gene_count"] = annot["gene_count"]

        # g.add_node(annot_node_label, attr_dict=annot_node_attrs)
        merge_node(g, annot_node_label, annot_node_attrs)

        edge_attrs = GENE_PATHWAY_EDGE_ATTRS.copy()
        edge_attrs["datasource"] = KEGG

        edge_hash = hash(frozenset(edge_attrs.items()))
        edge_attrs["edge_hash"] = edge_hash
        edge_data = g.get_edge_data(gene_node_label, annot_node_label)
        edge_data = {} if edge_data is None else edge_data
        node_exists = [x for x, y in edge_data.items() if y["attr_dict"]["edge_hash"] == edge_hash]

        if len(node_exists) == 0:
            g.add_edge(
                gene_node_label,
                annot_node_label,
                label=GENE_PATHWAY_EDGE_LABEL,
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
    for compound in compounds_list:
        if pd.isna(compound.get("name")):
            continue

        annot_node_label = compound["KEGG_identifier"]
        annot_node_attrs = KEGG_COMPOUND_NODE_ATTRS.copy()
        annot_node_attrs["id"] = compound["KEGG_identifier"]
        annot_node_attrs["label"] = compound["name"]

        merge_node(g, annot_node_label, annot_node_attrs)

        for _, path_row in combined_df.iterrows():
            pathways = path_row.get("KEGG_pathways", [])
            if isinstance(pathways, list) and pathways:
                for pathway in pathways:
                    if pathway_node_label == pathway.get("pathway_id", ""):
                        if "compounds" in pathway:
                            pathway_compounds = [
                                comp["KEGG_identifier"] for comp in pathway["compounds"]
                            ]
                            if compound["KEGG_identifier"] in pathway_compounds:
                                edge_attrs = KEGG_COMPOUND_EDGE_ATTRS.copy()
                                edge_hash = hash(frozenset(edge_attrs.items()))
                                edge_attrs["edge_hash"] = edge_hash
                                edge_data = g.get_edge_data(pathway_node_label, annot_node_label)
                                edge_data = {} if edge_data is None else edge_data
                                node_exists = [
                                    x
                                    for x, y in edge_data.items()
                                    if "attr_dict" in y
                                    and y["attr_dict"].get("edge_hash") == edge_hash
                                ]

                                if len(node_exists) == 0:
                                    g.add_edge(
                                        pathway_node_label,
                                        annot_node_label,
                                        label=KEGG_COMPOUND_EDGE_LABEL,
                                        attr_dict=edge_attrs,
                                    )

    return g


def process_kegg_pathway_compound(g, kegg_pathway_compound, combined_df):
    """Process pathway-compound relationships from KEGG and add them to the graph.

    :param g: the input graph to extend with pathway-compound relationships.
    :param kegg_pathway_compound: DataFrame containing pathway-compound relationships.
    :param combined_df: DataFrame containing KEGG pathway data.
    """
    for _, row in kegg_pathway_compound.iterrows():
        compound_info = row[KEGG_COMPOUND_COL]

        if isinstance(compound_info, dict):
            compounds_list = [compound_info]
        elif isinstance(compound_info, list):
            compounds_list = compound_info
        else:
            compounds_list = []

        for compound in compounds_list:
            compound_id = compound["KEGG_identifier"]

            for _, pathway_row in combined_df.iterrows():
                pathway_data = pathway_row["KEGG_pathways"]

                if isinstance(pathway_data, list):
                    for pathway in pathway_data:
                        pathway_id = pathway.get("pathway_id")
                        pathway_compounds = pathway.get("compounds", [])

                        if any(c.get("KEGG_identifier") == compound_id for c in pathway_compounds):
                            add_kegg_compounds_subgraph(g, pathway_id, compounds_list, combined_df)


def add_opentargets_gene_reactome_pathway_subgraph(g, gene_node_label, annot_list):
    """Construct part of the graph by linking the gene to Reactome pathways.

    :param g: the input graph to extend with new nodes and edges.
    :param gene_node_label: the gene node to be linked to Reactome pathways.
    :param annot_list: list of Reactome pathways from OpenTargets.
    :returns: a NetworkX MultiDiGraph
    """
    for annot in annot_list:
        if pd.isna(annot["pathway_id"]):
            continue

        annot_node_label = annot[PATHWAY_NODE_MAIN_LABEL]
        annot_node_attrs = PATHWAY_NODE_ATTRS.copy()
        annot_node_attrs["datasource"] = OPENTARGETS
        annot_node_attrs["name"] = annot["pathway_label"]
        annot_node_attrs["id"] = annot["pathway_id"]

        g.add_node(annot_node_label, attr_dict=annot_node_attrs)

        edge_attrs = GENE_PATHWAY_EDGE_ATTRS.copy()
        edge_attrs["datasource"] = OPENTARGETS

        edge_hash = hash(frozenset(edge_attrs.items()))
        edge_attrs["edge_hash"] = edge_hash
        edge_data = g.get_edge_data(gene_node_label, annot_node_label)
        edge_data = {} if edge_data is None else edge_data
        node_exists = [x for x, y in edge_data.items() if y["attr_dict"]["edge_hash"] == edge_hash]

        if len(node_exists) == 0:
            g.add_edge(
                gene_node_label,
                annot_node_label,
                label=GENE_PATHWAY_EDGE_LABEL,
                attr_dict=edge_attrs,
            )

    return g


def add_opentargets_gene_go_subgraph(g, gene_node_label, annot_list):
    """Construct part of the graph by linking the gene to gene ontologies.

    :param g: the input graph to extend with new nodes and edges.
    :param gene_node_label: the gene node to be linked to gene ontologies.
    :param annot_list: list of gene ontologies from OpenTargets.
    :returns: a NetworkX MultiDiGraph
    """
    for annot in annot_list:
        if pd.isna(annot["go_id"]):
            continue

        annot_node_label = annot[GO_NODE_MAIN_LABEL]
        annot_node_attrs = GO_NODE_ATTRS.copy()
        annot_node_attrs["name"] = annot["go_name"]
        annot_node_attrs["id"] = annot["go_id"]
        if annot["go_type"] == "P":
            annot_node_attrs["labels"] = GO_BP_NODE_LABELS
        elif annot["go_type"] == "F":
            annot_node_attrs["labels"] = GO_MF_NODE_LABELS
        elif annot["go_type"] == "C":
            annot_node_attrs["labels"] = GO_CC_NODE_LABELS

        g.add_node(annot_node_label, attr_dict=annot_node_attrs)

        edge_attrs = GENE_GO_EDGE_ATTRS.copy()

        edge_hash = hash(frozenset(edge_attrs.items()))
        edge_attrs["edge_hash"] = edge_hash
        edge_data = g.get_edge_data(gene_node_label, annot_node_label)
        edge_data = {} if edge_data is None else edge_data
        node_exists = [x for x, y in edge_data.items() if y["attr_dict"]["edge_hash"] == edge_hash]

        if len(node_exists) == 0:
            g.add_edge(
                gene_node_label,
                annot_node_label,
                label=GENE_GO_EDGE_LABEL,
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
    if isinstance(side_effects_list, list):
        for effect in side_effects_list:
            if pd.isna(effect["name"]):
                continue

            effect_node_label = effect["name"]
            effect_node_attrs = SIDE_EFFECT_NODE_ATTRS.copy()
            effect_node_attrs["name"] = effect["name"]

            g.add_node(effect_node_label, attr_dict=effect_node_attrs)

            edge_attrs = COMPOUND_SIDE_EFFECT_EDGE_ATTRS.copy()
            edge_hash = hash(frozenset(edge_attrs.items()))
            edge_attrs["edge_hash"] = edge_hash
            edge_data = g.get_edge_data(compound_node_label, effect_node_label)
            edge_data = {} if edge_data is None else edge_data
            node_exists = [
                x
                for x, y in edge_data.items()
                if "attr_dict" in y and y["attr_dict"].get("edge_hash") == edge_hash
            ]

            if len(node_exists) == 0:
                g.add_edge(
                    compound_node_label,
                    effect_node_label,
                    label=COMPOUND_SIDE_EFFECT_EDGE_LABEL,
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
    for annot in annot_list:
        if pd.isna(annot["relation"]):
            continue

        if not pd.isna(annot[COMPOUND_NODE_MAIN_LABEL]):
            annot_node_label = annot[COMPOUND_NODE_MAIN_LABEL]
        else:
            annot_node_label = annot["chembl_id"]
        annot_node_attrs = OPENTARGETS_COMPOUND_NODE_ATTRS.copy()
        annot_node_attrs["name"] = annot["compound_name"]
        if not pd.isna(annot[COMPOUND_NODE_MAIN_LABEL]):
            annot_node_attrs["id"] = annot[COMPOUND_NODE_MAIN_LABEL]
        else:
            annot_node_attrs["id"] = annot["chembl_id"]
        annot_node_attrs["chembl_id"] = annot["chembl_id"]
        if not pd.isna(annot["drugbank_id"]):
            annot_node_attrs["drugbank_id"] = annot["drugbank_id"]
        if not pd.isna(annot["compound_cid"]):
            annot_node_attrs["compound_cid"] = annot["compound_cid"]
        if not pd.isna(annot["clincal_trial_phase"]):
            annot_node_attrs["clincal_trial_phase"] = annot["clincal_trial_phase"]
        annot_node_attrs["is_approved"] = annot["is_approved"]
        if not pd.isna(annot["adverse_effect_count"]):
            annot_node_attrs["adverse_effect_count"] = annot["adverse_effect_count"]

        # g.add_node(annot_node_label, attr_dict=annot_node_attrs)
        merge_node(g, annot_node_label, annot_node_attrs)

        edge_attrs = OPENTARGETS_GENE_COMPOUND_EDGE_ATTRS.copy()
        edge_attrs["label"] = annot["relation"]
        edge_hash = hash(frozenset(edge_attrs.items()))
        edge_attrs["edge_hash"] = edge_hash
        edge_data = g.get_edge_data(annot_node_label, gene_node_label)
        edge_data = {} if edge_data is None else edge_data
        node_exists = [x for x, y in edge_data.items() if y["attr_dict"]["edge_hash"] == edge_hash]

        if len(node_exists) == 0:
            g.add_edge(
                annot_node_label,
                gene_node_label,
                label=annot["relation"],
                attr_dict=edge_attrs,
            )

        # Add side effects
        if annot["adverse_effect"]:
            add_opentargets_compound_side_effect_subgraph(
                g, annot_node_label, annot[SIDE_EFFECT_NODE_MAIN_LABEL]
            )

    return g


def add_molmedb_gene_inhibitor_subgraph(g, gene_node_label, annot_list):
    """Construct part of the graph by linking the gene to its inhibitors.

    :param g: the input graph to extend with new nodes and edges.
    :param gene_node_label: the gene node to be linked to its inhibitors.
    :param annot_list: list of gene inhibitors from MolMeDB.
    :returns: a NetworkX MultiDiGraph
    """
    for annot in annot_list:
        if pd.isna(annot["compound_name"]):
            continue

        if not pd.isna(annot[COMPOUND_NODE_MAIN_LABEL]):
            annot_node_label = annot[COMPOUND_NODE_MAIN_LABEL]
        else:
            annot_node_label = annot["molmedb_id"]

        annot_node_attrs = MOLMEDB_COMPOUND_NODE_ATTRS.copy()
        annot_node_attrs["name"] = annot["compound_name"]

        if not pd.isna(annot[COMPOUND_NODE_MAIN_LABEL]):
            annot_node_attrs["id"] = annot[COMPOUND_NODE_MAIN_LABEL]
        else:
            annot_node_attrs["id"] = annot["molmedb_id"]

        annot_node_attrs["molmedb_id"] = annot["molmedb_id"]
        if not pd.isna(annot["inchikey"]):
            annot_node_attrs["inchikey"] = annot["inchikey"]
        if not pd.isna(annot["smiles"]):
            annot_node_attrs["smiles"] = annot["smiles"]
        if not pd.isna(annot["compound_cid"]):
            annot_node_attrs["compound_cid"] = annot["compound_cid"]
        if not pd.isna(annot["chebi_id"]):
            annot_node_attrs["chebi_id"] = annot["chebi_id"]
        if not pd.isna(annot["drugbank_id"]):
            annot_node_attrs["drugbank_id"] = annot["drugbank_id"]
        # if not pd.isna(annot["pdb_ligand_id"]):
        #     annot_node_attrs["pdb_ligand_id"] = annot["pdb_ligand_id"]
        if not pd.isna(annot["source_pmid"]):
            annot_node_attrs["source_pmid"] = annot["source_pmid"]
        if not pd.isna(annot["uniprot_trembl_id"]):
            annot_node_attrs["uniprot_trembl_id"] = annot["uniprot_trembl_id"]

        merge_node(g, annot_node_label, annot_node_attrs)

        edge_attrs = MOLMEDB_PROTEIN_COMPOUND_EDGE_ATTRS.copy()

        edge_hash = hash(frozenset(edge_attrs.items()))
        edge_attrs["edge_hash"] = edge_hash
        edge_data = g.get_edge_data(gene_node_label, annot_node_label)
        edge_data = {} if edge_data is None else edge_data
        node_exists = [x for x, y in edge_data.items() if y["attr_dict"]["edge_hash"] == edge_hash]

        if len(node_exists) == 0:
            g.add_edge(
                annot_node_label,
                gene_node_label,
                label=MOLMEDB_PROTEIN_COMPOUND_EDGE_LABEL,
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
    for annot in annot_list:
        if pd.isna(annot["pubchem_assay_id"]):
            continue

        annot_node_label = annot[COMPOUND_NODE_MAIN_LABEL]
        annot_node_attrs = PUBCHEM_COMPOUND_NODE_ATTRS.copy()
        annot_node_attrs["name"] = annot["compound_name"]
        annot_node_attrs["id"] = annot["compound_cid"]
        annot_node_attrs["inchi"] = annot["inchi"]
        if not pd.isna(annot["smiles"]):
            annot_node_attrs["smiles"] = annot["smiles"]

        # g.add_node(annot_node_label, attr_dict=annot_node_attrs)
        merge_node(g, annot_node_label, annot_node_attrs)

        edge_attrs = PUBCHEM_GENE_COMPOUND_EDGE_ATTRS.copy()
        edge_attrs["assay_type"] = annot["assay_type"]
        edge_attrs["pubchem_assay_id"] = annot["pubchem_assay_id"]
        edge_attrs["outcome"] = annot["outcome"]
        edge_attrs["label"] = annot["outcome"]

        edge_hash = hash(frozenset(edge_attrs.items()))
        edge_attrs["edge_hash"] = edge_hash
        edge_data = g.get_edge_data(gene_node_label, annot_node_label)
        edge_data = {} if edge_data is None else edge_data
        node_exists = [x for x, y in edge_data.items() if y["attr_dict"]["edge_hash"] == edge_hash]

        if len(node_exists) == 0:
            g.add_edge(
                annot_node_label,
                gene_node_label,
                label=annot["outcome"],
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
    for ppi in annot_list:
        edge_attrs = STRING_PPI_EDGE_ATTRS.copy()
        edge_attrs["score"] = ppi["score"]

        edge_hash = hash(frozenset(edge_attrs.items()))
        edge_attrs["edge_hash"] = edge_hash
        edge_data = g.get_edge_data(gene_node_label, ppi[STRING_PPI_EDGE_MAIN_LABEL])

        edge_data = {} if edge_data is None else edge_data
        node_exists = [x for x, y in edge_data.items() if y["attr_dict"]["edge_hash"] == edge_hash]
        if len(node_exists) == 0 and not pd.isna(ppi[STRING_PPI_EDGE_MAIN_LABEL]):
            g.add_edge(
                gene_node_label,
                ppi[STRING_PPI_EDGE_MAIN_LABEL],
                label=STRING_PPI_EDGE_LABEL,
                attr_dict=edge_attrs,
            )
            g.add_edge(
                ppi[STRING_PPI_EDGE_MAIN_LABEL],
                gene_node_label,
                label=STRING_PPI_EDGE_LABEL,
                attr_dict=edge_attrs,
            )

    return g


def add_opentargets_disease_compound_subgraph(g, disease_node_label, annot_list):
    """Construct part of the graph by linking the disease to compounds.

    :param g: the input graph to extend with new nodes and edges.
    :param disease_node_label: the disease node to be linked to compounds.
    :param annot_list: list of compounds from OpenTargets.
    :returns: a NetworkX MultiDiGraph
    """
    for annot in annot_list:
        if pd.isna(annot["relation"]):
            continue

        if pd.isna(annot[COMPOUND_NODE_MAIN_LABEL]):
            annot_node_label = annot["chembl_id"]
        else:
            annot_node_label = annot[COMPOUND_NODE_MAIN_LABEL]

        # create compound node and merge with existing node if it exists
        annot_node_attrs = OPENTARGETS_COMPOUND_NODE_ATTRS.copy()
        annot_node_attrs["name"] = annot["compound_name"]
        annot_node_attrs["id"] = annot_node_label
        annot_node_attrs["chembl_id"] = annot["chembl_id"]
        if not pd.isna(annot["drugbank_id"]):
            annot_node_attrs["drugbank_id"] = annot["drugbank_id"]
        if not pd.isna(annot["compound_cid"]):
            annot_node_attrs["compound_cid"] = annot["compound_cid"]
        if not pd.isna(annot["clincal_trial_phase"]):
            annot_node_attrs["clincal_trial_phase"] = annot["clincal_trial_phase"]
        annot_node_attrs["is_approved"] = annot["is_approved"]
        if not pd.isna(annot["adverse_effect_count"]):
            annot_node_attrs["adverse_effect_count"] = annot["adverse_effect_count"]

        merge_node(g, annot_node_label, annot_node_attrs)

        edge_attrs = OPENTARGETS_DISEASE_COMPOUND_EDGE_ATTRS.copy()
        edge_attrs["label"] = annot["relation"]
        edge_hash = hash(frozenset(edge_attrs.items()))
        edge_attrs["edge_hash"] = edge_hash

        edge_data = g.get_edge_data(annot_node_label, disease_node_label)
        edge_data = {} if edge_data is None else edge_data
        node_exists = [x for x, y in edge_data.items() if y["attr_dict"]["edge_hash"] == edge_hash]

        if len(node_exists) == 0:
            g.add_edge(
                annot_node_label,
                disease_node_label,
                label=annot["relation"],
                attr_dict=edge_attrs,
            )

        # Add side effects
        if annot["adverse_effect"]:
            add_opentargets_compound_side_effect_subgraph(
                g, annot_node_label, annot[SIDE_EFFECT_NODE_MAIN_LABEL]
            )

    return g


def add_ensembl_homolog_subgraph(g, gene_node_label, annot_list):
    """Construct part of the graph by linking the gene to genes.

    :param g: the input graph to extend with new nodes and edges.
    :param gene_node_label: the gene node to be linked to other genes entities.
    :param annot_list: list of homologs from Ensembl.
    :returns: a NetworkX MultiDiGraph
    """
    for hl in annot_list:
        edge_attrs = ENSEMBL_HOMOLOG_EDGE_ATTRS.copy()

        edge_hash = hash(frozenset(edge_attrs.items()))
        edge_attrs["edge_hash"] = edge_hash
        edge_data = g.get_edge_data(gene_node_label, hl[ENSEMBL_HOMOLOG_MAIN_LABEL])

        edge_data = {} if edge_data is None else edge_data
        node_exists = [x for x, y in edge_data.items() if y["attr_dict"]["edge_hash"] == edge_hash]
        if len(node_exists) == 0 and not pd.isna(hl[ENSEMBL_HOMOLOG_MAIN_LABEL]):
            g.add_edge(
                gene_node_label,
                hl[ENSEMBL_HOMOLOG_MAIN_LABEL],
                label=ENSEMBL_HOMOLOG_EDGE_LABEL,
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
    gene_node_label = row["identifier"]
    gene_node_attrs = {
        "datasource": BRIDGEDB,
        "name": f"{row['identifier.source']}:{row['identifier']}",
        "id": row["target"],
        "labels": GENE_NODE_LABELS,
        row["target.source"]: row["target"],
    }

    for c in dea_columns:
        gene_node_attrs[c[:-4]] = row[c]

    g.add_node(gene_node_label, attr_dict=gene_node_attrs)
    return gene_node_label


def process_annotations(g, gene_node_label, row, func_dict):
    """Process the annotations for gene node from each row of the combined_df to the graph.

    :param g: the input graph to extend with gene nodes.
    :param gene_node_label: the gene node to be linked to annotation entities.
    :param row: row in the combined DataFrame.
    :param func_dict: dictionary of subgraph function.
    """
    for annot_key in func_dict:
        if annot_key in row:
            annot_list = row[annot_key]
            if not isinstance(annot_list, list):
                annot_list = []

            func_dict[annot_key](g, gene_node_label, annot_list)


def process_disease_compound(g, disease_compound):
    """Process disease-compound relationships and add them to the graph.

    :param g: the input graph to extend with gene nodes.
    :param disease_compound: the input DataFrame containing disease_compound relationships.
    """
    for _i, row in disease_compound.iterrows():
        disease_node_label = row["identifier"].replace("_", ":")
        disease_annot_list = json.loads(json.dumps(row[OPENTARGETS_DISEASE_COMPOUND_COL]))

        if isinstance(disease_annot_list, float):
            disease_annot_list = []

        add_opentargets_disease_compound_subgraph(g, disease_node_label, disease_annot_list)


def process_ppi(g, gene_node_label, row):
    """Process protein-protein interactions and add them to the graph.
    
    :param g: the input graph to extend with gene nodes.
    :param gene_node_label: the gene node to be linked to annotation entities.
    :param row: row in the combined DataFrame.
    """
    if STRING_PPI_COL in row and row[STRING_PPI_COL] is not None:
        try:
            ppi_list = json.loads(json.dumps(row[STRING_PPI_COL]))
        except (ValueError, TypeError):
            ppi_list = []

        if isinstance(ppi_list, list) and len(ppi_list) > 0:
            valid_ppi_list = [item for item in ppi_list if pd.notna(item.get("stringdb_link_to"))]
            if valid_ppi_list:
                add_stringdb_ppi_subgraph(g, gene_node_label, valid_ppi_list)

        if not isinstance(ppi_list, float):
            add_stringdb_ppi_subgraph(g, gene_node_label, ppi_list)

    if INTACT_INTERACT_COL in row and row[INTACT_INTERACT_COL] is not None:
        try:
            intact_ppi_list = json.loads(json.dumps(row[INTACT_INTERACT_COL]))
        except (ValueError, TypeError):
            intact_ppi_list = []

        if isinstance(intact_ppi_list, list) and len(intact_ppi_list) > 0:
            valid_intact_ppi_list = [item for item in intact_ppi_list if pd.notna(item.get("id_A")) and pd.notna(item.get("id_B"))]
            if valid_intact_ppi_list:
                add_intact_ppi_subgraph(g, gene_node_label, valid_intact_ppi_list)


def process_homologs(g, combined_df, homolog_df_list, func_dict, dea_columns):
    """Process homolog dataframes and combined df and add them to the graph.

    :param g: the input graph to extend with gene nodes.
    :param combined_df: dataframe without homolog information.
    :param homolog_df_list: list of dataframes from homolog queries.
    :param func_dict: list of functions for node generation.
    :param dea_columns: columns ending with _dea
    """
    func_dict_hl = {}

    for homolog_df in homolog_df_list:
        last_col = homolog_df.columns[-1]
        for key, func in func_dict.items():
            if last_col == key and last_col in combined_df.columns:
                func_dict_hl[last_col] = func

    for _i, row in tqdm(combined_df.iterrows(), total=combined_df.shape[0], desc="Building graph"):
        if pd.isna(row["identifier"]) or pd.isna(row["target"]):
            continue
        gene_node_label = add_gene_node(g, row, dea_columns)
        func_dict_non_hl = {key: func for key, func in func_dict.items() if key not in func_dict_hl}
        process_annotations(g, gene_node_label, row, func_dict_non_hl)
        process_ppi(g, gene_node_label, row)

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
                    annot_node_attrs = ENSEMBL_HOMOLOG_NODE_ATTRS.copy()
                    annot_node_attrs["id"] = homolog_node_label
                    annot_node_attrs["labels"] = HOMOLOG_NODE_LABELS
                    g.add_node(homolog_node_label, attr_dict=annot_node_attrs)

                    process_annotations(g, homolog_node_label, row, func_dict_hl)


def process_homologs(g, combined_df, homolog_df_list, func_dict, dea_columns):
    """Process homolog dataframes and combined df and add them to the graph.

    :param g: the input graph to extend with gene nodes.
    :param combined_df: dataframe without homolog information.
    :param homolog_df_list: list of dataframes from homolog queries.
    :param func_dict: list of functions for node generation.
    :param dea_columns: columns ending with _dea
    """
    func_dict_hl = {}

    for homolog_df in homolog_df_list:
        last_col = homolog_df.columns[-1]
        for key, func in func_dict.items():
            if last_col == key and last_col in combined_df.columns:
                func_dict_hl[last_col] = func

    for _i, row in tqdm(combined_df.iterrows(), total=combined_df.shape[0], desc="Building graph"):
        if pd.isna(row["identifier"]) or pd.isna(row["target"]):
            continue
        gene_node_label = add_gene_node(g, row, dea_columns)
        func_dict_non_hl = {key: func for key, func in func_dict.items() if key not in func_dict_hl}
        process_annotations(g, gene_node_label, row, func_dict_non_hl)
        process_ppi(g, gene_node_label, row)

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
                    annot_node_attrs = ENSEMBL_HOMOLOG_NODE_ATTRS.copy()
                    annot_node_attrs["id"] = homolog_node_label
                    annot_node_attrs["labels"] = HOMOLOG_NODE_LABELS
                    g.add_node(homolog_node_label, attr_dict=annot_node_attrs)

                    process_annotations(g, homolog_node_label, row, func_dict_hl)


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
    """
    g = nx.MultiDiGraph()
    combined_df = combined_df[(combined_df["target.source"] == "Ensembl")]

    dea_columns = [c for c in combined_df.columns if c.endswith("_dea")]

    func_dict = {
        BGEE_GENE_EXPRESSION_LEVELS_COL: add_gene_bgee_subgraph,
        DISGENET_DISEASE_COL: add_disgenet_gene_disease_subgraph,
        LITERATURE_DISEASE_COL: add_literature_gene_disease_subgraph,
        MINERVA: add_minerva_gene_pathway_subgraph,
        WIKIPATHWAYS: add_wikipathways_gene_pathway_subgraph,
        INTACT_COMPOUND_INTERACT_COL: add_intact_compound_subgraph,
        KEGG_COL: add_kegg_gene_pathway_subgraph,
        OPENTARGETS_REACTOME_COL: add_opentargets_gene_reactome_pathway_subgraph,
        OPENTARGETS_GO_COL: add_opentargets_gene_go_subgraph,
        OPENTARGETS_GENE_COMPOUND_COL: add_opentargets_gene_compound_subgraph,
        MOLMEDB_PROTEIN_COMPOUND_COL: add_molmedb_gene_inhibitor_subgraph,
        PUBCHEM_COMPOUND_ASSAYS_COL: add_pubchem_assay_subgraph,
        ENSEMBL_HOMOLOG_COL: add_ensembl_homolog_subgraph,
    }

    if homolog_df_list is not None:
        process_homologs(g, combined_df, homolog_df_list, func_dict, dea_columns)

    if homolog_df_list is None:
        for _i, row in tqdm(
            combined_df.iterrows(), total=combined_df.shape[0], desc="Building graph"
        ):
            if pd.isna(row["identifier"]) or pd.isna(row["target"]):
                continue
            gene_node_label = add_gene_node(g, row, dea_columns)
            process_annotations(g, gene_node_label, row, func_dict)
            process_ppi(g, gene_node_label, row)

    if disease_compound is not None:
        process_disease_compound(g, disease_compound)

    if pathway_compound is not None:
        process_kegg_pathway_compound(g, pathway_compound, combined_df)

    normalize_node_attributes(g)
    normalize_edge_attributes(g)

    return g


def save_graph(
    combined_df: pd.DataFrame,
    combined_metadata: Dict[Any, Any],
    disease_compound: pd.DataFrame = None,
    graph_name: str = "combined",
    graph_dir: str = "examples/usecases/",
):
    """Save the graph to a file.

    :param combined_df: the input DataFrame to be converted into a graph.
    :param combined_metadata: the metadata of the graph.
    :param disease_compound: the input DataFrame containing disease-compound relationships.
    :param graph_name: the name of the graph.
    :param graph_dir: the directory to save the graph.
    :returns: a NetworkX MultiDiGraph

    """
    graph_path = f"{graph_dir}/{graph_name}"
    os.makedirs(graph_path, exist_ok=True)

    df_path = f"{graph_path}/{graph_name}_df.pkl"
    metadata_path = f"{graph_path}/{graph_name}_metadata.pkl"
    graph_path_pickle = f"{graph_path}/{graph_name}_graph.pkl"
    graph_path_gml = f"{graph_path}/{graph_name}_graph.gml"

    # Save the combined DataFrame
    combined_df.to_pickle(df_path)
    logger.warning(f"Combined DataFrame saved in {df_path}")

    # Save the metadata
    with open(metadata_path, "wb") as file:
        pickle.dump(combined_metadata, file)
    logger.warning(f"Metadata saved in {metadata_path}")

    # Save the graph
    g = build_networkx_graph(combined_df, disease_compound)
    logger.warning("Graph is built successfully")

    with open(graph_path_pickle, "wb") as f:
        pickle.dump(g, f)
    nx.write_gml(g, graph_path_gml)
    logger.warning(f"Graph saved in {graph_path_pickle} and {graph_path_gml}")

    return g
