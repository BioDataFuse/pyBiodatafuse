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
        if pd.isna(annot["pathway_label"]):
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
    if STRING_PPI_COL in row:
        ppi_list = json.loads(json.dumps(row[STRING_PPI_COL]))
        for item in ppi_list:
            if pd.isna(item["stringdb_link_to"]):
                ppi_list = []

        if not isinstance(ppi_list, float):
            add_stringdb_ppi_subgraph(g, gene_node_label, ppi_list)


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
) -> nx.MultiDiGraph:
    """Construct a NetWorkX graph from a Pandas DataFrame of genes and their multi-source annotations.

    :param combined_df: the input DataFrame to be converted into a graph.
    :param disease_compound: the input DataFrame containing disease-compound relationships.
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
        OPENTARGETS_REACTOME_COL: add_opentargets_gene_reactome_pathway_subgraph,
        OPENTARGETS_GO_COL: add_opentargets_gene_go_subgraph,
        OPENTARGETS_GENE_COMPOUND_COL: add_opentargets_gene_compound_subgraph,
        MOLMEDB_PROTEIN_COMPOUND_COL: add_molmedb_gene_inhibitor_subgraph,
        PUBCHEM_COMPOUND_ASSAYS_COL: add_pubchem_assay_subgraph,
    }

    for _i, row in tqdm(combined_df.iterrows(), total=combined_df.shape[0], desc="Building graph"):
        if pd.isna(row["identifier"]) or pd.isna(row["target"]):
            continue
        gene_node_label = add_gene_node(g, row, dea_columns)
        process_annotations(g, gene_node_label, row, func_dict)
        process_ppi(g, gene_node_label, row)

    if disease_compound is not None:
        process_disease_compound(g, disease_compound)

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
