# coding: utf-8

"""Python module to construct a NetworkX graph from the annotated data frame."""

import json
import pickle

import networkx as nx
import pandas as pd

from pyBiodatafuse.constants import (
    BGEE,
    BGEE_EDGE_ATTRS,
    BGEE_EDGE_LABEL,
    BGEE_NODE_ATTRS,
    BGEE_NODE_MAIN_LABEL,
    DISGENET,
    DISGENET_EDGE_ATTRS,
    DISGENET_EDGE_LABEL,
    DISGENET_NODE_ATTRS,
    DISGENET_NODE_MAIN_LABEL,
    GENE_NODE_LABELS,
    MINERVA,
    MINERVA_EDGE_ATTRS,
    MINERVA_EDGE_LABEL,
    MINERVA_NODE_ATTRS,
    MINERVA_NODE_MAIN_LABEL,
    MOLMEDB_INHIBITOR_COL,
    MOLMEDB_PROTEIN_COMPOUND_EDGE_ATTRS,
    MOLMEDB_PROTEIN_COMPOUND_EDGE_LABEL,
    MOLMEDB_PROTEIN_COMPOUND_NODE_ATTRS,
    MOLMEDB_PROTEIN_COMPOUND_NODE_MAIN_LABEL,
    OPENTARGETS_COMPOUND_COL,
    OPENTARGETS_COMPOUND_EDGE_ATTRS,
    OPENTARGETS_COMPOUND_NODE_ATTRS,
    OPENTARGETS_COMPOUND_NODE_MAIN_LABEL,
    OPENTARGETS_DISEASE_COL,
    OPENTARGETS_DISEASE_EDGE_ATTRS,
    OPENTARGETS_DISEASE_EDGE_LABEL,
    OPENTARGETS_DISEASE_NODE_ATTRS,
    OPENTARGETS_DISEASE_NODE_MAIN_LABEL,
    OPENTARGETS_GO_BP_NODE_LABELS,
    OPENTARGETS_GO_CC_NODE_LABELS,
    OPENTARGETS_GO_COL,
    OPENTARGETS_GO_EDGE_ATTRS,
    OPENTARGETS_GO_EDGE_LABEL,
    OPENTARGETS_GO_MF_NODE_LABELS,
    OPENTARGETS_GO_NODE_ATTRS,
    OPENTARGETS_GO_NODE_MAIN_LABEL,
    OPENTARGETS_REACTOME_COL,
    OPENTARGETS_REACTOME_EDGE_ATTRS,
    OPENTARGETS_REACTOME_EDGE_LABEL,
    OPENTARGETS_REACTOME_NODE_ATTRS,
    OPENTARGETS_REACTOME_NODE_MAIN_LABEL,
    OPENTARGETS_SIDE_EFFECT_EDGE_LABEL,
    OPENTARGETS_SIDE_EFFECT_NODE_ATTRS,
    OPENTARGETS_SIDE_EFFECT_NODE_MAIN_LABEL,
    PUBCHEM_ASSAYS_COL,
    PUBCHEM_EDGE_ATTRS,
    PUBCHEM_NODE_ATTRS,
    PUBCHEM_NODE_MAIN_LABEL,
    STRING,
    STRING_EDGE_ATTRS,
    STRING_EDGE_LABEL,
    STRING_EDGE_MAIN_LABEL,
    WIKIPATHWAYS,
    WIKIPATHWAYS_EDGE_ATTRS,
    WIKIPATHWAYS_EDGE_LABEL,
    WIKIPATHWAYS_NODE_ATTRS,
    WIKIPATHWAYS_NODE_MAIN_LABEL,
)


def load_dataframe_from_pickle(pickle_path: str) -> pd.DataFrame:
    """Load a previously annotated dataframe from a pickle file.

    :param pickle_path: the path to a previously obtained annotation dataframe dumped as a pickle file.
    :returns: a Pandas dataframe.
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


def add_bgee_subgraph(g, gene_node_label, annot_list):
    """Construct part of the graph by linking the gene to a list of annotation entities (disease, drug ..etc).

    :param g: the input graph to extend with new nodes and edges.
    :param gene_node_label: the gene node to be linked to annotation entities.
    :param annot_list: list of annotations from a specific source (e.g. DisGeNET, WikiPathways ..etc).
    :returns: a NetworkX MultiDiGraph
    """
    for annot in annot_list:
        if not pd.isna(annot["anatomical_entity_name"]):
            annot_node_label = annot[BGEE_NODE_MAIN_LABEL]
            annot_node_attrs = BGEE_NODE_ATTRS.copy()
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
            node_exists = [
                x for x, y in edge_data.items() if y["attr_dict"]["edge_hash"] == edge_hash
            ]

            if len(node_exists) == 0:
                g.add_edge(
                    gene_node_label, annot_node_label, label=BGEE_EDGE_LABEL, attr_dict=edge_attrs
                )

    return g


def add_disgenet_disease_subgraph(g, gene_node_label, annot_list):
    """Construct part of the graph by linking the gene to a list of annotation entities (disease, drug ..etc).

    :param g: the input graph to extend with new nodes and edges.
    :param gene_node_label: the gene node to be linked to annotation entities.
    :param annot_list: list of annotations from a specific source (e.g. DisGeNET, WikiPathways ..etc).
    :returns: a NetworkX MultiDiGraph
    """
    for annot in annot_list:
        if not pd.isna(annot["disease_name"]):
            annot_node_label = annot[DISGENET_NODE_MAIN_LABEL]
            annot_node_attrs = DISGENET_NODE_ATTRS.copy()
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
            edge_attrs["score"] = edge_attrs["score"]

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
                    label=DISGENET_EDGE_LABEL,
                    attr_dict=edge_attrs,
                )

    return g


# TODO: The annotations are not curated and will be used again when the OpenTarget annotation improves.
# def add_opentargets_disease_subgraph(g, gene_node_label, annot_list):  # TODO: should be updated
#     """Construct part of the graph by linking the gene to a list of annotation entities (disease, drug ..etc).

#     :param g: the input graph to extend with new nodes and edges.
#     :param gene_node_label: the gene node to be linked to annotation entities.
#     :param annot_list: list of annotations from a specific source (e.g. DisGeNET, WikiPathways ..etc).
#     :returns: a NetworkX MultiDiGraph
#     """
#     for annot in annot_list:
#         if not pd.isna(annot["disease_name"]):
#             annot_node_label = annot[OPENTARGETS_DISEASE_NODE_MAIN_LABEL]
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
#                     label=OPENTARGETS_DISEASE_EDGE_LABEL,
#                     attr_dict=edge_attrs,
#                 )

#     return g


def add_minerva_subgraph(g, gene_node_label, annot_list):
    """Construct part of the graph by linking the gene to a list of annotation entities (disease, drug ..etc).

    :param g: the input graph to extend with new nodes and edges.
    :param gene_node_label: the gene node to be linked to annotation entities.
    :param annot_list: list of annotations from a specific source (e.g. DisGeNET, WikiPathways ..etc).
    :returns: a NetworkX MultiDiGraph
    """
    for annot in annot_list:
        if not pd.isna(annot["pathway_label"]):
            annot_node_label = annot[MINERVA_NODE_MAIN_LABEL]
            annot_node_attrs = MINERVA_NODE_ATTRS.copy()
            annot_node_attrs["name"] = annot["pathway_label"]
            annot_node_attrs["id"] = annot["pathway_id"]
            annot_node_attrs["gene_count"] = annot["pathway_gene_count"]

            g.add_node(annot_node_label, attr_dict=annot_node_attrs)

            edge_attrs = MINERVA_EDGE_ATTRS.copy()

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
                    label=MINERVA_EDGE_LABEL,
                    attr_dict=edge_attrs,
                )

    return g


def add_wikipathways_subgraph(g, gene_node_label, annot_list):
    """Construct part of the graph by linking the gene to a list of annotation entities (disease, drug ..etc).

    :param g: the input graph to extend with new nodes and edges.
    :param gene_node_label: the gene node to be linked to annotation entities.
    :param annot_list: list of annotations from a specific source (e.g. DisGeNET, WikiPathways ..etc).
    :returns: a NetworkX MultiDiGraph
    """
    for annot in annot_list:
        if not pd.isna(annot["pathway_label"]):
            annot_node_label = annot[WIKIPATHWAYS_NODE_MAIN_LABEL]
            annot_node_attrs = WIKIPATHWAYS_NODE_ATTRS.copy()
            annot_node_attrs["name"] = annot["pathway_label"]
            annot_node_attrs["id"] = annot["pathway_id"]
            annot_node_attrs["gene_count"] = annot["pathway_gene_count"]

            g.add_node(annot_node_label, attr_dict=annot_node_attrs)

            edge_attrs = WIKIPATHWAYS_EDGE_ATTRS.copy()

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
                    label=WIKIPATHWAYS_EDGE_LABEL,
                    attr_dict=edge_attrs,
                )

    return g


def add_opentargets_reactome_pathway_subgraph(g, gene_node_label, annot_list):
    """Construct part of the graph by linking the gene to a list of annotation entities (disease, drug ..etc).

    :param g: the input graph to extend with new nodes and edges.
    :param gene_node_label: the gene node to be linked to annotation entities.
    :param annot_list: list of annotations from a specific source (e.g. DisGeNET, WikiPathways ..etc).
    :returns: a NetworkX MultiDiGraph
    """
    for annot in annot_list:
        if not pd.isna(annot["pathway_id"]):
            annot_node_label = annot[OPENTARGETS_REACTOME_NODE_MAIN_LABEL]
            annot_node_attrs = OPENTARGETS_REACTOME_NODE_ATTRS.copy()
            annot_node_attrs["name"] = annot["pathway_label"]
            annot_node_attrs["id"] = annot["pathway_id"]

            g.add_node(annot_node_label, attr_dict=annot_node_attrs)

            edge_attrs = OPENTARGETS_REACTOME_EDGE_ATTRS.copy()

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
                    label=OPENTARGETS_REACTOME_EDGE_LABEL,
                    attr_dict=edge_attrs,
                )

    return g


def add_opentargets_go_subgraph(g, gene_node_label, annot_list):
    """Construct part of the graph by linking the gene to a list of annotation entities (disease, drug ..etc).

    :param g: the input graph to extend with new nodes and edges.
    :param gene_node_label: the gene node to be linked to annotation entities.
    :param annot_list: list of annotations from a specific source (e.g. DisGeNET, WikiPathways ..etc).
    :returns: a NetworkX MultiDiGraph
    """
    for annot in annot_list:
        annot_node_label = annot[OPENTARGETS_GO_NODE_MAIN_LABEL]
        annot_node_attrs = OPENTARGETS_GO_NODE_ATTRS.copy()
        annot_node_attrs["name"] = annot["go_name"]
        annot_node_attrs["id"] = annot["go_id"]
        if annot["go_type"] == "BP":
            annot_node_attrs["labels"] = OPENTARGETS_GO_BP_NODE_LABELS
        elif annot["go_type"] == "MF":
            annot_node_attrs["labels"] = OPENTARGETS_GO_MF_NODE_LABELS
        elif annot["go_type"] == "CC":
            annot_node_attrs["labels"] = OPENTARGETS_GO_CC_NODE_LABELS

        g.add_node(annot_node_label, attr_dict=annot_node_attrs)

        edge_attrs = OPENTARGETS_GO_EDGE_ATTRS.copy()

        edge_hash = hash(frozenset(edge_attrs.items()))
        edge_attrs["edge_hash"] = edge_hash
        edge_data = g.get_edge_data(gene_node_label, annot_node_label)
        edge_data = {} if edge_data is None else edge_data
        node_exists = [x for x, y in edge_data.items() if y["attr_dict"]["edge_hash"] == edge_hash]

        if len(node_exists) == 0:
            g.add_edge(
                gene_node_label,
                annot_node_label,
                label=OPENTARGETS_GO_EDGE_LABEL,
                attr_dict=edge_attrs,
            )

    return g


def add_opentargets_compound_subgraph(g, gene_node_label, annot_list):
    """Construct part of the graph by linking the gene to a list of annotation entities (disease, drug ..etc).

    :param g: the input graph to extend with new nodes and edges.
    :param gene_node_label: the gene node to be linked to annotation entities.
    :param annot_list: list of annotations from a specific source (e.g. DisGeNET, WikiPathways ..etc).
    :returns: a NetworkX MultiDiGraph
    """
    for annot in annot_list:
        if not pd.isna(annot["relation"]):
            if not pd.isna(annot[OPENTARGETS_COMPOUND_NODE_MAIN_LABEL]):
                annot_node_label = annot[OPENTARGETS_COMPOUND_NODE_MAIN_LABEL]
            else:
                annot_node_label = annot["chembl_id"]
            annot_node_attrs = OPENTARGETS_COMPOUND_NODE_ATTRS.copy()
            annot_node_attrs["name"] = annot["compound_name"]
            if not pd.isna(annot[OPENTARGETS_COMPOUND_NODE_MAIN_LABEL]):
                annot_node_attrs["id"] = annot[OPENTARGETS_COMPOUND_NODE_MAIN_LABEL]
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

            edge_attrs = OPENTARGETS_COMPOUND_EDGE_ATTRS.copy()
            edge_attrs["label"] = annot["relation"]
            edge_hash = hash(frozenset(edge_attrs.items()))
            edge_attrs["edge_hash"] = edge_hash
            edge_data = g.get_edge_data(annot_node_label, gene_node_label)
            edge_data = {} if edge_data is None else edge_data
            node_exists = [
                x for x, y in edge_data.items() if y["attr_dict"]["edge_hash"] == edge_hash
            ]

            if len(node_exists) == 0:
                g.add_edge(
                    annot_node_label,
                    gene_node_label,
                    label=annot["relation"],
                    attr_dict=edge_attrs,
                )

            # Add side effects
            if annot["adverse_effect"]:
                for effect in annot[OPENTARGETS_SIDE_EFFECT_NODE_MAIN_LABEL]:
                    effect_node_label = effect["name"]
                    effect_node_attrs = OPENTARGETS_SIDE_EFFECT_NODE_ATTRS.copy()
                    effect_node_attrs["name"] = effect["name"]

                    g.add_node(effect_node_label, attr_dict=effect_node_attrs)

                    new_edge = (annot_node_label, effect_node_label)

                    # Check if the edge already exists
                    if not g.has_edge(*new_edge):
                        # Add the edge to the graph
                        g.add_edge(
                            annot_node_label,
                            effect_node_label,
                            label=OPENTARGETS_SIDE_EFFECT_EDGE_LABEL,
                        )
    return g


def add_molmedb_gene_inhibitor(g, gene_node_label, annot_list):
    """Construct part of the graph by linking the gene to a list of annotation entities (disease, compound ..etc).

    :param g: the input graph to extend with new nodes and edges.
    :param gene_node_label: the gene node to be linked to annotation entities.
    :param annot_list: list of annotations from a specific source (e.g. DisGeNET, WikiPathways ..etc).
    :returns: a NetworkX MultiDiGraph
    """
    for annot in annot_list:
        if not pd.isna(annot["compound_name"]):
            if not pd.isna(annot[MOLMEDB_PROTEIN_COMPOUND_NODE_MAIN_LABEL]):
                annot_node_label = annot[MOLMEDB_PROTEIN_COMPOUND_NODE_MAIN_LABEL]
            else:
                annot_node_label = annot["molmedb_id"]

            annot_node_attrs = MOLMEDB_PROTEIN_COMPOUND_NODE_ATTRS.copy()
            annot_node_attrs["name"] = annot["compound_name"]

            if not pd.isna(annot[MOLMEDB_PROTEIN_COMPOUND_NODE_MAIN_LABEL]):
                annot_node_attrs["id"] = annot[MOLMEDB_PROTEIN_COMPOUND_NODE_MAIN_LABEL]
            else:
                annot_node_attrs["id"] = annot["molmedb_id"]

            annot_node_attrs["molmedb_id"] = annot["molmedb_id"]
            if not pd.isna(annot["chebi_id"]):
                annot_node_attrs["chebi_id"] = annot["chebi_id"]
            if not pd.isna(annot["drugbank_id"]):
                annot_node_attrs["drugbank_id"] = annot["drugbank_id"]
            if not pd.isna(annot["compound_cid"]):
                annot_node_attrs["compound_cid"] = annot["compound_cid"]
            # if not pd.isna(annot["pdb_ligand_id"]):
            #     annot_node_attrs["pdb_ligand_id"] = annot["pdb_ligand_id"]
            if not pd.isna(annot["inchikey"]):
                annot_node_attrs["inchikey"] = annot["inchikey"]
            if not pd.isna(annot["smiles"]):
                annot_node_attrs["smiles"] = annot["smiles"]
            if not pd.isna(annot["source_doi"]):
                annot_node_attrs["source_doi"] = annot["source_doi"]
            if not pd.isna(annot["source_pmid"]):
                annot_node_attrs["source_pmid"] = annot["source_pmid"]

            merge_node(g, annot_node_label, annot_node_attrs)

            edge_attrs = MOLMEDB_PROTEIN_COMPOUND_EDGE_ATTRS.copy()

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
                    label=MOLMEDB_PROTEIN_COMPOUND_EDGE_LABEL,
                    attr_dict=edge_attrs,
                )

    return g


def add_pubchem_assay(g, gene_node_label, annot_list):
    """Construct part of the graph by linking the gene to a list of annotation entities (disease, drug ..etc).

    :param g: the input graph to extend with new nodes and edges.
    :param gene_node_label: the gene node to be linked to annotation entities.
    :param annot_list: list of annotations from a specific source (e.g. DisGeNET, WikiPathways ..etc).
    :returns: a NetworkX MultiDiGraph
    """
    for annot in annot_list:
        if not pd.isna(annot["pubchem_assay_id"]):
            annot_node_label = annot[PUBCHEM_NODE_MAIN_LABEL]
            annot_node_attrs = PUBCHEM_NODE_ATTRS.copy()
            annot_node_attrs["name"] = annot["compound_name"]
            annot_node_attrs["id"] = annot["compound_cid"]
            annot_node_attrs["inchi"] = annot["inchi"]
            if not pd.isna(annot["smiles"]):
                annot_node_attrs["smiles"] = annot["smiles"]

            # g.add_node(annot_node_label, attr_dict=annot_node_attrs)
            merge_node(g, annot_node_label, annot_node_attrs)

            edge_attrs = PUBCHEM_EDGE_ATTRS.copy()
            edge_attrs["assay_type"] = annot["assay_type"]
            edge_attrs["pubchem_assay_id"] = annot["pubchem_assay_id"]
            edge_attrs["outcome"] = annot["outcome"]
            edge_attrs["label"] = annot["outcome"]

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
                    label=annot["outcome"],
                    attr_dict=edge_attrs,
                )

    return g


def add_ppi_subgraph(g, gene_node_label, annot_list):
    """Construct part of the graph by linking the gene to a list of annotation entities (disease, drug ..etc).

    :param g: the input graph to extend with new nodes and edges.
    :param gene_node_label: the gene node to be linked to annotation entities.
    :param annot_list: list of annotations from a specific source (e.g. DisGeNET, WikiPathways ..etc).
    :returns: a NetworkX MultiDiGraph
    """
    for ppi in annot_list:
        edge_attrs = STRING_EDGE_ATTRS.copy()
        edge_attrs["score"] = ppi["score"]

        edge_hash = hash(frozenset(edge_attrs.items()))
        edge_attrs["edge_hash"] = edge_hash
        edge_data = g.get_edge_data(gene_node_label, ppi[STRING_EDGE_MAIN_LABEL])

        edge_data = {} if edge_data is None else edge_data
        node_exists = [x for x, y in edge_data.items() if y["attr_dict"]["edge_hash"] == edge_hash]
        if len(node_exists) == 0:
            g.add_edge(
                gene_node_label,
                ppi[STRING_EDGE_MAIN_LABEL],
                label=STRING_EDGE_LABEL,
                attr_dict=edge_attrs,
            )

    return g


# def add_drug_disease_subgraph(g, drug_node_label, annot_list):
#     """Construct part of the graph by linking the gene to a list of annotation entities (disease, drug ..etc).

#     :param g: the input graph to extend with new nodes and edges.
#     :param drug_node_label: the gene node to be linked to annotation entities.
#     :param annot_list: list of annotations from a specific source (e.g. DisGeNET, WikiPathways ..etc).
#     :returns: a NetworkX MultiDiGraph
#     """
#     for ddi in annot_list:
#         edge_attrs = {"source": OPENTARGETS, "label": "treated_with"}

#         edge_hash = hash(frozenset(edge_attrs.items()))
#         edge_attrs["edge_hash"] = edge_hash

#         new_edge = (ddi["umls"], drug_node_label)

#         # Check if the edge already exists
#         if not g.has_edge(*new_edge):
#             # Add the edge to the graph
#             g.add_edge(ddi["umls"], drug_node_label, label="treated_with", attr_dict=edge_attrs)

#     return g


def networkx_graph(fuse_df: pd.DataFrame, drug_disease=None):
    """Construct a NetWorkX graph from a Pandas DataFrame of genes and their multi-source annotations.

    :param fuse_df: the input dataframe to be converted into a graph.
    :param drug_disease: the input dataframe containing drug_disease relationships
    :returns: a NetworkX MultiDiGraph
    """
    g = nx.MultiDiGraph()

    dea_columns = [c for c in fuse_df.columns if c.endswith("_dea")]

    func_dict = {
        BGEE: add_bgee_subgraph,
        DISGENET: add_disgenet_disease_subgraph,
        OPENTARGETS_DISEASE_COL: add_opentargets_disease_subgraph,
        MINERVA: add_minerva_subgraph,
        WIKIPATHWAYS: add_wikipathways_subgraph,
        OPENTARGETS_REACTOME_COL: add_opentargets_reactome_pathway_subgraph,
        OPENTARGETS_GO_COL: add_opentargets_go_subgraph,
        OPENTARGETS_COMPOUND_COL: add_opentargets_compound_subgraph,
        MOLMEDB_INHIBITOR_COL: add_molmedb_gene_inhibitor,
        PUBCHEM_ASSAYS_COL: add_pubchem_assay,
    }

    for _i, row in fuse_df.iterrows():
        if pd.notna(row["identifier"]) and pd.notna(row["target"]):
            gene_node_label = row["identifier"]
            gene_node_attrs = {
                "source": "BridgeDB",
                "name": row["identifier"],
                "id": row["target"],
                "labels": GENE_NODE_LABELS,
                row["target.source"]: row["target"],
            }

            for c in dea_columns:
                gene_node_attrs[c[:-4]] = row[c]

            g.add_node(gene_node_label, attr_dict=gene_node_attrs)

            for annot_key in func_dict:
                if annot_key in row:
                    annot_list = json.loads(json.dumps(row[annot_key]))
                    if not isinstance(annot_list, list):
                        annot_list = []

                    func_dict[annot_key](g, gene_node_label, annot_list)

            if STRING in row:
                gene_node_label = row["identifier"]
                ppi_list = json.loads(json.dumps(row[STRING]))

                if ppi_list is None:
                    ppi_list = []

                if not isinstance(ppi_list, float):
                    add_ppi_subgraph(g, gene_node_label, ppi_list)
    # TODO:
    # if drug_disease is not None:
    #     fuse_df = pd.concat(
    #         [fuse_df, drug_disease[["identifier", "drug_diseases"]]], ignore_index=True
    #     )
    #     if "drug_diseases" in row:
    #         for _i, row in fuse_df.iterrows():
    #             gene_node_label_2 = row["identifier"]
    #             ddi_list = json.loads(json.dumps(row["drug_diseases"]))

    #             if type(ddi_list) == float:
    #                 ddi_list = []

    #             add_drug_disease_subgraph(g, gene_node_label_2, ddi_list)

    for node in g.nodes():
        if "attr_dict" in g.nodes[node]:
            for k, v in g.nodes[node]["attr_dict"].items():
                if v is not None:
                    g.nodes[node][k] = v

            del g.nodes[node]["attr_dict"]

    for u, v, k in g.edges(keys=True):
        if "attr_dict" in g[u][v][k]:
            for x, y in g[u][v][k]["attr_dict"].items():
                if y is not None and x != "edge_hash":
                    g[u][v][k][x] = y

            del g[u][v][k]["attr_dict"]

    return g
