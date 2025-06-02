# coding: utf-8

"""Python module to construct a NetworkX graph from the annotated data frame."""

import json
import pickle

import networkx as nx
import pandas as pd
from tqdm import tqdm

import pyBiodatafuse.constants as Cons


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
    """Merge and/or adding nodes to the graph.

    :param g: the graph to which the node will be added.
    :param node_label: node label.
    :param node_attrs: dictionary of node attributes.
    """
    # Ensure all node attributes are strings or set to None
    node_attrs = {k: (v if isinstance(v, str) else "") for k, v in node_attrs.items()}
    if node_label not in g.nodes():
        g.add_node(node_label, attr_dict=node_attrs)
    else:
        if "attr_dict" in g.nodes[node_label]:
            for k, v in node_attrs.items():
                node_attr_dict = g.nodes[node_label]["attr_dict"]
                if k in node_attr_dict:
                    if node_attr_dict[k] is not None:
                        if isinstance(v, str):
                            v_list = node_attr_dict[k].split("|")
                            v_list.append(v)
                            node_attr_dict[k] = "|".join(list(set(v_list)))
                    else:
                        node_attr_dict[k] = v
                else:
                    node_attr_dict[k] = v
        else:
            g.add_node(node_label, attr_dict=node_attrs)


def add_gene_bgee_subgraph(g: nx.MultiDiGraph, gene_node_label: str, annot_list: list):
    """Construct part of the graph by linking the gene to a list of anatomical entities.

    :param g: the input graph to extend with new nodes and edges.
    :param gene_node_label: the gene node to be linked to annotation entities.
    :param annot_list: list of anatomical entities from Bgee with gene expression levels.
    :returns: a NetworkX MultiDiGraph
    """
    for annot in annot_list:
        if pd.isna(annot["anatomical_entity_name"]):
            continue

        annot_node_label = annot[Cons.BGEE_ANATOMICAL_NODE_MAIN_LABEL]
        entity_attrs = Cons.BGEE_ANATOMICAL_NODE_ATTRS.copy()
        entity_attrs.update(
            {
                "name": annot["anatomical_entity_name"],
                "id": annot["anatomical_entity_id"],
                "datasource": Cons.BGEE,
            }
        )

        g.add_node(annot_node_label, attr_dict=entity_attrs)

        edge_attrs = Cons.BGEE_EDGE_ATTRS.copy()
        fields = {
            "confidence_level_id": ("confidence_level_name", "confidence_level_id"),
            "expression_level": ("expression_level",),
            "developmental_stage_name": ("developmental_stage_name",),
            "developmental_stage_id": ("developmental_stage_id",),
        }

        for field, attrs in fields.items():
            if pd.notna(annot[field]):
                for attr in attrs:
                    edge_attrs[attr] = annot[attr]

        edge_hash = hash(frozenset(edge_attrs.items()))
        edge_attrs["edge_hash"] = edge_hash  # type: ignore
        edge_data = g.get_edge_data(gene_node_label, annot_node_label)
        edge_data = {} if edge_data is None else edge_data
        edge_exists = any(
            y.get("attr_dict", {}).get("edge_hash") == edge_hash for y in edge_data.values()
        )

        if not edge_exists:
            g.add_edge(
                gene_node_label,
                annot_node_label,
                label=Cons.BGEE_GENE_ANATOMICAL_EDGE_LABEL,
                attr_dict=edge_attrs,
            )

    return g


def add_disgenet_gene_disease_subgraph(g: nx.MultiDiGraph, gene_node_label: str, annot_list: list):
    """Construct part of the graph by linking the gene to diseases.

    :param g: the input graph to extend with new nodes and edges.
    :param gene_node_label: the gene node to be linked to diseases.
    :param annot_list: list of diseases from DisGeNET.
    :returns: a NetworkX MultiDiGraph
    """
    for annot in annot_list:
        if pd.isna(annot["disease_name"]):
            continue

        annot_node_label = annot[Cons.DISEASE_NODE_MAIN_LABEL]
        annot_node_attrs = Cons.DISGENET_DISEASE_NODE_ATTRS.copy()
        annot_node_attrs.update(
            {
                "name": annot["disease_name"],
                "id": annot["UMLS"],
                "datasource": Cons.DISGENET,
            }
        )

        other_ids = {
            "HPO": annot["HPO"],
            "NCI": annot["NCI"],
            "OMIM": annot["OMIM"],
            "MONDO": annot["MONDO"],
            "ORDO": annot["ORDO"],
            "EFO": annot["EFO"],
            "DO": annot["DO"],
            "MESH": annot["MESH"],
            "UMLS": annot["UMLS"],
            "disease_type": annot["disease_type"],
        }

        for key, value in other_ids.items():
            if pd.notna(value):
                annot_node_attrs[key] = value

        g.add_node(annot_node_label, attr_dict=annot_node_attrs)

        edge_attrs = Cons.DISGENET_EDGE_ATTRS.copy()
        edge_attrs["score"] = annot["score"]

        if pd.notna(annot["ei"]):
            edge_attrs["ei"] = annot["ei"]
        if pd.notna(annot["el"]):
            edge_attrs["el"] = annot["el"]

        edge_hash = hash(frozenset(edge_attrs.items()))
        edge_attrs["edge_hash"] = edge_hash  # type: ignore
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


def add_literature_gene_disease_subgraph(
    g: nx.MultiDiGraph, gene_node_label: str, annot_list: list
):
    """Construct part of the graph by linking the gene to diseases form literature.

    :param g: the input graph to extend with new nodes and edges.
    :param gene_node_label: the gene node to be linked to diseases.
    :param annot_list: list of diseases from DisGeNET.
    :returns: a NetworkX MultiDiGraph
    """
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
        edge_attrs["edge_hash"] = edge_hash  # type: ignore
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
#             annot_node_attrs.update(
#                 {
#                     "name": annot["disease_name"],
#                     "id": annot["disease_id"],
#                     "therapeutic_areas": annot["therapeutic_areas"],
#                 }
#             )

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


def add_minerva_gene_pathway_subgraph(g: nx.MultiDiGraph, gene_node_label: str, annot_list: list):
    """Construct part of the graph by linking the gene to MINERVA pathways.

    :param g: the input graph to extend with new nodes and edges.
    :param gene_node_label: the gene node to be linked to MINERVA pathways.
    :param annot_list: list of MINERVA pathways from MINERVA.
    :returns: a NetworkX MultiDiGraph
    """
    for annot in annot_list:
        if pd.isna(annot.get("pathway_label")):
            continue

        annot_node_label = annot[Cons.PATHWAY_NODE_MAIN_LABEL]
        annot_node_attrs = Cons.PATHWAY_NODE_ATTRS.copy()
        annot_node_attrs.update(
            {
                "datasource": Cons.MINERVA,
                "name": annot["pathway_label"],
                "id": annot["pathway_id"],
                "gene_count": annot["pathway_gene_count"],
            }
        )

        g.add_node(annot_node_label, attr_dict=annot_node_attrs)

        edge_attrs = Cons.GENE_PATHWAY_EDGE_ATTRS.copy()
        edge_attrs["datasource"] = Cons.MINERVA

        edge_hash = hash(frozenset(edge_attrs.items()))
        edge_attrs["edge_hash"] = edge_hash  # type: ignore
        edge_data = g.get_edge_data(gene_node_label, annot_node_label)
        edge_data = {} if edge_data is None else edge_data
        node_exists = [x for x, y in edge_data.items() if y["attr_dict"]["edge_hash"] == edge_hash]

        if len(node_exists) == 0:
            g.add_edge(
                gene_node_label,
                annot_node_label,
                label=Cons.GENE_PATHWAY_EDGE_LABEL,
                attr_dict=edge_attrs,
            )

    return g


def add_wikipathways_gene_pathway_subgraph(
    g: nx.MultiDiGraph, gene_node_label: str, annot_list: list
):
    """Construct part of the graph by linking the gene to pathways from WikiPathways.

    :param g: the input graph to extend with new nodes and edges.
    :param gene_node_label: the gene node to be linked to pathways from WikiPathways.
    :param annot_list: list of pathways from WikiPathways.
    :returns: a NetworkX MultiDiGraph
    """
    for annot in annot_list:
        if pd.isna(annot["pathway_label"]):
            continue

        annot_node_label = annot[Cons.PATHWAY_NODE_MAIN_LABEL]
        annot_node_attrs = Cons.PATHWAY_NODE_ATTRS.copy()
        annot_node_attrs.update(
            {
                "datasource": Cons.WIKIPATHWAYS,
                "name": annot["pathway_label"],
                "id": annot["pathway_id"],
                "gene_count": annot["pathway_gene_count"],
            }
        )

        g.add_node(annot_node_label, attr_dict=annot_node_attrs)

        edge_attrs = Cons.GENE_PATHWAY_EDGE_ATTRS.copy()
        edge_attrs["datasource"] = Cons.WIKIPATHWAYS

        edge_hash = hash(frozenset(edge_attrs.items()))
        edge_attrs["edge_hash"] = edge_hash  # type: ignore
        edge_data = g.get_edge_data(gene_node_label, annot_node_label)
        edge_data = {} if edge_data is None else edge_data
        node_exists = [x for x, y in edge_data.items() if y["attr_dict"]["edge_hash"] == edge_hash]

        if len(node_exists) == 0:
            g.add_edge(
                gene_node_label,
                annot_node_label,
                label=Cons.GENE_PATHWAY_EDGE_LABEL,
                attr_dict=edge_attrs,
            )

    return g


def add_kegg_gene_pathway_subgraph(g: nx.MultiDiGraph, gene_node_label: str, annot_list: list):
    """Construct part of the graph by linking the gene to pathways from KEGG.

    :param g: the input graph to extend with new nodes and edges.
    :param gene_node_label: the gene node to be linked to pathways from KEGG.
    :param annot_list: list of pathways from KEGG.
    :returns: a NetworkX MultiDiGraph
    """
    for annot in annot_list:
        if pd.isna(annot.get("pathway_label")):
            continue

        annot_node_label = annot[Cons.PATHWAY_NODE_MAIN_LABEL]
        annot_node_attrs = Cons.PATHWAY_NODE_ATTRS.copy()
        annot_node_attrs.update(
            {
                "datasource": Cons.KEGG,
                "name": annot["pathway_label"],
                "id": annot["pathway_id"],
                "gene_count": annot["gene_count"],
            }
        )

        # g.add_node(annot_node_label, attr_dict=annot_node_attrs)
        merge_node(g, annot_node_label, annot_node_attrs)

        edge_attrs = Cons.GENE_PATHWAY_EDGE_ATTRS.copy()
        edge_attrs["datasource"] = Cons.KEGG

        edge_hash = hash(frozenset(edge_attrs.items()))
        edge_attrs["edge_hash"] = edge_hash  # type: ignore
        edge_data = g.get_edge_data(gene_node_label, annot_node_label)
        edge_data = {} if edge_data is None else edge_data
        node_exists = [x for x, y in edge_data.items() if y["attr_dict"]["edge_hash"] == edge_hash]

        if len(node_exists) == 0:
            g.add_edge(
                gene_node_label,
                annot_node_label,
                label=Cons.GENE_PATHWAY_EDGE_LABEL,
                attr_dict=edge_attrs,
            )

    return g


def add_kegg_compounds_subgraph(
    g: nx.MultiDiGraph, pathway_node_label: str, compounds_list: list, combined_df: pd.DataFrame
):
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
        annot_node_attrs = Cons.KEGG_COMPOUND_NODE_ATTRS.copy()
        annot_node_attrs.update(
            {
                "id": compound["KEGG_identifier"],
                "label": compound["name"],
            }
        )

        merge_node(g, annot_node_label, annot_node_attrs)

        for _, path_row in combined_df.iterrows():
            pathways = path_row.get("KEGG_pathways", [])

            for pathway in pathways:
                if pathway_node_label != pathway.get("pathway_id"):
                    continue

                if "compounds" not in pathway:
                    continue

                pathway_compounds = [comp["KEGG_identifier"] for comp in pathway["compounds"]]
                if compound["KEGG_identifier"] not in pathway_compounds:
                    continue

                edge_attrs = Cons.KEGG_COMPOUND_EDGE_ATTRS.copy()
                edge_hash = hash(frozenset(edge_attrs.items()))
                edge_attrs["edge_hash"] = edge_hash  # type: ignore
                edge_data = g.get_edge_data(pathway_node_label, annot_node_label)
                edge_data = {} if edge_data is None else edge_data
                node_exists = [
                    x
                    for x, y in edge_data.items()
                    if "attr_dict" in y and y["attr_dict"].get("edge_hash") == edge_hash
                ]

                if len(node_exists) == 0:
                    g.add_edge(
                        pathway_node_label,
                        annot_node_label,
                        label=Cons.KEGG_COMPOUND_EDGE_LABEL,
                        attr_dict=edge_attrs,
                    )

    return g


def process_kegg_pathway_compound(
    g: nx.MultiDiGraph, kegg_pathway_compound: pd.DataFrame, combined_df: pd.DataFrame
):
    """Process pathway-compound relationships from KEGG and add them to the graph.

    :param g: the input graph to extend with pathway-compound relationships.
    :param kegg_pathway_compound: DataFrame containing pathway-compound relationships.
    :param combined_df: DataFrame containing KEGG pathway data.
    """
    for _, row in kegg_pathway_compound.iterrows():
        compound_info = row[Cons.KEGG_COMPOUND_COL]

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

                if not isinstance(pathway_data, list):
                    continue

                for pathway in pathway_data:
                    pathway_id = pathway.get("pathway_id")
                    pathway_compounds = pathway.get("compounds", [])

                    if any(c.get("KEGG_identifier") == compound_id for c in pathway_compounds):
                        add_kegg_compounds_subgraph(g, pathway_id, compounds_list, combined_df)


def add_opentargets_gene_reactome_pathway_subgraph(
    g: nx.MultiDiGraph, gene_node_label: str, annot_list: list
):
    """Construct part of the graph by linking the gene to Reactome pathways.

    :param g: the input graph to extend with new nodes and edges.
    :param gene_node_label: the gene node to be linked to Reactome pathways.
    :param annot_list: list of Reactome pathways from OpenTargets.
    :returns: a NetworkX MultiDiGraph
    """
    for annot in annot_list:
        if pd.isna(annot["pathway_id"]):
            continue

        annot_node_label = annot[Cons.PATHWAY_NODE_MAIN_LABEL]
        annot_node_attrs = Cons.PATHWAY_NODE_ATTRS.copy()
        annot_node_attrs.update(
            {
                "datasource": Cons.OPENTARGETS,
                "name": annot["pathway_label"],
                "id": annot["pathway_id"],
            }
        )

        g.add_node(annot_node_label, attr_dict=annot_node_attrs)

        edge_attrs = Cons.GENE_PATHWAY_EDGE_ATTRS.copy()
        edge_attrs["datasource"] = Cons.OPENTARGETS

        edge_hash = hash(frozenset(edge_attrs.items()))
        edge_attrs["edge_hash"] = edge_hash  # type: ignore
        edge_data = g.get_edge_data(gene_node_label, annot_node_label)
        edge_data = {} if edge_data is None else edge_data
        node_exists = [x for x, y in edge_data.items() if y["attr_dict"]["edge_hash"] == edge_hash]

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
    for annot in annot_list:
        if pd.isna(annot["go_id"]):
            continue

        annot_node_label = annot[Cons.GO_NODE_MAIN_LABEL]
        annot_node_attrs = Cons.GO_NODE_ATTRS.copy()
        annot_node_attrs.update(
            {
                "name": annot["go_name"],
                "id": annot["go_id"],
                "datasource": Cons.OPENTARGETS,
            }
        )

        if annot["go_type"] == "P":
            annot_node_attrs["labels"] = Cons.GO_BP_NODE_LABELS
        elif annot["go_type"] == "F":
            annot_node_attrs["labels"] = Cons.GO_MF_NODE_LABELS
        elif annot["go_type"] == "C":
            annot_node_attrs["labels"] = Cons.GO_CC_NODE_LABELS
        else:
            raise ValueError(f"Invalid GO type: {annot['go_type']}")

        g.add_node(annot_node_label, attr_dict=annot_node_attrs)

        edge_attrs = Cons.GENE_GO_EDGE_ATTRS.copy()

        edge_hash = hash(frozenset(edge_attrs.items()))
        edge_attrs["edge_hash"] = edge_hash  # type: ignore
        edge_data = g.get_edge_data(gene_node_label, annot_node_label)
        edge_data = {} if edge_data is None else edge_data
        node_exists = [x for x, y in edge_data.items() if y["attr_dict"]["edge_hash"] == edge_hash]

        if len(node_exists) == 0:
            g.add_edge(
                gene_node_label,
                annot_node_label,
                label=Cons.GENE_GO_EDGE_LABEL,
                attr_dict=edge_attrs,
            )

    return g


def add_opentargets_compound_side_effect_subgraph(
    g: nx.MultiDiGraph, compound_node_label: str, side_effects_list: list
):
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
            effect_node_attrs = Cons.SIDE_EFFECT_NODE_ATTRS.copy()
            effect_node_attrs.update(
                {
                    "name": effect["name"],
                    "datasource": Cons.OPENTARGETS,
                }
            )

            g.add_node(effect_node_label, attr_dict=effect_node_attrs)

            edge_attrs = Cons.COMPOUND_SIDE_EFFECT_EDGE_ATTRS.copy()
            edge_hash = hash(frozenset(edge_attrs.items()))
            edge_attrs["edge_hash"] = edge_hash  # type: ignore
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
    for annot in annot_list:
        if pd.isna(annot["relation"]):
            continue

        if not pd.isna(annot[Cons.COMPOUND_NODE_MAIN_LABEL]):
            annot_node_label = annot[Cons.COMPOUND_NODE_MAIN_LABEL]
        else:
            annot_node_label = annot["chembl_id"]

        annot_node_attrs = Cons.OPENTARGETS_COMPOUND_NODE_ATTRS.copy()
        annot_node_attrs.update(
            {
                "name": annot["compound_name"],
                "id": annot["chembl_id"],
                "datasource": Cons.OPENTARGETS,
            }
        )

        if not pd.isna(annot[Cons.COMPOUND_NODE_MAIN_LABEL]):
            annot_node_attrs["id"] = annot[Cons.COMPOUND_NODE_MAIN_LABEL]
        else:
            annot_node_attrs["id"] = annot["chembl_id"]
        annot_node_attrs["chembl_id"] = annot["chembl_id"]

        other_info = {
            "drugbank_id": annot["drugbank_id"],
            "compound_cid": annot["compound_cid"],
            "clincal_trial_phase": annot["clincal_trial_phase"],
            "is_approved": annot["is_approved"],
            "adverse_effect_count": annot["adverse_effect_count"],
        }

        for key, value in other_info.items():
            if not pd.isna(value):
                annot_node_attrs[key] = value

        merge_node(g, annot_node_label, annot_node_attrs)

        edge_attrs = Cons.OPENTARGETS_GENE_COMPOUND_EDGE_ATTRS.copy()
        edge_attrs["label"] = annot["relation"]
        edge_hash = hash(frozenset(edge_attrs.items()))
        edge_attrs["edge_hash"] = edge_hash  # type: ignore
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
                g, annot_node_label, annot[Cons.SIDE_EFFECT_NODE_MAIN_LABEL]
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
        edge_data = g.get_edge_data(gene_node_label, annot_node_label)
        edge_data = {} if edge_data is None else edge_data
        node_exists = [x for x, y in edge_data.items() if y["attr_dict"]["edge_hash"] == edge_hash]

        if len(node_exists) == 0:
            g.add_edge(
                annot_node_label,
                gene_node_label,
                label=Cons.MOLMEDB_PROTEIN_COMPOUND_EDGE_LABEL,
                attr_dict=edge_attrs,
            )

    return g


def add_molmedb_compound_gene_subgraph(g, compound_node_label, annot_list):
    """Construct part of the graph by linking the compound to inhibited genes.

    :param g: the input graph to extend with new nodes and edges.
    :param compound_node_label: the compound node to be linked to its inhibitors.
    :param annot_list: list of gene inhibitors from MolMeDB.
    :returns: a NetworkX MultiDiGraph
    """
    for annot in annot_list:
        if pd.isna(annot["compound_name"]):
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
    for annot in annot_list:
        if pd.isna(annot["pubchem_assay_id"]):
            continue

        annot_node_label = annot[Cons.COMPOUND_NODE_MAIN_LABEL]
        annot_node_attrs = Cons.PUBCHEM_COMPOUND_NODE_ATTRS.copy()
        annot_node_attrs.update(
            {
                "name": annot["compound_name"],
                "id": annot["compound_cid"],
                "inchi": annot["inchi"],
                "datasource": Cons.PUBCHEM,
            }
        )
        if not pd.isna(annot["smiles"]):
            annot_node_attrs["smiles"] = annot["smiles"]

        merge_node(g, annot_node_label, annot_node_attrs)

        edge_attrs = Cons.PUBCHEM_GENE_COMPOUND_EDGE_ATTRS.copy()
        edge_attrs["assay_type"] = annot["assay_type"]
        edge_attrs["pubchem_assay_id"] = annot["pubchem_assay_id"]
        edge_attrs["outcome"] = annot["outcome"]
        edge_attrs["label"] = annot["outcome"]

        edge_hash = hash(frozenset(edge_attrs.items()))
        edge_attrs["edge_hash"] = edge_hash  # type: ignore
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
        edge_attrs = Cons.STRING_PPI_EDGE_ATTRS.copy()
        edge_attrs["score"] = ppi["score"]

        edge_hash = hash(frozenset(edge_attrs.items()))
        edge_attrs["edge_hash"] = edge_hash  # type: ignore
        edge_data = g.get_edge_data(gene_node_label, ppi[Cons.STRING_PPI_EDGE_MAIN_LABEL])

        edge_data = {} if edge_data is None else edge_data
        node_exists = [x for x, y in edge_data.items() if y["attr_dict"]["edge_hash"] == edge_hash]
        if len(node_exists) == 0 and not pd.isna(ppi[Cons.STRING_PPI_EDGE_MAIN_LABEL]):
            g.add_edge(
                gene_node_label,
                ppi[Cons.STRING_PPI_EDGE_MAIN_LABEL],
                label=Cons.STRING_PPI_EDGE_LABEL,
                attr_dict=edge_attrs,
            )
            g.add_edge(
                ppi[Cons.STRING_PPI_EDGE_MAIN_LABEL],
                gene_node_label,
                label=Cons.STRING_PPI_EDGE_LABEL,
                attr_dict=edge_attrs,
            )

    return g


def add_opentargets_disease_compound_subgraph(
    g: nx.MultiDiGraph, disease_node_label: str, annot_list: list
):
    """Construct part of the graph by linking the disease to compounds.

    :param g: the input graph to extend with new nodes and edges.
    :param disease_node_label: the disease node to be linked to compounds.
    :param annot_list: list of compounds from OpenTargets.
    :returns: a NetworkX MultiDiGraph
    """
    for annot in annot_list:
        if pd.isna(annot["relation"]):
            continue

        if pd.isna(annot[Cons.COMPOUND_NODE_MAIN_LABEL]):
            annot_node_label = annot["chembl_id"]
        else:
            annot_node_label = annot[Cons.COMPOUND_NODE_MAIN_LABEL]

        # create compound node and merge with existing node if it exists
        annot_node_attrs = Cons.OPENTARGETS_COMPOUND_NODE_ATTRS.copy()
        annot_node_attrs.update(
            {
                "name": annot["compound_name"],
                "id": annot_node_label,
                "datasource": Cons.OPENTARGETS,
            }
        )

        other_info = {
            "drugbank_id": annot["drugbank_id"],
            "compound_cid": annot["compound_cid"],
            "clincal_trial_phase": annot["clincal_trial_phase"],
            "is_approved": annot["is_approved"],
            "adverse_effect_count": annot["adverse_effect_count"],
        }

        for key, value in other_info.items():
            if not pd.isna(value):
                annot_node_attrs[key] = value

        merge_node(g, annot_node_label, annot_node_attrs)

        edge_attrs = Cons.OPENTARGETS_DISEASE_COMPOUND_EDGE_ATTRS.copy()
        edge_attrs["label"] = annot["relation"]
        edge_hash = hash(frozenset(edge_attrs.items()))
        edge_attrs["edge_hash"] = edge_hash  # type: ignore

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
                g, annot_node_label, annot[Cons.SIDE_EFFECT_NODE_MAIN_LABEL]
            )

    return g


def add_compoundwiki_subgraph(g, compound_node_label, annot_list):
    """
    Enrich an existing compound node in the graph with CompoundWiki annotations.

    :param g: the input NetworkX graph (MultiDiGraph).
    :param compound_node_label: the identifier (e.g., PubChem CID) of the compound node in the graph.
    :param annot_list: list of annotation dictionaries from CompoundWiki (one per compound).
    :return: the enriched NetworkX MultiDiGraph.
    """
    for annot in annot_list:
        target_cid = str(annot.get("target"))
        if not target_cid or target_cid != compound_node_label:
            continue

        annotations = annot.get("CompoundWiki_compounds", {})

        if g.has_node(compound_node_label):
            existing_attrs = g.nodes[compound_node_label].get("attr_dict", {})
            existing_attrs.update(annotations)
            existing_attrs["source"] = Cons.COMPOUNDWIKI
            g.nodes[compound_node_label]["attr_dict"] = existing_attrs
        else:
            node_attrs = Cons.COMPOUNDWIKI_NODE_ATTRS.copy()
            node_attrs.update(annotations)
            node_attrs["id"] = compound_node_label
            node_attrs["source"] = Cons.COMPOUNDWIKI
            g.merge_node(compound_node_label, attr_dict=node_attrs)

    return g


def add_wikipathways_molecular_subgraph(g, gene_node_label, annot_list):
    """Construct part of the graph by linking molecular entities from WP with MIMs.

    :param g: the input graph to extend with new nodes and edges.
    :param gene_node_label: the disease node to be linked to compounds.
    :param annot_list: result of querying WP for molecular interactions.
    :returns: a NetworkX MultiDiGraph
    """
    for annot in annot_list:
        for target_key in ["targetGene", "targetMetabolite"]:
            target = annot.get(target_key)
            target_node_label = None
            if target and target != gene_node_label:  # No interactions with self
                target_node_label = str(target)

            if target_node_label is not None:
                interaction_type = annot.get("mimtype", "Interaction")
                edge_attrs = Cons.MOLECULAR_INTERACTION_EDGE_ATTRS.copy()
                edge_attrs["interaction_type"] = interaction_type
                edge_attrs["rhea_id"] = annot.get("rhea_id", "")
                edge_attrs["pathway_id"] = annot.get("pathway_id", "")
                edge_attrs["edge_hash"] = hash(frozenset(edge_attrs.items()))  # type: ignore

                if not g.has_node(target_node_label):
                    node_attrs = Cons.MOLECULAR_PATHWAY_NODE_ATTRS.copy()
                    node_attrs.update(
                        {
                            "pathway_id": annot.get("pathway_id", ""),
                            "pathway_label": annot.get("pathway_label", ""),
                            "id": target_node_label,
                        }
                    )
                    g.add_node(target_node_label, attr_dict=node_attrs)

                edge_exists = False
                if g.has_edge(gene_node_label, target_node_label):
                    edge_data = g.get_edge_data(gene_node_label, target_node_label)
                    for edge_key in edge_data:
                        if (
                            edge_data[edge_key].get("attr_dict", {}).get("edge_hash")
                            == edge_attrs["edge_hash"]
                        ):
                            edge_exists = True
                            break

                if not edge_exists:
                    g.add_edge(
                        gene_node_label,
                        target_node_label,
                        label=interaction_type.upper(),
                        attr_dict=edge_attrs,
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


def add_gene_node(g, row, dea_columns):
    """Add gene node from each row of the combined_df to the graph.

    :param g: the input graph to extend with gene nodes.
    :param row: row in the combined DataFrame.
    :param dea_columns: list of dea_columns.
    :returns: label for gene node
    """
    gene_node_label = row["identifier"]
    gene_node_attrs = {
        "datasource": Cons.BRIDGEDB,
        "name": f"{row['identifier.source']}:{row['identifier']}",
        "id": row["target"],
        "labels": Cons.GENE_NODE_LABELS,
        row["target.source"]: row["target"],
    }

    for c in dea_columns:
        gene_node_attrs[c[:-4]] = row[c]

    g.add_node(gene_node_label, attr_dict=gene_node_attrs)
    return


def add_compound_node(g, row):
    """Add compound node from each row of the combined_df to the graph.

    :param g: the input graph to extend with compound nodes.
    :param row: row in the combined DataFrame.
    :returns: label for compound node
    """
    compound_node_label = row["identifier"]
    compound_node_attrs = {
        "datasource": Cons.BRIDGEDB,
        "name": f"{row['identifier.source']}:{row['identifier']}",
        "id": row["target"],
        "labels": Cons.COMPOUND_NODE_LABELS,
        row["target.source"]: row["target"],
    }

    g.add_node(compound_node_label, attr_dict=compound_node_attrs)
    return compound_node_label


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
    :param disease_compound: the input DataFrame containing disease-compound relationships.
    """
    for _i, row in disease_compound.iterrows():
        disease_node_label = row["identifier"].replace("_", ":")
        disease_annot_list = json.loads(json.dumps(row[Cons.OPENTARGETS_DISEASE_COMPOUND_COL]))

        if isinstance(disease_annot_list, float):
            disease_annot_list = []

        add_opentargets_disease_compound_subgraph(g, disease_node_label, disease_annot_list)


def process_ppi(g, gene_node_label, row):
    """Process protein-protein interactions and add them to the graph.

    :param g: the input graph to extend with gene nodes.
    :param gene_node_label: the gene node to be linked to annotation entities.
    :param row: row in the combined DataFrame.
    """
    if Cons.STRING_PPI_COL in row and row[Cons.STRING_PPI_COL] is not None:
        try:
            ppi_list = json.loads(json.dumps(row[Cons.STRING_PPI_COL]))
        except (ValueError, TypeError):
            ppi_list = []

        if isinstance(ppi_list, list) and len(ppi_list) > 0:
            valid_ppi_list = [
                item for item in ppi_list if pd.notna(item.get(Cons.STRING_PPI_EDGE_MAIN_LABEL))
            ]
            if valid_ppi_list:
                add_stringdb_ppi_subgraph(g, gene_node_label, valid_ppi_list)

        if not isinstance(ppi_list, float):
            add_stringdb_ppi_subgraph(g, gene_node_label, ppi_list)


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
                    annot_node_attrs = Cons.ENSEMBL_HOMOLOG_NODE_ATTRS.copy()
                    annot_node_attrs["id"] = homolog_node_label
                    annot_node_attrs["labels"] = Cons.HOMOLOG_NODE_LABELS
                    g.add_node(homolog_node_label, attr_dict=annot_node_attrs)

                    process_annotations(g, homolog_node_label, row, func_dict_hl)


def normalize_node_attributes(g):
    """Normalize node attributes by flattening the 'attr_dict'.

    :param g: the input graph to extend with gene nodes.
    """
    for node in g.nodes():
        if "attr_dict" not in g.nodes[node]:
            continue

        for k, v in g.nodes[node]["attr_dict"].items():
            if v is not None:
                g.nodes[node][k] = v

        del g.nodes[node]["attr_dict"]


def normalize_edge_attributes(g):
    """Normalize edge attributes by flattening the 'attr_dict'.

    :param g: the input graph to extend with gene nodes.
    """
    for u, v, k in g.edges(keys=True):
        if "attr_dict" not in g[u][v][k]:
            continue

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
    combined_df = combined_df[(combined_df["target.source"] == "Ensembl")]

    dea_columns = [c for c in combined_df.columns if c.endswith("_dea")]

    func_dict = {
        Cons.BGEE_GENE_EXPRESSION_LEVELS_COL: add_gene_bgee_subgraph,
        Cons.DISGENET_DISEASE_COL: add_disgenet_gene_disease_subgraph,
        Cons.LITERATURE_DISEASE_COL: add_literature_gene_disease_subgraph,
        Cons.MINERVA: add_minerva_gene_pathway_subgraph,
        Cons.WIKIPATHWAYS: add_wikipathways_gene_pathway_subgraph,
        Cons.KEGG_COL: add_kegg_gene_pathway_subgraph,
        Cons.OPENTARGETS_REACTOME_COL: add_opentargets_gene_reactome_pathway_subgraph,
        Cons.OPENTARGETS_GO_COL: add_opentargets_gene_go_subgraph,
        Cons.OPENTARGETS_GENE_COMPOUND_COL: add_opentargets_gene_compound_subgraph,
        Cons.MOLMEDB_PROTEIN_COMPOUND_COL: add_molmedb_gene_inhibitor_subgraph,
        Cons.PUBCHEM_COMPOUND_ASSAYS_COL: add_pubchem_assay_subgraph,
        Cons.WIKIPATHWAYS_MOLECULAR_COL: add_wikipathways_molecular_subgraph,
        Cons.ENSEMBL_HOMOLOG_COL: add_ensembl_homolog_subgraph,
    }

    if homolog_df_list is not None:
        process_homologs(g, combined_df, homolog_df_list, func_dict, dea_columns)

    else:
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


def _built_compound_based_graph(
    g: nx.MultiDiGraph,
    combined_df: pd.DataFrame,
    disease_compound=None,
    pathway_compound=None,
    homolog_df_list=None,
):
    """Build a gene-based graph."""
    # combined_df = combined_df[(combined_df["target.source"] == "Ensembl")]

    dea_columns = [c for c in combined_df.columns if c.endswith("_dea")]

    func_dict = {Cons.COMPOUNDWIKI_COL: add_compoundwiki_subgraph}

    if homolog_df_list is not None:
        process_homologs(g, combined_df, homolog_df_list, func_dict, dea_columns)

    else:
        for _i, row in tqdm(
            combined_df.iterrows(), total=combined_df.shape[0], desc="Building graph"
        ):
            if pd.isna(row["identifier"]) or pd.isna(row["target"]):
                continue
            compound_node_label = add_compound_node(g, row)
            process_annotations(g, compound_node_label, row, func_dict)
            process_ppi(g, compound_node_label, row)

    if disease_compound is not None:
        process_disease_compound(g, disease_compound)

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

    if main_target_type == "Ensembl":
        return _built_gene_based_graph(
            g, combined_df, disease_compound, pathway_compound, homolog_df_list
        )
    elif main_target_type == "PubChem Compound":
        return _built_compound_based_graph(
            g, combined_df, disease_compound, pathway_compound, homolog_df_list
        )
    else:
        raise ValueError(f"Unsupported target type: {main_target_type}")
