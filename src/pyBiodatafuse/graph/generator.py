# coding: utf-8

"""Python module to construct a NetworkX graph from the annotated data frame."""

import json
import pickle

import networkx as nx
import pandas as pd

from pyBiodatafuse.constants import (
    BGEE,
    BGEE_EDGE_LABEL,
    BGEE_NODE_LABELS,
    DISGENET,
    DISGENET_EDGE_LABEL,
    DISGENET_NODE_LABELS,
    MINERVA,
    MINERVA_EDGE_LABEL,
    MINERVA_NODE_LABELS,
    MOLMEDB,
    MOLMEDB_EDGE_LABEL,
    MOLMEDB_INHIBITOR_COL,
    MOLMEDB_NODE_LABELS,
    OPENTARGETS,
    OPENTARGETS_COMPOUND_COL,
    OPENTARGETS_COMPOUND_NODE_LABELS,
    OPENTARGETS_DISEASE_COL,
    OPENTARGETS_DISEASE_EDGE_LABEL,
    OPENTARGETS_DISEASE_NODE_LABELS,
    OPENTARGETS_GO_COL,
    OPENTARGETS_GO_EDGE_LABEL,
    OPENTARGETS_GO_NODE_LABELS,
    OPENTARGETS_LOCATION_COL,
    OPENTARGETS_LOCATION_EDGE_LABEL,
    OPENTARGETS_LOCATION_NODE_LABELS,
    OPENTARGETS_REACTOME_COL,
    OPENTARGETS_REACTOME_EDGE_LABEL,
    OPENTARGETS_REACTOME_NODE_LABELS,
    PUBCHEM,
    PUBCHEM_NODE_LABELS,
    STRING,
    STRING_EDGE_LABEL,
    WIKIPATHWAYS,
    WIKIPATHWAYS_EDGE_LABEL,
    WIKIPATHWAYS_NODE_LABELS,
)


def load_dataframe_from_pickle(pickle_path: str) -> pd.DataFrame:
    """Load a previously annotated dataframe from a pickle file.

    :param pickle_path: the path to a previously obtained annotation dataframe dumped as a pickle file.
    :returns: a Pandas dataframe.
    """
    with open(pickle_path, "rb") as rin:
        df = pickle.load(rin)

    return df


def add_bgee_subgraph(g, gene_node_label, annot_list):
    """Construct part of the graph by linking the gene to a list of annotation entities (disease, drug ..etc).

    :param g: the input graph to extend with new nodes and edges.
    :param gene_node_label: the gene node to be linked to annotation entities.
    :param annot_list: list of annotations from a specific source (e.g. DisGeNET, WikiPathways ..etc).
    :returns: a NetworkX MultiDiGraph
    """
    for annot in annot_list:
        if not pd.isna(annot["anatomical_entity_name"]):
            annot_node_label = annot["anatomical_entity_name"]
            annot_node_attrs = {
                "source": BGEE,
                "name": annot["anatomical_entity_name"],
                "id": annot["anatomical_entity_id"],
                "labels": BGEE_NODE_LABELS,
            }

            if not pd.isna(annot["developmental_stage_id"]):
                annot_node_attrs["developmental_stage_id"] = annot["developmental_stage_id"]
            if not pd.isna(annot["developmental_stage_name"]):
                annot_node_attrs["developmental_stage_name"] = annot["developmental_stage_name"]
            if not pd.isna(annot["expression_level"]):
                annot_node_attrs["expression_level"] = annot["expression_level"]
            if not pd.isna(annot["confidence_level_id"]):
                annot_node_attrs["confidence_level_id"] = annot["confidence_level_id"]

            g.add_node(annot_node_label, attr_dict=annot_node_attrs)

            edge_attrs = {
                "source": BGEE,
                "label": BGEE_EDGE_LABEL,
            }

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
    for dg in annot_list:
        if not pd.isna(dg["disease_name"]):
            dg_node_label = dg["disease_name"]
            dg_node_attrs = {
                "source": DISGENET,
                "name": dg["disease_name"],
                "id": dg["disease_id"],
                "labels": DISGENET_NODE_LABELS,
                "disease_source": dg["source"],
            }

            g.add_node(dg_node_label, attr_dict=dg_node_attrs)

            edge_attrs = {
                "source": DISGENET,
                "label": DISGENET_EDGE_LABEL,
                "score": dg["score"],
            }

            edge_hash = hash(frozenset(edge_attrs.items()))
            edge_attrs["edge_hash"] = edge_hash
            edge_data = g.get_edge_data(gene_node_label, dg_node_label)
            edge_data = {} if edge_data is None else edge_data
            node_exists = [
                x for x, y in edge_data.items() if y["attr_dict"]["edge_hash"] == edge_hash
            ]

            if len(node_exists) == 0:
                g.add_edge(
                    gene_node_label, dg_node_label, label=DISGENET_EDGE_LABEL, attr_dict=edge_attrs
                )

    return g


def add_minerva_subgraph(g, gene_node_label, annot_list):
    """Construct part of the graph by linking the gene to a list of annotation entities (disease, drug ..etc).

    :param g: the input graph to extend with new nodes and edges.
    :param gene_node_label: the gene node to be linked to annotation entities.
    :param annot_list: list of annotations from a specific source (e.g. DisGeNET, WikiPathways ..etc).
    :returns: a NetworkX MultiDiGraph
    """
    for pathway in annot_list:
        if not pd.isna(pathway["pathway_label"]):
            pathway_node_label = pathway["pathway_label"]
            pathway_node_attrs = {
                "source": MINERVA,
                "name": pathway["pathway_label"],
                "id": pathway["pathway_id"],
                "labels": MINERVA_NODE_LABELS,
                "gene_count": pathway["pathway_gene_count"],
            }

            g.add_node(pathway_node_label, attr_dict=pathway_node_attrs)

            edge_attrs = {"source": MINERVA, "label": MINERVA_EDGE_LABEL}

            edge_hash = hash(frozenset(edge_attrs.items()))
            edge_attrs["edge_hash"] = edge_hash
            edge_data = g.get_edge_data(gene_node_label, pathway_node_label)
            edge_data = {} if edge_data is None else edge_data
            node_exists = [
                x for x, y in edge_data.items() if y["attr_dict"]["edge_hash"] == edge_hash
            ]

            if len(node_exists) == 0:
                g.add_edge(
                    gene_node_label,
                    pathway_node_label,
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
    for pathway in annot_list:
        if not pd.isna(pathway["pathway_label"]):
            pathway_node_label = pathway["pathway_label"]
            pathway_node_attrs = {
                "source": WIKIPATHWAYS,
                "name": pathway["pathway_label"],
                "id": pathway["pathway_id"],
                "labels": WIKIPATHWAYS_NODE_LABELS,
                "gene_count": pathway["pathway_gene_count"],
            }

            g.add_node(pathway_node_label, attr_dict=pathway_node_attrs)

            edge_attrs = {"source": WIKIPATHWAYS, "label": WIKIPATHWAYS_EDGE_LABEL}

            edge_hash = hash(frozenset(edge_attrs.items()))
            edge_attrs["edge_hash"] = edge_hash
            edge_data = g.get_edge_data(gene_node_label, pathway_node_label)
            edge_data = {} if edge_data is None else edge_data
            node_exists = [
                x for x, y in edge_data.items() if y["attr_dict"]["edge_hash"] == edge_hash
            ]

            if len(node_exists) == 0:
                g.add_edge(
                    gene_node_label,
                    pathway_node_label,
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
    for pathway in annot_list:
        if not pd.isna(pathway["pathway_id"]):
            pathway_node_label = pathway["pathway_label"]
            pathway_node_attrs = {
                "source": OPENTARGETS,
                "name": pathway["pathway_label"],
                "id": pathway["pathway_id"],
                "labels": OPENTARGETS_REACTOME_NODE_LABELS,
            }

            g.add_node(pathway_node_label, attr_dict=pathway_node_attrs)

            edge_attrs = {"source": OPENTARGETS, "label": OPENTARGETS_REACTOME_EDGE_LABEL}

            edge_hash = hash(frozenset(edge_attrs.items()))
            edge_attrs["edge_hash"] = edge_hash
            edge_data = g.get_edge_data(gene_node_label, pathway_node_label)
            edge_data = {} if edge_data is None else edge_data
            node_exists = [
                x for x, y in edge_data.items() if y["attr_dict"]["edge_hash"] == edge_hash
            ]

            if len(node_exists) == 0:
                g.add_edge(
                    gene_node_label,
                    pathway_node_label,
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
    for go in annot_list:
        go_node_label = go["go_name"]
        go_node_attrs = {
            "source": OPENTARGETS,
            "name": go["go_name"],
            "id": go["go_id"],
            "labels": OPENTARGETS_GO_NODE_LABELS,
        }

        g.add_node(go_node_label, attr_dict=go_node_attrs)

        edge_attrs = {"source": OPENTARGETS, "label": OPENTARGETS_GO_EDGE_LABEL}

        edge_hash = hash(frozenset(edge_attrs.items()))
        edge_attrs["edge_hash"] = edge_hash
        edge_data = g.get_edge_data(gene_node_label, go_node_label)
        edge_data = {} if edge_data is None else edge_data
        node_exists = [x for x, y in edge_data.items() if y["attr_dict"]["edge_hash"] == edge_hash]

        if len(node_exists) == 0:
            g.add_edge(
                gene_node_label,
                go_node_label,
                label=OPENTARGETS_GO_EDGE_LABEL,
                attr_dict=edge_attrs,
            )

    return g


def add_opentargets_location_subgraph(g, gene_node_label, annot_list):
    """Construct part of the graph by linking the gene to a list of annotation entities (disease, drug ..etc).

    :param g: the input graph to extend with new nodes and edges.
    :param gene_node_label: the gene node to be linked to annotation entities.
    :param annot_list: list of annotations from a specific source (e.g. DisGeNET, WikiPathways ..etc).
    :returns: a NetworkX MultiDiGraph
    """
    for loc in annot_list:
        if not pd.isna(loc["location"]) and not pd.isna(loc["subcellular_location"]):
            loc_node_label = loc["location"]
            loc_node_attrs = {
                "source": OPENTARGETS,
                "name": loc["location"],
                "id": loc["location_id"],
                "subcellular_location": loc["subcellular_location"],
                "labels": OPENTARGETS_LOCATION_NODE_LABELS,
            }

            g.add_node(loc_node_label, attr_dict=loc_node_attrs)

            edge_attrs = {"source": OPENTARGETS, "label": OPENTARGETS_LOCATION_EDGE_LABEL}

            edge_hash = hash(frozenset(edge_attrs.items()))
            edge_attrs["edge_hash"] = edge_hash
            edge_data = g.get_edge_data(gene_node_label, loc_node_label)
            edge_data = {} if edge_data is None else edge_data
            node_exists = [
                x for x, y in edge_data.items() if y["attr_dict"]["edge_hash"] == edge_hash
            ]

            if len(node_exists) == 0:
                g.add_edge(
                    gene_node_label,
                    loc_node_label,
                    label=OPENTARGETS_LOCATION_EDGE_LABEL,
                    attr_dict=edge_attrs,
                )

    return g


def add_opentargets_disease_subgraph(g, gene_node_label, annot_list):
    """Construct part of the graph by linking the gene to a list of annotation entities (disease, drug ..etc).

    :param g: the input graph to extend with new nodes and edges.
    :param gene_node_label: the gene node to be linked to annotation entities.
    :param annot_list: list of annotations from a specific source (e.g. DisGeNET, WikiPathways ..etc).
    :returns: a NetworkX MultiDiGraph
    """
    for dg in annot_list:
        if not pd.isna(dg["disease_name"]):
            dg_node_label = dg["disease_name"]
            dg_node_attrs = {
                "source": OPENTARGETS,
                "name": dg["disease_name"],
                "id": dg["disease_id"],
                "labels": OPENTARGETS_DISEASE_NODE_LABELS,
                "therapeutic_areas": dg["therapeutic_areas"],
            }

            g.add_node(dg_node_label, attr_dict=dg_node_attrs)

            edge_attrs = {"source": OPENTARGETS, "label": OPENTARGETS_DISEASE_EDGE_LABEL}

            edge_hash = hash(frozenset(edge_attrs.items()))
            edge_attrs["edge_hash"] = edge_hash
            edge_data = g.get_edge_data(gene_node_label, dg_node_label)
            edge_data = {} if edge_data is None else edge_data
            node_exists = [
                x for x, y in edge_data.items() if y["attr_dict"]["edge_hash"] == edge_hash
            ]

            if len(node_exists) == 0:
                g.add_edge(
                    gene_node_label,
                    dg_node_label,
                    label=OPENTARGETS_DISEASE_EDGE_LABEL,
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
    for compound in annot_list:
        if not pd.isna(compound["relation"]):
            compound_node_label = compound["compound_name"]
            compound_node_attrs = {
                "source": OPENTARGETS,
                "name": compound["compound_name"],
                "id": compound["chembl_id"],
                "is_approved": compound["is_approved"],
                "labels": OPENTARGETS_COMPOUND_NODE_LABELS,
            }

            g.add_node(compound_node_label, attr_dict=compound_node_attrs)

            edge_attrs = {"source": OPENTARGETS, "label": compound["relation"]}

            edge_hash = hash(frozenset(edge_attrs.items()))
            edge_attrs["edge_hash"] = edge_hash
            edge_data = g.get_edge_data(compound_node_label, gene_node_label)
            edge_data = {} if edge_data is None else edge_data
            node_exists = [
                x for x, y in edge_data.items() if y["attr_dict"]["edge_hash"] == edge_hash
            ]

            if len(node_exists) == 0:
                g.add_edge(
                    compound_node_label,
                    gene_node_label,
                    label=compound["relation"],
                    attr_dict=edge_attrs,
                )

    return g


def add_molmedb_gene_inhibitor(g, gene_node_label, annot_list):
    """Construct part of the graph by linking the gene to a list of annotation entities (disease, compound ..etc).

    :param g: the input graph to extend with new nodes and edges.
    :param gene_node_label: the gene node to be linked to annotation entities.
    :param annot_list: list of annotations from a specific source (e.g. DisGeNET, WikiPathways ..etc).
    :returns: a NetworkX MultiDiGraph
    """
    for inhibitor in annot_list:
        if not pd.isna(inhibitor["InChIKey"]):
            inhibitor_node_label = inhibitor["compound_name"]
            inhibitor_node_attrs = {
                "source": MOLMEDB,
                "id": inhibitor["compound_cid"],
                "name": inhibitor["compound_name"],
                "InChIKey": inhibitor["InChIKey"],
                "MolMeDB_id": inhibitor["molmedb_id"],
                "labels": MOLMEDB_NODE_LABELS,
            }

            if not pd.isna(inhibitor["SMILES"]):
                inhibitor_node_attrs["SMILES"] = inhibitor["SMILES"]
            if not pd.isna(inhibitor["compound_cid"]):
                inhibitor_node_attrs["compound_cid"] = inhibitor["compound_cid"]
            if not pd.isna(inhibitor["chebi_id"]):
                inhibitor_node_attrs["ChEBI_id"] = inhibitor["chebi_id"]
            if not pd.isna(inhibitor["drugbank_id"]):
                inhibitor_node_attrs["DrugBank_id"] = inhibitor["drugbank_id"]
            if not pd.isna(inhibitor["source_doi"]):
                inhibitor_node_attrs["source_doi"] = inhibitor["source_doi"]
            if not pd.isna(inhibitor["source_pmid"]):
                inhibitor_node_attrs["source_pmid"] = inhibitor["source_pmid"]
            if not pd.isna(inhibitor["pdb_ligand_id"]):
                inhibitor_node_attrs["pdb_ligand_id"] = inhibitor["pdb_ligand_id"]

            g.add_node(inhibitor_node_label, attr_dict=inhibitor_node_attrs)

            edge_attrs = {
                "source": MOLMEDB,
                "label": MOLMEDB_EDGE_LABEL,
            }

            edge_hash = hash(frozenset(edge_attrs.items()))
            edge_attrs["edge_hash"] = edge_hash
            edge_data = g.get_edge_data(gene_node_label, inhibitor_node_label)
            edge_data = {} if edge_data is None else edge_data
            node_exists = [
                x for x, y in edge_data.items() if y["attr_dict"]["edge_hash"] == edge_hash
            ]

            if len(node_exists) == 0:
                g.add_edge(
                    gene_node_label,
                    inhibitor_node_label,
                    label=MOLMEDB_EDGE_LABEL,
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
    for assay in annot_list:
        if not pd.isna(assay["InChIKey"]):
            inhibitor_node_label = assay["compound_name"]
            inhibitor_node_attrs = {
                "source": PUBCHEM,
                "id": assay["compound_cid"],
                "name": assay["compound_name"],
                "InChI": assay["InChI"],
                "labels": PUBCHEM_NODE_LABELS,
            }

            if not pd.isna(assay["SMILES"]):
                inhibitor_node_attrs["SMILES"] = assay["SMILES"]

            g.add_node(inhibitor_node_label, attr_dict=inhibitor_node_attrs)

            edge_attrs = {
                "source": "PubChem",
                "assay_type": assay["assay_type"],
                "label": assay["outcome"],
            }

            edge_hash = hash(frozenset(edge_attrs.items()))
            edge_attrs["edge_hash"] = edge_hash
            edge_data = g.get_edge_data(gene_node_label, inhibitor_node_label)
            edge_data = {} if edge_data is None else edge_data
            node_exists = [
                x for x, y in edge_data.items() if y["attr_dict"]["edge_hash"] == edge_hash
            ]

            if len(node_exists) == 0:
                g.add_edge(
                    gene_node_label,
                    inhibitor_node_label,
                    label=assay["outcome"],
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
        edge_attrs = {"source": STRING, "label": STRING_EDGE_LABEL, "score": ppi.get("score")}

        edge_hash = hash(frozenset(edge_attrs.items()))
        edge_attrs["edge_hash"] = edge_hash
        edge_data = g.get_edge_data(gene_node_label, ppi["stringdb_link_to"])
        edge_data = {} if edge_data is None else edge_data
        node_exists = [x for x, y in edge_data.items() if y["attr_dict"]["edge_hash"] == edge_hash]

        if len(node_exists) == 0:
            g.add_edge(
                gene_node_label,
                ppi["stringdb_link_to"],
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
        DISGENET: add_disgenet_disease_subgraph,
        OPENTARGETS_LOCATION_COL: add_opentargets_location_subgraph,
        OPENTARGETS_GO_COL: add_opentargets_go_subgraph,
        OPENTARGETS_REACTOME_COL: add_opentargets_reactome_pathway_subgraph,
        OPENTARGETS_COMPOUND_COL: add_opentargets_compound_subgraph,
        OPENTARGETS_DISEASE_COL: add_opentargets_disease_subgraph,
        WIKIPATHWAYS: add_wikipathways_subgraph,
        MINERVA: add_minerva_subgraph,
        MOLMEDB_INHIBITOR_COL: add_molmedb_gene_inhibitor,
        BGEE: add_bgee_subgraph,
    }

    for _i, row in fuse_df.iterrows():
        gene_node_label = row["identifier"]
        gene_node_attrs = {
            "source": "BridgeDB",
            "name": row["identifier"],
            "id": row["target"],
            "labels": "Gene",
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
        for _i, row in fuse_df.iterrows():
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
