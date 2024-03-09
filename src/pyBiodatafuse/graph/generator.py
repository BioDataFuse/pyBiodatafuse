# coding: utf-8

"""Python module to construct a NetworkX graph from the annotated data frame."""

import json
import pickle

import networkx as nx
import pandas as pd


def load_dataframe_from_pickle(pickle_path: str) -> pd.DataFrame:
    """Load a previously annotated dataframe from a pickle file.

    :param pickle_path: the path to a previously obtained annotation dataframe dumped as a pickle file.
    :returns: a Pandas dataframe.
    """
    with open(pickle_path, "rb") as rin:
        df = pickle.load(rin)

    return df


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
                "source": "DisGeNET",
                "labels": dg["disease_name"],
                "id": dg["diseaseid"],
                "node_type": "disease",
                "disease_id": dg["diseaseid"],
                "disease_class": dg["disease_class"],
                "disease_class_name": dg["disease_class_name"],
                "disease_type": dg["disease_type"],
                "disease_semantic_type": dg["disease_semantic_type"],
            }

            g.add_node(dg_node_label, attr_dict=dg_node_attrs)

            edge_attrs = {
                "source": "DisGeNET",
                "label": "associated_with",
                "score": dg["score"],
                "year_initial": dg["year_initial"] if not pd.isna(dg["year_initial"]) else "",
                "year_final": dg["year_final"] if not pd.isna(dg["year_final"]) else "",
                "ei": dg["ei"] if not pd.isna(dg["ei"]) else "",
                "el": dg["el"] if not pd.isna(dg["el"]) else "",
            }

            edge_hash = hash(frozenset(edge_attrs.items()))
            edge_attrs["edge_hash"] = edge_hash
            edge_data = g.get_edge_data(gene_node_label, dg_node_label)
            edge_data = {} if edge_data is None else edge_data
            node_exists = [
                x for x, y in edge_data.items() if y["attr_dict"]["edge_hash"] == edge_hash
            ]

            if len(node_exists) == 0:
                g.add_edge(gene_node_label, dg_node_label, attr_dict=edge_attrs)

    return g


def add_opentargets_location_subgraph(g, gene_node_label, annot_list):
    """Construct part of the graph by linking the gene to a list of annotation entities (disease, drug ..etc).

    :param g: the input graph to extend with new nodes and edges.
    :param gene_node_label: the gene node to be linked to annotation entities.
    :param annot_list: list of annotations from a specific source (e.g. DisGeNET, WikiPathways ..etc).
    :returns: a NetworkX MultiDiGraph
    """
    for loc in annot_list:
        if not pd.isna(loc["location"]) and not pd.isna(loc["subcellular_loc"]):
            loc_node_label = loc["location"]
            loc_node_attrs = {
                "source": "OpenTargets",
                "labels": loc["location"],
                "id": loc["loc_identifier"],
                "node_type": "location",
                "subcellular_loc": loc["subcellular_loc"],
            }

            g.add_node(loc_node_label, attr_dict=loc_node_attrs)

            edge_attrs = {"source": "OpenTargets", "label": "localized_in"}

            edge_hash = hash(frozenset(edge_attrs.items()))
            edge_attrs["edge_hash"] = edge_hash
            edge_data = g.get_edge_data(gene_node_label, loc_node_label)
            edge_data = {} if edge_data is None else edge_data
            node_exists = [
                x for x, y in edge_data.items() if y["attr_dict"]["edge_hash"] == edge_hash
            ]

            if len(node_exists) == 0:
                g.add_edge(gene_node_label, loc_node_label, attr_dict=edge_attrs)

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
            "source": "OpenTargets",
            "labels": go["go_name"],
            "id": go["go_id"],
            "node_type": "gene ontology",
        }

        g.add_node(go_node_label, attr_dict=go_node_attrs)

        edge_attrs = {"source": "OpenTargets", "label": "part_of_go"}

        edge_hash = hash(frozenset(edge_attrs.items()))
        edge_attrs["edge_hash"] = edge_hash
        edge_data = g.get_edge_data(gene_node_label, go_node_label)
        edge_data = {} if edge_data is None else edge_data
        node_exists = [x for x, y in edge_data.items() if y["attr_dict"]["edge_hash"] == edge_hash]

        if len(node_exists) == 0:
            g.add_edge(gene_node_label, go_node_label, attr_dict=edge_attrs)

    return g


def add_opentargets_pathway_subgraph(g, gene_node_label, annot_list):
    """Construct part of the graph by linking the gene to a list of annotation entities (disease, drug ..etc).

    :param g: the input graph to extend with new nodes and edges.
    :param gene_node_label: the gene node to be linked to annotation entities.
    :param annot_list: list of annotations from a specific source (e.g. DisGeNET, WikiPathways ..etc).
    :returns: a NetworkX MultiDiGraph
    """
    for pathway in annot_list:
        if not pd.isna(pathway["pathway_id"]):
            pathway_node_label = pathway["pathway_name"]
            pathway_node_attrs = {
                "source": "OpenTargets",
                "labels": pathway["pathway_name"],
                "id": pathway["pathway_id"],
                "node_type": "reactome pathways",
            }

            g.add_node(pathway_node_label, attr_dict=pathway_node_attrs)

            edge_attrs = {"source": "OpenTargets", "label": "part_of_pathway"}

            edge_hash = hash(frozenset(edge_attrs.items()))
            edge_attrs["edge_hash"] = edge_hash
            edge_data = g.get_edge_data(gene_node_label, pathway_node_label)
            edge_data = {} if edge_data is None else edge_data
            node_exists = [
                x for x, y in edge_data.items() if y["attr_dict"]["edge_hash"] == edge_hash
            ]

            if len(node_exists) == 0:
                g.add_edge(gene_node_label, pathway_node_label, attr_dict=edge_attrs)

    return g


def add_opentargets_drug_subgraph(g, gene_node_label, annot_list):
    """Construct part of the graph by linking the gene to a list of annotation entities (disease, drug ..etc).

    :param g: the input graph to extend with new nodes and edges.
    :param gene_node_label: the gene node to be linked to annotation entities.
    :param annot_list: list of annotations from a specific source (e.g. DisGeNET, WikiPathways ..etc).
    :returns: a NetworkX MultiDiGraph
    """
    for drug in annot_list:
        if not pd.isna(drug["relation"]):
            drug_node_label = drug["drug_name"]
            drug_node_attrs = {
                "source": "OpenTargets",
                "labels": drug["drug_name"],
                "id": drug["chembl_id"],
                "node_type": "drug interactions",
            }

            g.add_node(drug_node_label, attr_dict=drug_node_attrs)

            edge_attrs = {"source": "OpenTargets", "label": drug["relation"]}

            edge_hash = hash(frozenset(edge_attrs.items()))
            edge_attrs["edge_hash"] = edge_hash
            edge_data = g.get_edge_data(drug_node_label, gene_node_label)
            edge_data = {} if edge_data is None else edge_data
            node_exists = [
                x for x, y in edge_data.items() if y["attr_dict"]["edge_hash"] == edge_hash
            ]

            if len(node_exists) == 0:
                g.add_edge(drug_node_label, gene_node_label, attr_dict=edge_attrs)

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
                "source": "OpenTargets",
                "labels": dg["disease_name"],
                "id": dg["disease_id"],
                "node_type": "disease",
                "therapeutic_areas": dg["therapeutic_areas"],
            }

            g.add_node(dg_node_label, attr_dict=dg_node_attrs)

            edge_attrs = {"source": "OpenTargets", "label": "associated_with"}

            edge_hash = hash(frozenset(edge_attrs.items()))
            edge_attrs["edge_hash"] = edge_hash
            edge_data = g.get_edge_data(gene_node_label, dg_node_label)
            edge_data = {} if edge_data is None else edge_data
            node_exists = [
                x for x, y in edge_data.items() if y["attr_dict"]["edge_hash"] == edge_hash
            ]

            if len(node_exists) == 0:
                g.add_edge(gene_node_label, dg_node_label, attr_dict=edge_attrs)

    return g


def add_wikipathways_subgraph(g, gene_node_label, annot_list):
    """Construct part of the graph by linking the gene to a list of annotation entities (disease, drug ..etc).

    :param g: the input graph to extend with new nodes and edges.
    :param gene_node_label: the gene node to be linked to annotation entities.
    :param annot_list: list of annotations from a specific source (e.g. DisGeNET, WikiPathways ..etc).
    :returns: a NetworkX MultiDiGraph
    """
    for pathway in annot_list:
        if not pd.isna(pathway["pathwayLabel"]):
            pathway_node_label = pathway["pathwayLabel"]
            pathway_node_attrs = {
                "source": "WikiPathways",
                "labels": pathway["pathwayLabel"],
                "id": pathway["pathwayId"],
                "node_type": "wikipathways pathway",
                "gene_count": pathway["pathwayGeneCount"],
            }

            g.add_node(pathway_node_label, attr_dict=pathway_node_attrs)

            edge_attrs = {"source": "WikiPathways", "label": "part_of_pathway"}

            edge_hash = hash(frozenset(edge_attrs.items()))
            edge_attrs["edge_hash"] = edge_hash
            edge_data = g.get_edge_data(gene_node_label, pathway_node_label)
            edge_data = {} if edge_data is None else edge_data
            node_exists = [
                x for x, y in edge_data.items() if y["attr_dict"]["edge_hash"] == edge_hash
            ]

            if len(node_exists) == 0:
                g.add_edge(gene_node_label, pathway_node_label, attr_dict=edge_attrs)

    return g


def add_ppi_subgraph(g, gene_node_label, annot_list):
    """Construct part of the graph by linking the gene to a list of annotation entities (disease, drug ..etc).

    :param g: the input graph to extend with new nodes and edges.
    :param gene_node_label: the gene node to be linked to annotation entities.
    :param annot_list: list of annotations from a specific source (e.g. DisGeNET, WikiPathways ..etc).
    :returns: a NetworkX MultiDiGraph
    """
    for ppi in annot_list:
        edge_attrs = {"source": "STRING", "label": "interacts_with", "score": ppi["score"]}

        edge_hash = hash(frozenset(edge_attrs.items()))
        edge_attrs["edge_hash"] = edge_hash
        edge_data = g.get_edge_data(gene_node_label, ppi["stringdb_link_to"])
        edge_data = {} if edge_data is None else edge_data
        node_exists = [x for x, y in edge_data.items() if y["attr_dict"]["edge_hash"] == edge_hash]

        if len(node_exists) == 0:
            g.add_edge(gene_node_label, ppi["stringdb_link_to"], attr_dict=edge_attrs)

    return g


def add_gene_inhibitor(g, gene_node_label, annot_list):
    """Construct part of the graph by linking the gene to a list of annotation entities (disease, drug ..etc).

    :param g: the input graph to extend with new nodes and edges.
    :param gene_node_label: the gene node to be linked to annotation entities.
    :param annot_list: list of annotations from a specific source (e.g. DisGeNET, WikiPathways ..etc).
    :returns: a NetworkX MultiDiGraph
    """
    for inhibitor in annot_list:
        if not pd.isna(inhibitor["InChIKey"]):
            inhibitor_node_label = inhibitor["label"]
            inhibitor_node_attrs = {
                "source": "MolMeDB",
                "labels": inhibitor["label"],
                "InChIKey": inhibitor["InChIKey"],
                "MolMeDB_id": inhibitor["molmedb_id"],
                "node_type": "inhibitor",
            }

            if not pd.isna(inhibitor["SMILES"]):
                inhibitor_node_attrs["SMILES"] = inhibitor["SMILES"]
            if not pd.isna(inhibitor["pubchem_compound_id"]):
                inhibitor_node_attrs["PubChem_compound_id"] = inhibitor["pubchem_compound_id"]
            if not pd.isna(inhibitor["chebi_id"]):
                inhibitor_node_attrs["ChEBI_id"] = inhibitor["chebi_id"]
            if not pd.isna(inhibitor["drugbank_id"]):
                inhibitor_node_attrs["DrugBank_id"] = inhibitor["drugbank_id"]

            g.add_node(inhibitor_node_label, attr_dict=inhibitor_node_attrs)

            edge_attrs = {
                "source": "MolMeDB",
                "label": "is_inhibited_by",
            }

            if not pd.isna(inhibitor["source_doi"]):
                inhibitor_node_attrs["reference_doi"] = inhibitor["source_doi"]
            if not pd.isna(inhibitor["source_pmid"]):
                inhibitor_node_attrs["reference_pmid"] = inhibitor["source_pmid"]

            edge_hash = hash(frozenset(edge_attrs.items()))
            edge_attrs["edge_hash"] = edge_hash
            edge_data = g.get_edge_data(gene_node_label, inhibitor_node_label)
            edge_data = {} if edge_data is None else edge_data
            node_exists = [
                x for x, y in edge_data.items() if y["attr_dict"]["edge_hash"] == edge_hash
            ]

            if len(node_exists) == 0:
                g.add_edge(gene_node_label, inhibitor_node_label, attr_dict=edge_attrs)

    return g


def generate_networkx_graph(fuse_df: pd.DataFrame):
    """Construct a NetWorkX graph from a Pandas DataFrame of genes and their multi-source annotations.

    :param fuse_df: the input dataframe to be converted into a graph.
    :returns: a NetworkX MultiDiGraph
    """
    g = nx.MultiDiGraph()

    dea_columns = [c for c in fuse_df.columns if c.endswith("_dea")]

    func_dict = {
        "DisGeNET": add_disgenet_disease_subgraph,
        "OpenTargets_Location": add_opentargets_location_subgraph,
        "GO_Process": add_opentargets_go_subgraph,
        "Reactome_Pathways": add_opentargets_pathway_subgraph,
        "ChEMBL_Drugs": add_opentargets_drug_subgraph,
        "OpenTargets_Diseases": add_opentargets_disease_subgraph,
        "WikiPathways": add_wikipathways_subgraph,
        "transporter_inhibitor": add_gene_inhibitor,
    }

    for _i, row in fuse_df.iterrows():
        gene_node_label = row["identifier"]
        gene_node_attrs = {
            "source": "BridgeDB",
            "labels": row["identifier"],
            "id": row["target"],
            "node_type": "gene",
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

    if "stringdb" in row:
        for _i, row in fuse_df.iterrows():
            ppi_list = json.loads(json.dumps(row["stringdb"]))

            if ppi_list is None:
                ppi_list = []

            add_ppi_subgraph(g, gene_node_label, ppi_list)

    for node in g.nodes():
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
