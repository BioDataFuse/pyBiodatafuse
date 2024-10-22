"""Codes taken from DreamWalk repository: https://github.com/eugenebang/DREAMwalk"""

import pandas as pd
from tqdm import tqdm
import networkx as nx
from pyBiodatafuse.analyzer.summarize import BioGraph


def get_drug_disease_file(graph, output_dir: str):
    """Generate DDA files."""
    drug_nodes = [node for node, label in graph.nodes(data="labels") if label == "Compound"]

    disease_nodes = [node for node, label in graph.nodes(data="labels") if label == "Disease"]

    drug_disease_edges = []

    for drug in tqdm(drug_nodes):
        for disease in disease_nodes:
            if graph.has_edge(drug, disease):
                drug_disease_edges.append({"drug": drug, "disease": disease, "label": 1})
            else:
                drug_disease_edges.append({"drug": drug, "disease": disease, "label": 0})

    tmp = pd.DataFrame(drug_disease_edges)

    known_dda = tmp[tmp["label"] == 1]
    unknown_dda = tmp[tmp["label"] == 0]

    print(f"Known DDA: {len(known_dda)}, Unknown DDA: {len(unknown_dda)}")

    # Sampling same number of negative samples as positive samples
    for i in range(1, 11):
        sampled = unknown_dda.sample(n=len(known_dda))
        sampled.to_csv(f"{output_dir}/dda_{i}.tsv", sep="\t", index=False)


def create_files(graph_obj: BioGraph, output_dir: str = "./dreamwalk_data"):
    """Generate files for the DREAMwalk algorithm from the BioGraph object."""
    graph = graph_obj.graph

    # Node type label file - TSV files with two columns: node_id, node_type
    nodes_of_interest = [
        node
        for node, label in graph.nodes(data="labels")
        if label in ["Gene", "Disease", "Compound"]
    ]

    subgraph = graph_obj.get_subgraph(nodes_of_interest)
    print(f"Subgraph nodes: {len(subgraph.nodes())}, edges: {len(subgraph.edges())}")

    nx.write_gml(subgraph, f"{output_dir}/subgraph_graph.gml")

    tmp = []
    for node, node_label in subgraph.nodes(data="labels"):
        tmp.append({"node": node, "type": node_label})

    tmp_df = pd.DataFrame(tmp)
    tmp_df["type"] = tmp_df["type"].map({"Gene": "gene", "Disease": "disease", "Compound": "drug"})
    tmp_df = tmp_df.sort_values(by="type")
    tmp_df.to_csv(f"{output_dir}/nodetypes.tsv", sep="\t", index=False)

    # Network file - TXT file with four columns: source, target, edgetype, weight
    rel_to_id = {"activates": 1, "inhibits": 1, "associated_with": 2, "interacts_with": 3}

    graph_data = []

    for source, target, edge in subgraph.edges(data=True):
        if edge["label"] not in rel_to_id:
            continue

        edge_id = rel_to_id[edge["label"]]

        graph_data.append(
            {
                "source": source,
                "target": target,
                "edgetype": edge_id,
                "weight": 1,
            }
        )

    output_graph = pd.DataFrame(graph_data)
    output_graph["edge_counter"] = range(1, len(output_graph) + 1)
    output_graph.to_csv(f"{output_dir}/graph.txt", sep="\t", index=False, header=False)

    # Hierarchy file - CSV file with two columns: child, parent
    drug_hierarchy_df = get_drug_hierarchy(subgraph, output_dir)
    # disease_hierarchy_df = get_disease_hierarchy(subgraph)
    # hierarchy_df = pd.concat([drug_hierarchy_df, disease_hierarchy_df], ignore_index=True)

    # hierarchy_df["parent"] = hierarchy_df["parent"].apply(lambda x: x.lower() if x in ["Drug", "Disease"] else x)
    # hierarchy_df.to_csv("./dreamwalk_data/hierarchy.txt", sep="\t", index=False)

    d_nodes = [node for node, label in subgraph.nodes(data="labels") if label == "Disease"]
    t = []
    for node in d_nodes:
        t.append({"child": node, "parent": "disease"})

    hierarchy_df = pd.concat([drug_hierarchy_df, pd.DataFrame(t)], ignore_index=True)
    hierarchy_df["parent"] = hierarchy_df["parent"].apply(
        lambda x: x.lower() if x in ["Drug", "Disease"] else x
    )
    hierarchy_df.to_csv(f"{output_dir}/hierarchy.csv", index=False)

    # Positive/negative drug-disease association file - TSV file with three columns: drug, disease, label
    get_drug_disease_file(subgraph, output_dir)


def get_drug_hierarchy(graph: BioGraph, output_dir: str):
    """Generating the drug hierarchy using ATC classification."""
    drug_classes = []

    atc_hierarchy = pd.read_csv(
        f"{output_dir}/atc_hierarchy.csv", usecols=["dbID", "atcClassification", "id"]
    )

    atc_hierarchy = atc_hierarchy.rename(
        columns={"dbID": "drug", "atcClassification": "atc", "id": "chembl_id"}
    )
    atc_hierarchy["drug"] = "DrugBank:" + atc_hierarchy["drug"]
    atc_hierarchy["chembl_id"] = "ChEMBL:" + atc_hierarchy["chembl_id"]

    # With the main ATC hierarchy
    for atc_classes in atc_hierarchy["atc"].unique():
        for class_record in atc_classes.split(";"):
            cnames = class_record.split(",")
            cnames.reverse()

            parents = cnames[:-1]
            children = cnames[1:]

            for parent, child in zip(parents, children):
                drug_classes.append({"child": child, "parent": parent})

    # With hierarch of drugs in KG
    drug_nodes = [node for node, label in graph.nodes(data="labels") if label == "Compound"]

    skipped = 0
    for drug in drug_nodes:
        node_info = graph.nodes[drug]
        if "drugbank_id" in node_info:
            atc_classes = atc_hierarchy[atc_hierarchy["drug"] == node_info["drugbank_id"]][
                "atc"
            ].values
        elif "chembl_id" in node_info:
            atc_classes = atc_hierarchy[atc_hierarchy["chembl_id"] == node_info["chembl_id"]][
                "atc"
            ].values

        if len(atc_classes) == 0:
            skipped += 1
            continue

        if len(atc_classes) > 1:
            print(f"Multiple ATC classes for {drug}: {atc_classes}")

        atc_classes = atc_classes[0]
        cnames = atc_classes.split(",")
        parent = cnames[0]
        drug_classes.append({"child": drug, "parent": parent})

    print(f"Skipped {skipped} drugs out of {len(drug_nodes)}")
    return pd.DataFrame(drug_classes)


def generate_disease_hierarchy(disease_df, kg_data) -> pd.DataFrame:
    # filter rows with ':Drug' in '_labels'
    diseases_filtered = kg_data[kg_data["_labels"].isin([":Disease"])][["diseaseID", "class"]]
    # disease_hierarchy = generate_disease_hierarchy(diseases_filtered)

    # writing the merged hierarchy to a csv file
    # hierarchy_frames = [drug_hierarchy, disease_hierarchy]
    # hierarchy_result = pd.concat(hierarchy_frames)
    # hierarchy_result.to_csv('hierarchy.csv', sep=",", index = False)
    print("Hierarchy file is saved!")

    disease_hierarchy_df = pd.DataFrame(columns=["child", "parent"])
    disease_hierarchy_dict = {}
    for index, row in disease_df.iterrows():
        diseaseID = row["diseaseID"]
        mesh_classifications = row["class"].split(";")
        for mesh_classification in mesh_classifications:
            disease_hierarchy_df.loc[len(disease_hierarchy_df) + 1] = {
                "child": diseaseID,
                "parent": mesh_classification,
            }
            disease_hierarchy_dict[mesh_classification] = mesh_classification[0:1]
            disease_hierarchy_dict[mesh_classification[0:1]] = "disease"

    disease_hierarchy_df2 = pd.DataFrame(
        disease_hierarchy_dict.items(), columns=["child", "parent"]
    )
    disease_hierarchy_df = pd.concat([disease_hierarchy_df, disease_hierarchy_df2])
    return disease_hierarchy_df
