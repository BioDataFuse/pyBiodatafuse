# Code to get the disease hierarchy from the MONDO ontology

import pandas as pd
from owlready2 import *
from tqdm import tqdm

from pyBiodatafuse.analyzer.summarize import BioGraph


def get_direct_parent(label_idx):
    """Get the direct parent of the disease.
    :param label_idx: Disease label identifier
    """
    sqarl_query = """
        PREFIX oboinowl: <http://www.geneontology.org/formats/oboInOwl#>
        SELECT distinct ?mondo_id ?l
        { 
            ?x oboinowl:id ?dis_name .
            ?x rdfs:subClassOf ?y
            ?y rdfs:subClassOf ?z
            ?y oboinowl:id ?mondo_id
            ?y rdfs:label ?l
        }
        group by ?x
        order by count(?y)
        """
    sparql_query = sqarl_query.replace("?dis_name", f'"{label_idx}"')

    k = list(default_world.sparql(sparql_query))
    return k


def get_lineage(idx, lineage_list):
    """Get the lineage of the disease.
    :param idx: Disease label identifier
    :param disease_lineage: List to store the lineage of the disease
    """
    for parent in get_direct_parent(idx):
        lineage_list.append(parent[1])
        get_lineage(parent[0], lineage_list)

    return lineage_list


if __name__ == "__main__":
    graph_obj = BioGraph(graph_path="PCS_graph.gml", graph_format="gml")
    graph = graph_obj.get_subgraph(node_types=["Gene", "Disease", "Compound"])

    d_nodes2 = []

    for node, data_label in graph.nodes(data=True):
        if data_label["labels"] != "Disease":
            continue
        data_label["id"] = node
        d_nodes2.append(data_label)

    d_nodes2 = pd.DataFrame(d_nodes2)
    tmp = d_nodes2[["id", "name", "MONDO", "MESH", "DO"]]

    onto = get_ontology("http://purl.obolibrary.org/obo/mondo.owl").load()

    parent_list = []

    for umls, mondo in tqdm(tmp[["id", "MONDO"]].values):
        if pd.isna(mondo):
            parent_list.append({"umls": umls, "mondo": mondo, "parent": "disease"})
            continue

        multiple_mondo = mondo.split(",")
        parents = get_lineage(multiple_mondo[0], lineage_list=[])
        if len(parents) == 0:
            parent_list.append({"umls": umls, "mondo": mondo, "parent": "disease"})
        else:
            parent_list.append({"umls": umls, "mondo": mondo, "parent": ",".join(parents)})

    parent_df = pd.DataFrame(parent_list)
    parent_df.to_csv("dreamwalk_data/mondo_hierarchy.tsv", sep="\t", index=False)
