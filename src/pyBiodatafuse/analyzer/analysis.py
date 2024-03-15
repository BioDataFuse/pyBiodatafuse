import networkx as nx
import pandas as pd
from generator import generate_networkx_graph


def main():
    df_test = pd.read_pickle("./combined_df.pkl")
    print(df_test.head())
    g = generate_networkx_graph(df_test)
    print("Number of nodes in graph: {}".format(len(g.nodes)))
    print(len(g.edges))

    max_out_degree = 0  # node associated to most diseases
    node_max = None
    for node in g.nodes:
        if g.out_degree(node) > max_out_degree:
            max_out_degree = g.out_degree(node)
            node_max = node

    print(
        "Node with most diseases associated: {} with {} disease associations known".format(
            node_max, max_out_degree
        )
    )

    # print(list(g.edges(data=True))[:10])
    # print(g.get_edge_data("VAMP1", "CHRNE"))
    # print(list(g.edges(data="label", default=1, keys=True))[0])

    # fetch back all labels in the graph, this can help filter by edge type later
    labels = set()
    for u, v, k in g.edges(data=True):
        labels.add(k["label"])

    print("Labels: {}".format(labels))

    # e.g. retrieving back gene - disease association links only:
    gene_disease = ((u, v) for u, v, d in g.edges(data=True) if d["label"] == "associated_with")

    # how many edges are of type associated_with?
    print(len(list(gene_disease)))

    # compute an overview of edge types
    for label_type in labels:
        subgraph = ((u, v) for u, v, d in g.edges(data=True) if d["label"] == label_type)
        print("For label type {} there are {} edges.".format(label_type, len(list(subgraph))))

    # extract interaction subgraph
    ppi_edges = ((u, v) for u, v, d in g.edges(data=True) if d["label"] == "interacts_with")

    ppi_subgraph = nx.DiGraph()
    nodes = set()
    for u, v in ppi_edges:
        ppi_subgraph.add_node(u)
        ppi_subgraph.add_node(v)
        nodes.add(u)
        nodes.add(v)
        ppi_subgraph.add_edge(u, v)

    communities_generator = nx.community.girvan_newman(ppi_subgraph)
    top_level_communities = next(communities_generator)
    print("Top level communities: {}".format(sorted(map(sorted, top_level_communities))))

    gene_disease_edges = (
        (u, v) for u, v, d in g.edges(data=True) if d["label"] == "associated_with"
    )

    gene_disease_subgraph = nx.DiGraph()
    nodes = set()
    for u, v in gene_disease_edges:
        gene_disease_subgraph.add_node(u)
        gene_disease_subgraph.add_node(v)
        nodes.add(u)
        nodes.add(v)
        gene_disease_subgraph.add_edge(u, v)

    print(len(nodes))
    communities_generator = nx.community.girvan_newman(gene_disease_subgraph.to_undirected())

    top_level_communities = next(communities_generator)
    print("Top level communities: {}".format(sorted(map(sorted, top_level_communities))))

    # basic link prediction using Jaccard
    predictions = list(nx.jaccard_coefficient(gene_disease_subgraph.to_undirected()))
    non_zero_predictions = [
        (gene, disease, predicted_score)
        for (gene, disease, predicted_score) in predictions
        if predicted_score > 0 and predicted_score != 1
    ]

    # sort by prediction score descending
    non_zero_predictions.sort(key=lambda x: x[2], reverse=True)

    # print top 10 predictions
    # ... actually predicts disease - disease links !
    print("Predicted links: {}".format(non_zero_predictions[:10]))


if __name__ == "__main__":
    main()
