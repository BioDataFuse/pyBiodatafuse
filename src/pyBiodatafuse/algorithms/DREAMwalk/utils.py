"""Codes taken from DreamWalk repository: https://github.com/eugenebang/DREAMwalk."""

import os
import random

import networkx as nx
import numpy as np
import pandas as pd


def read_graph(
    edgelist_path: str, weighted: bool = True, directed: bool = False, delimiter: str = "\t"
) -> nx.MultiDiGraph:
    """Read the input network in networkx.

    :param edgelist_path: Edge list file path
    :param weighted: True if the network is weighted
    :param directed: True if the network is directed
    :param delimiter: Delimiter for the edge list file
    :returns: networkx graph
    """
    df = pd.read_csv(edgelist_path, sep=delimiter)
    if weighted:
        graph = nx.from_pandas_edgelist(
            df,
            source="source",
            target="target",
            edge_attr=True,
            create_using=nx.MultiDiGraph(),
        )
    else:
        graph = nx.from_pandas_edgelist(
            df,
            source="source",
            target="target",
            edge_attr=True,
            create_using=nx.MultiDiGraph(),
        )
        for source, target in graph.edges():
            edge = graph[source][target]
            for i in range(len(edge)):
                edge[i]["weight"] = 1.0

    if not directed:
        graph = graph.to_undirected()

    return graph


def set_seed(seed: int = 42) -> None:
    """Set seed for reproducibility.

    :param seed: Seed value
    """
    os.environ["PYTHONHASHSEED"] = str(seed)
    random.seed(seed)
    np.random.seed(seed)
