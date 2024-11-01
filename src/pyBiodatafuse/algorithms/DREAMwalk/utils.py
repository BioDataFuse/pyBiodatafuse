"""Codes taken from DreamWalk repository: https://github.com/eugenebang/DREAMwalk."""

import os
import random

import networkx as nx
import numpy as np


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
    if weighted:
        graph = nx.read_edgelist(
            edgelist_path,
            nodetype=str,
            data=(("type", int), ("weight", float), ("id", int)),
            create_using=nx.MultiDiGraph(),
            delimiter=delimiter,
        )
    else:
        graph = nx.read_edgelist(
            edgelist_path,
            nodetype=str,
            data=(("type", int)),
            create_using=nx.MultiDiGraph(),
            delimiter=delimiter,
        )
        for edge in graph.edges():
            edge = graph[edge[0]][edge[1]]
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
