"""Codes taken from DreamWalk repository: https://github.com/eugenebang/DREAMwalk"""

import os
import random

import networkx as nx
import numpy as np


def read_graph(edgeList, weighted=True, directed=False, delimiter="\t") -> nx.MultiDiGraph:
    """Reads the input network in networkx.
    :param edgeList: Edge list file
    :param weighted: True if the network is weighted
    :param directed: True if the network is directed
    :param delimiter: Delimiter for the edge list file
    :return: networkx graph
    """
    if weighted:
        G = nx.read_edgelist(
            edgeList,
            nodetype=str,
            data=(("type", int), ("weight", float), ("id", int)),
            create_using=nx.MultiDiGraph(),
            delimiter=delimiter,
        )
    else:
        G = nx.read_edgelist(
            edgeList,
            nodetype=str,
            data=(("type", int)),
            create_using=nx.MultiDiGraph(),
            delimiter=delimiter,
        )
        for edge in G.edges():
            edge = G[edge[0]][edge[1]]
            for i in range(len(edge)):
                edge[i]["weight"] = 1.0

    if not directed:
        G = G.to_undirected()

    return G


def set_seed(seed=42):
    """Set seed for reproducibility.
    :param seed: Seed value
    """
    os.environ["PYTHONHASHSEED"] = str(seed)
    random.seed(seed)
    np.random.seed(seed)
