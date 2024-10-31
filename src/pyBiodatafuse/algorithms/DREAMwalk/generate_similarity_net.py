"""Codes taken from DreamWalk repository: https://github.com/eugenebang/DREAMwalk"""

import logging
import math
from collections import Counter, defaultdict

import pandas as pd

from pyBiodatafuse.algorithms.DREAMwalk.utils import read_graph

logger = logging.getLogger(__name__)


def _generate_tree(hier_df: pd.DataFrame, nodes: list) -> dict:
    """Generate tree from hierarchy dataframe.
    :param hier_df: hierarchy dataframe
    :param nodes: list of nodes
    :return: tree
    """
    i = 1
    tempdf = hier_df[hier_df["child"].isin(nodes)]
    tempdf.columns = [0, 1]

    while True:
        tempdf = pd.merge(tempdf, hier_df, left_on=i, right_on="child", how="left").drop(
            columns=["child"]
        )
        tempdf.columns = list(range(len(tempdf.columns)))
        i += 1
        if sum(~tempdf[i].isna()) == 0:
            break
    rows = [tempdf.iloc[i].dropna().tolist() for i in range(len(tempdf))]

    tree = defaultdict(dict)  # type: dict
    for row in rows:
        ntype = row[-1]
        try:
            tree[ntype][row[0]].append(row)
        except:
            tree[ntype][row[0]] = []
            tree[ntype][row[0]].append(row)
    return tree


def generate_sim_graph(hier_df: str, nodes: list, cutoff: float, directed: bool = True):
    """Generate similarity graph.
    :param hier_df: hierarchy dataframe
    :param nodes: list of nodes
    :param cutoff: cutoff value for similarity
    :param directed: directed or undirected graph
    :return: similarity values
    """
    tree = _generate_tree(hier_df, nodes)
    ic_values = _ic_from_tree(tree, nodes)

    sim_values = {}  # type: dict
    for ntype in tree.keys():
        sim_values[ntype] = []
        ids = list(tree[ntype].keys())
        for i in range(len(ids)):
            id1 = ids[i]
            for j in range(i + 1, len(ids)):
                id2 = ids[j]
                sim = _simJC_from_tree(id1, id2, tree[ntype], ic_values[ntype])
                if sim > cutoff:
                    sim_values[ntype].append((id1, id2, sim))
                    if directed:
                        sim_values[ntype].append((id2, id1, sim))
    return sim_values


def _calculate_ic(total_counts: dict, max_wn: int, nodes: list) -> dict:
    """Calculate information content for each node.
    :param total_counts: total counts of nodes
    :param max_wn: maximum count of nodes
    :param nodes: list of nodes
    """
    ic_values = {}
    for entity in list(total_counts.keys()):
        if entity in nodes:  # if node is in graph
            ic_value = 1
        else:
            count = total_counts[entity]
            ic_value = 1 - math.log(count) / math.log(max_wn)  # type: float
        ic_values[entity] = ic_value
    return ic_values


def _ic_from_tree(tree: dict, nodes: list):
    ic_values = {}
    for ntype in tree.keys():
        total_counts = Counter()  # type: Counter
        for rows in tree[ntype].values():
            for row in rows:
                total_counts += Counter(row)
        max_wn = total_counts[ntype]
        ic_values[ntype] = _calculate_ic(total_counts, max_wn, nodes)
    return ic_values


def _simJC_from_tree(id1: str, id2: str, tree, ic_values):
    """Calculate similarity between two nodes using Jiang-Conrath similarity.
    :param id1: node 1
    :param id2: node 2
    :param tree: tree
    :param ic_values: information content values
    :return: similarity value"""
    if id1 == id2:
        return 1.0

    try:
        paths1 = tree[id1]
        paths2 = tree[id2]
    except KeyError:  # if leaf does not have a path
        return 0

    path_pairs = []
    for path1 in paths1:
        for path2 in paths2:
            path_pairs.append([path1, path2])

    # get all common ancecstors
    com_ancs = set()  # type: set
    for pair in path_pairs:
        com_ancs = com_ancs | (set(pair[0]) & set(pair[1]))

    # get max_ic among all common ancestors
    max_ic = 0
    max_anc = False
    for anc in com_ancs:
        anc_ic = ic_values[anc]
        if anc_ic > max_ic:
            max_ic = anc_ic
            max_anc = anc
    if not max_anc:  # no common ancestor
        return 0

    # calculate similarity

    simJC = 1 - (ic_values[id1] + ic_values[id2] - 2 * max_ic) / 2
    return simJC


def save_sim_graph(
    networkf: str,
    hierf: str,
    outputf: str,
    cutoff: float,
    directed: bool = False,
    net_delimiter: str = "\t",
):
    """Generate similarity graph for the algorithm.
    :param networkf: Path to network file
    :param hierf: hierarchy file
    :param outputf: output file name
    :param cutoff: cutoff value for similarity
    :param directed: directed or undirected graph
    :param net_delimiter: delimiter of networks file; default = tab
    """
    G = read_graph(networkf, weighted=True, directed=directed, delimiter=net_delimiter)
    nodes = list(G.nodes())
    hier_df = pd.read_csv(hierf)

    logger.warning("Generating similarity graph...")
    sim_values = generate_sim_graph(hier_df, nodes, cutoff, directed)
    logger.warning("Similarity graph generated!")

    with open(outputf, "w") as fw:
        index = 0
        type_number = 1
        for ntype in sim_values.keys():
            for row in sim_values[ntype]:
                id1, id2, weight = row
                fw.write(
                    "{}{delim}{}{delim}{}{delim}{}{delim}{}\n".format(
                        id1, id2, type_number, weight, index, delim=net_delimiter
                    )
                )
                index += 1
            type_number += 1
    logger.warning(f"Similarity graph saved: {outputf}")
