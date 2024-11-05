"""Codes taken from DreamWalk repository: https://github.com/eugenebang/DREAMwalk."""

import logging
import math
from collections import Counter, defaultdict
from itertools import combinations_with_replacement
from typing import Dict, List, Set

import pandas as pd

from pyBiodatafuse.algorithms.DREAMwalk.utils import read_graph

logger = logging.getLogger(__name__)


def _generate_tree(hier_df: pd.DataFrame, nodes: list) -> Dict[str, Dict[str, List]]:
    """Generate tree from hierarchy dataframe.

    :param hier_df: hierarchy dataframe
    :param nodes: list of nodes
    :returns: tree
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

    tree = defaultdict(lambda: defaultdict(list))  # type: Dict[str, Dict[str, List]]
    for row in rows:
        ntype = row[-1]
        n = row[0]
        tree[ntype][n].append(row)
    return tree


def _calculate_ic(total_counts: dict, max_wn: int, nodes: list) -> Dict[str, float]:
    """Calculate information content for each node.

    :param total_counts: total counts of nodes
    :param max_wn: maximum count of nodes
    :param nodes: list of nodes
    :returns: information content values
    """
    ic_values = defaultdict(float)  # type: Dict[str, float]
    for entity in total_counts:
        if entity in nodes:  # if node is in graph
            ic_value = 1.0
        else:
            count = total_counts[entity]
            ic_value = round(1 - math.log(count) / math.log(max_wn), 3)
        ic_values[entity] = ic_value
    return ic_values


def _ic_from_tree(tree: Dict[str, Dict[str, List]], nodes: list) -> Dict[str, Dict[str, float]]:
    """Calculate information content values from tree.

    :param tree: tree
    :param nodes: list of nodes
    :returns: information content values
    """
    ic_values = defaultdict(dict)  # type: Dict[str, Dict[str, float]]
    for ntype, n_vals in tree.items():
        total_counts = Counter()  # type: Counter
        for rows in n_vals.values():
            for row in rows:
                total_counts += Counter(row)
        max_wn = total_counts[ntype]
        ic_values[ntype] = _calculate_ic(total_counts, max_wn, nodes)
    return ic_values


def _sim_jc_from_tree(
    id1: str, id2: str, tree: Dict[str, List], ic_values: Dict[str, float]
) -> float:
    """Calculate similarity between two nodes using Jiang-Conrath similarity.

    :param id1: node 1
    :param id2: node 2
    :param tree: tree
    :param ic_values: information content values
    :returns: similarity value
    """
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
    com_ancs = set()  # type: Set[str]
    for pair in path_pairs:
        com_ancs = com_ancs | (set(pair[0]) & set(pair[1]))

    # get max_ic among all common ancestors
    max_ic = 0.0
    max_anc = ""
    for anc in com_ancs:
        anc_ic = ic_values[anc]
        if anc_ic > max_ic:
            max_ic = anc_ic
            max_anc = anc

    if max_anc == "":  # no common ancestor
        return 0

    # calculate similarity
    sim_jc = round(1 - (ic_values[id1] + ic_values[id2] - 2 * max_ic) / 2, 3)
    return sim_jc


def generate_sim_graph(
    hier_df: str, nodes: list, cutoff: float, directed: bool = True
) -> Dict[str, List]:
    """Generate similarity graph.

    :param hier_df: hierarchy dataframe
    :param nodes: list of nodes
    :param cutoff: cutoff value for similarity
    :param directed: directed or undirected graph
    :returns: similarity values
    """
    tree = _generate_tree(hier_df, nodes)
    ic_values = _ic_from_tree(tree, nodes)

    sim_values = defaultdict(list)
    for ntype, nvals in tree.items():
        ids = list(nvals.keys())
        all_pairs = list(combinations_with_replacement(ids, 2))

        for id1, id2 in all_pairs:
            sim = _sim_jc_from_tree(id1, id2, tree[ntype], ic_values[ntype])

            if sim < cutoff:
                continue
            sim_values[ntype].append((id1, id2, sim))

            if directed:
                sim_values[ntype].append((id2, id1, sim))
    return sim_values


def save_sim_graph(
    networkf: str,
    hierf: str,
    outputf: str,
    cutoff: float,
    weighted: bool = True,
    directed: bool = False,
    net_delimiter: str = "\t",
) -> None:
    """Generate similarity graph for the algorithm.

    :param networkf: Path to network file
    :param hierf: hierarchy file
    :param outputf: output file name
    :param cutoff: cutoff value for similarity
    :param weighted: weighted or unweighted graph
    :param directed: directed or undirected graph
    :param net_delimiter: delimiter of networks file; default = tab
    """
    graph = read_graph(networkf, weighted=weighted, directed=directed, delimiter=net_delimiter)
    nodes = list(graph.nodes())
    hier_df = pd.read_csv(hierf)

    logger.warning("Generating similarity graph...")
    sim_values = generate_sim_graph(hier_df, nodes, cutoff, directed)
    logger.warning("Similarity graph generated!")

    d = []
    for nvals in sim_values.values():
        index_num = 0
        type_number = 1
        for row in nvals:
            id1, id2, weight = row
            d.append([id1, id2, type_number, weight, index_num])

            index_num += 1
        type_number += 1

    df = pd.DataFrame(d, columns=["source", "target", "type", "weight", "edge_id"])
    df.to_csv(outputf, sep=net_delimiter, index=False)
    logger.warning(f"Similarity graph saved: {outputf}")
