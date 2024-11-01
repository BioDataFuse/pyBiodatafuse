"""Codes taken from DreamWalk repository: https://github.com/eugenebang/DREAMwalk."""

import logging
import math
import os
import pickle
import random
from collections import Counter, defaultdict
from typing import Union

import networkx as nx
import numpy as np
import parmap
from scipy import stats

from pyBiodatafuse.algorithms.DREAMwalk.HeterogeneousSG import run_heterogeneous_sg
from pyBiodatafuse.algorithms.DREAMwalk.utils import read_graph, set_seed

logger = logging.getLogger(__name__)


def _init_edge_transition_matrix(graph: nx.MultiDiGraph) -> np.ndarray:
    """Initialize edge type transition matrix.

    :param graph: NetworkX graph object
    :returns: Edge type transition matrix
    """
    edgetypes = set()

    for _, _, d in graph.edges(data=True):
        edgetypes.add(d["type"])

    type_count = len(edgetypes)  # Number of edge types

    matrix = np.ones((type_count, type_count))
    matrix = matrix / matrix.sum()
    return matrix


def _edge_transition_walk(
    edge: tuple, graph: nx.MultiDiGraph, matrix: np.ndarray, walk_length: int, p: float, q: float
) -> list:
    """Generate edge transition walk.

    :param edge: Edge tuple
    :param graph: NetworkX graph object
    :param matrix: Edge type transition matrix
    :param walk_length: Length of random walk
    :param p: Return parameter
    :param q: In-out parameter
    :returns: Edge path
    """
    edge = (edge[0], edge[1], edge[2]["type"])
    walk = [edge]
    edge_path = [edge[2]]

    while len(walk) < walk_length:
        cur_edge = walk[-1]
        prev_node = cur_edge[0]
        cur_node = cur_edge[1]
        cur_edge_type = cur_edge[2]

        nbrs = graph.neighbors(cur_node)
        # calculate edge weights for all neighbors
        nbrs_list = []
        weights_list = []

        for nbr in nbrs:
            for nbr_edge in graph[cur_node][nbr].values():
                nbr_edge_type = nbr_edge["type"]
                nbr_edge_weight = nbr_edge["weight"]
                trans_weight = matrix[cur_edge_type - 1][nbr_edge_type - 1]

                if graph.has_edge(nbr, prev_node) or graph.has_edge(prev_node, nbr):
                    nbrs_list.append((nbr, nbr_edge_type))
                    weights_list.append(trans_weight * nbr_edge_weight)

                elif nbr == prev_node:  # p: Return parameter
                    nbrs_list.append((nbr, nbr_edge_type))
                    weights_list.append(trans_weight * nbr_edge_weight / p)

                else:  # q: In-out parameter q
                    nbrs_list.append((nbr, nbr_edge_type))
                    weights_list.append(trans_weight * nbr_edge_weight / q)

        if sum(weights_list) > 0:  # the direction_node has next_link_end_node
            next_edge = random.choices(nbrs_list, weights=weights_list)[0]
            next_edge_updated = (cur_node, next_edge[0], next_edge[1])
            walk.append(next_edge_updated)
            edge_path.append(next_edge_updated[2])
        else:
            break
    return edge_path


def _sample_edge_paths(
    graph: nx.MultiDiGraph, trans_matrix: np.ndarray, walk_length: int, p: float, q: float
) -> list:
    """Sample edge paths from the network.

    :param graph: NetworkX graph object
    :param trans_matrix: Edge type transition matrix
    :param walk_length: Length of random walk
    :param p: Return parameter
    :param q: In-out parameter
    :returns: Sampled edge paths
    """
    edges = list(graph.edges(data=True))

    # sample 1% of edges from the original network
    sampled_edges = random.sample(edges, int(len(edges) * 0.01))

    edge_walks = []

    for edge in sampled_edges:
        edge_walks.append(_edge_transition_walk(edge, graph, trans_matrix, walk_length, p, q))

    return edge_walks


def sigmoid(x: float) -> float:
    """Calculate sigmoid function.

    :param x: Input value
    :returns: Sigmoid value
    """
    return 1 / (1 + math.exp(-x))


def pearsonr_test(v1: list, v2: list) -> float:
    """Calculate Pearson correlation coefficient.

    :param v1: List of values
    :param v2: List of values
    :returns: Pearson correlation coefficient
    """
    result = stats.mstats.pearsonr(v1, v2)[0]
    return sigmoid(result)


def _update_trans_matrix(walks: list, matrix: np.ndarray) -> np.ndarray:
    """Update edge type transition matrix.

    :param walks: List of edge walks
    :param matrix: Edge type transition matrix
    :returns: Updated edge type transition matrix
    """
    type_count = len(matrix)
    matrix = np.zeros(matrix.shape)
    repo = defaultdict(list)

    for walk in walks:
        walk = [i - 1 for i in walk]
        edge_count = Counter(walk)
        #         if edge_id in curr
        for i in range(type_count):
            repo[i].append(edge_count[i])
    for i in range(type_count):
        for j in range(type_count):
            sim_score = pearsonr_test(repo[i], repo[j])
            matrix[i][j] = sim_score
    return np.nan_to_num(matrix)


def train_edgetype_transition_matrix(
    em_max_iter: int, graph: nx.MultiDiGraph, walk_length: int, p: float, q: float
) -> np.ndarray:
    """Train edge type transition matrix using EM algorithm.

    :param em_max_iter: Maximum number of EM iterations
    :param graph: NetworkX graph object
    :param walk_length: Length of random walk
    :param p: Return parameter
    :param q: In-out parameter
    :returns: Edge type transition matrix
    """
    matrix_conv_rate = 0.01
    matrices = {0: _init_edge_transition_matrix(graph=graph)}

    for i in range(em_max_iter):  # EM iteration
        walks = _sample_edge_paths(graph, matrices[i], walk_length, p, q)  # M step
        matrices[i + 1] = _update_trans_matrix(walks, matrices[i])  # E step
        matrix_diff = np.nan_to_num(
            np.absolute((matrices[i + 1] - matrices[i]) / matrices[i])
        ).mean()
        if matrix_diff < matrix_conv_rate:
            break
    return matrices[i + 1]


def _network_traverse(
    cur: str, prev: tuple, graph: nx.MultiDiGraph, trans_matrix: np.ndarray, p: float, q: float
) -> tuple:
    """Traverse network to find next node and edge type.

    :param cur: Current node
    :param prev: Previous node and edge type
    :param graph: NetworkX graph object
    :param trans_matrix: Edge type transition matrix
    :param p: Return parameter
    :param q: In-out parameter
    :returns: Next node and edge type
    """
    prev_node, cur_edge_type = prev
    cur_nbrs = sorted(graph.neighbors(cur))
    nbrs_list = []
    weights_list = []

    # search for reachable edges and their weights
    for nbr in cur_nbrs:
        nbr_edges = graph[cur][nbr]
        for nbr_edge in nbr_edges.values():
            nbr_edge_type = nbr_edge["type"]
            nbr_edge_weight = nbr_edge["weight"]
            trans_weight = trans_matrix[cur_edge_type - 1][nbr_edge_type - 1]

            if graph.has_edge(nbr, prev_node) or graph.has_edge(prev_node, nbr):
                nbrs_list.append((nbr, nbr_edge_type))
                weights_list.append(trans_weight * nbr_edge_weight)

            elif nbr == prev_node:  # p: Return parameter
                nbrs_list.append((nbr, nbr_edge_type))
                weights_list.append(trans_weight * nbr_edge_weight / p)

            else:  # q: In-Out parameter
                nbrs_list.append((nbr, nbr_edge_type))
                weights_list.append(trans_weight * nbr_edge_weight / q)

    # sample next node and edge type from searched weights
    next_edge = random.choices(nbrs_list, weights=weights_list)[0]
    return next_edge


def _teleport_operation(cur: str, g_sim: nx.MultiDiGraph) -> str:
    """Perform teleport operation.

    :param cur: Current node
    :param g_sim: NetworkX graph object
    :returns: Next node
    """
    cur_nbrs = sorted(g_sim.neighbors(cur))
    random.shuffle(cur_nbrs)
    selected_nbrs = []
    distance_sum = 0
    for nbr in cur_nbrs:
        nbr_links = g_sim[cur][nbr]
        for i in nbr_links:
            nbr_link_weight = nbr_links[i]["weight"]
            distance_sum += nbr_link_weight
            selected_nbrs.append(nbr)

    rand = np.random.rand() * distance_sum
    threshold = 0

    for nbr in set(selected_nbrs):
        nbr_links = g_sim[cur][nbr]
        for i in nbr_links:
            nbr_link_weight = nbr_links[i]["weight"]
            threshold += nbr_link_weight
            if threshold >= rand:
                next = nbr
                break
    return next


def _dreamwalker(
    start_node: str,
    graph: nx.MultiDiGraph,
    g_sim: nx.MultiDiGraph,
    trans_matrix: np.ndarray,
    p: float,
    q: float,
    walk_length: int,
    tp_factor: float,
) -> list:
    """Generate DREAMwalk path.

    :param start_node: Start node
    :param graph: NetworkX graph object
    :param g_sim: NetworkX graph object
    :param trans_matrix: Edge type transition matrix
    :param p: Return parameter
    :param q: In-out parameter
    :param walk_length: Length of random walk
    :param tp_factor: Teleport factor
    :returns: DREAMwalk path
    """
    walk = [start_node]
    edge_walk = []

    # select first edge from any neighbors
    nbrs_list = []
    weights_list = []
    for nbr in sorted(graph.neighbors(start_node)):
        for cur_edge in graph[start_node][nbr].values():
            nbrs_list.append((nbr, cur_edge["type"]))
            weights_list.append(cur_edge["weight"])
    next_edge = random.choices(nbrs_list, weights=weights_list)[0]
    walk.append(next_edge[0])
    edge_walk.append(next_edge[1])

    while len(walk) < walk_length:
        prev = walk[-2]
        cur = walk[-1]
        cur_edge_type = edge_walk[-1]

        cur_nbrs = sorted(graph.neighbors(cur))
        if len(cur_nbrs) > 0:
            # perform DREAMwalk path generation
            if (cur in g_sim.nodes()) & (np.random.rand() < tp_factor):
                next_node = _teleport_operation(cur, g_sim)
            else:
                prev_data = (prev, cur_edge_type)
                next_node, next_edge_type = _network_traverse(
                    cur, prev_data, graph, trans_matrix, p, q
                )
                edge_walk.append(next_edge_type)
            walk.append(next_node)
        else:  # dead end
            break  # if start node has 0 neighbour : dead end
    return walk


# parallel walks
def _parmap_walks(
    _: list,
    nodes: list,
    graph: nx.MultiDiGraph,
    g_sim: nx.MultiDiGraph,
    trans_matrix: np.ndarray,
    p: float,
    q: float,
    walk_length: int,
    tp_factor: float,
) -> list:
    """Generate DREAMwalk paths in parallel.

    :param _: Dummy variable
    :param nodes: List of nodes
    :param graph: NetworkX graph object
    :param g_sim: NetworkX graph object
    :param trans_matrix: Edge type transition matrix
    :param p: Return parameter
    :param q: In-out parameter
    :param walk_length: Length of random walk
    :param tp_factor: Teleport factor
    :returns: DREAMwalk paths
    """
    walks = []
    random.shuffle(nodes)
    for node in nodes:
        walks.append(_dreamwalker(node, graph, g_sim, trans_matrix, p, q, walk_length, tp_factor))
    return walks


# Teleport guided random walk
def generate_dreamwalk_paths(
    graph: nx.MultiDiGraph,
    g_sim: nx.MultiDiGraph,
    trans_matrix: np.ndarray,
    p: float,
    q: float,
    num_walks: int,
    walk_length: int,
    tp_factor: float,
    workers: int,
) -> list:
    """Generate DREAMwalk paths.

    :param graph: NetworkX graph object
    :param g_sim: NetworkX graph object
    :param trans_matrix: Edge type transition matrix
    :param p: Return parameter
    :param q: In-out parameter
    :param num_walks: Number of walks
    :param walk_length: Length of random walk
    :param tp_factor: Teleport factor
    :param workers: Number of workers
    :returns: DREAMwalk paths
    """
    tot_walks = []
    nodes = list(graph.nodes())

    walks = parmap.map(
        _parmap_walks,
        range(num_walks),
        nodes,
        graph,
        g_sim,
        trans_matrix,
        p,
        q,
        walk_length,
        tp_factor,
        pm_pbar=False,
        pm_processes=min(num_walks, workers),
    )
    for walk in walks:
        tot_walks += walk

    return tot_walks


def save_embedding_files(
    netf: str,
    sim_netf: str,
    outputf: str,
    nodetypef: str,
    tp_factor: float = 0.5,
    directed: bool = False,
    weighted: bool = True,
    em_max_iter: int = 5,
    num_walks: int = 100,
    walk_length: int = 10,
    dimension: int = 128,
    window_size: int = 4,
    p: float = 1,
    q: float = 1,
    net_delimiter: str = "\t",
) -> None:
    """Generate node embeddings using DREAMwalk algorithm.

    :param netf: Network file path
    :param sim_netf: Similarity network file path
    :param outputf: Output file path
    :param nodetypef: Node type file path
    :param tp_factor: Teleport factor
    :param directed: Boolean for whether network is directed or not
    :param weighted: Boolean for whether network is weighted or not
    :param em_max_iter: Maximum number of EM iterations
    :param num_walks: Number of walks
    :param walk_length: Length of random walk
    :param dimension: Embedding dimension
    :param window_size: Window size
    :param p: Return parameter
    :param q: In-out parameter
    :param net_delimiter: Delimiter for network file
    :raises ValueError: If input similarity file is not detected
    """
    set_seed(42)
    cpus = os.cpu_count()
    if cpus is None:
        workers = 1
    else:
        workers = cpus - 2
    graph = read_graph(netf, weighted=weighted, directed=directed, delimiter=net_delimiter)
    if not os.path.exists(sim_netf):
        raise ValueError(
            f"No input similarity file detected! Please check the file path - {sim_netf} or pass the input explicitly with parameter `sim_netf`."
        )
    g_sim = read_graph(sim_netf, weighted=True, directed=False)

    logging.info("Training edge type transition matrix...")
    trans_matrix = train_edgetype_transition_matrix(em_max_iter, graph, walk_length, p, q)

    logging.info("Generating paths...")
    walks = generate_dreamwalk_paths(
        graph=graph,
        g_sim=g_sim,
        trans_matrix=trans_matrix,
        p=p,
        q=q,
        num_walks=num_walks,
        walk_length=walk_length,
        tp_factor=tp_factor,
        workers=workers,
    )

    logging.info("Generating node embeddings...")
    use_hetsg = True if nodetypef is not None else False
    embeddings = run_heterogeneous_sg(
        use_hetsg,
        walks,
        set(graph.nodes()),
        nodetypef=nodetypef,
        embedding_size=dimension,
        window_length=window_size,
        workers=workers,
    )
    with open(outputf, "wb") as fw:
        pickle.dump(embeddings, fw)

    logging.info(f"Node embeddings saved: {outputf}")
