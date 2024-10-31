"""Codes taken from DreamWalk repository: https://github.com/eugenebang/DREAMwalk"""

import logging
import math
import os
import pickle
import random
from collections import Counter, defaultdict

import networkx as nx
import numpy as np
import parmap
from scipy import stats

from pyBiodatafuse.algorithms.DREAMwalk.HeterogeneousSG import HeterogeneousSG
from pyBiodatafuse.algorithms.DREAMwalk.utils import read_graph, set_seed

logger = logging.getLogger(__name__)


def train_edgetype_transition_matrix(em_max_iter, G, networkf, net_delimiter, walk_length, p, q):
    """Train edge type transition matrix using EM algorithm.
    :param em_max_iter: Maximum number of EM iterations
    :param G: NetworkX graph object
    :param networkf: Network file path
    :param net_delimiter: Delimiter for network file
    :param walk_length: Length of random walk
    :param p: Return parameter
    :param q: In-out parameter
    :return: Edge type transition matrix
    """
    matrix_conv_rate = 0.01
    matrices = {0: _init_edge_transition_matrix(networkf, net_delimiter)}

    for i in range(em_max_iter):  # EM iteration
        walks = _sample_edge_paths(G, matrices[i], walk_length, p, q)  # M step
        matrices[i + 1] = _update_trans_matrix(walks, matrices[i])  # E step
        matrix_diff = np.nan_to_num(
            np.absolute((matrices[i + 1] - matrices[i]) / matrices[i])
        ).mean()
        if matrix_diff < matrix_conv_rate:
            break
    return matrices[i + 1]


def _sample_edge_paths(G, trans_matrix, walk_length, p, q):
    """Sample edge paths from the network.
    :param G: NetworkX graph object
    :param trans_matrix: Edge type transition matrix
    :param walk_length: Length of random walk
    :param p: Return parameter
    :param q: In-out parameter
    """
    edges = list(G.edges(data=True))
    sampled_edges = random.sample(
        edges, int(len(edges) * 0.01)
    )  # sample 1% of edges from the original network
    edge_walks = []
    for edge in sampled_edges:
        edge_walks.append(_edge_transition_walk(edge, G, trans_matrix, walk_length, p, q))
    return edge_walks


def _edge_transition_walk(edge, G, matrix, walk_length, p, q):
    """Generate edge transition walk.
    :param edge: Edge tuple
    :param G: NetworkX graph object
    :param matrix: Edge type transition matrix
    :param walk_length: Length of random walk
    :param p: Return parameter
    """
    edge = (edge[0], edge[1], edge[2]["type"])
    walk = [edge]
    edge_path = [edge[2]]

    while len(walk) < walk_length:
        cur_edge = walk[-1]
        prev_node = cur_edge[0]
        cur_node = cur_edge[1]
        cur_edge_type = cur_edge[2]

        nbrs = G.neighbors(cur_node)
        # calculate edge weights for all neighbors
        nbrs_list = []
        weights_list = []
        for nbr in nbrs:
            for nbr_edge in G[cur_node][nbr].values():
                nbr_edge_type = nbr_edge["type"]
                nbr_edge_weight = nbr_edge["weight"]
                trans_weight = matrix[cur_edge_type - 1][nbr_edge_type - 1]
                if G.has_edge(nbr, prev_node) or G.has_edge(prev_node, nbr):
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
            next_edge = (cur_node, next_edge[0], next_edge[1])
            walk.append(next_edge)
            edge_path.append(next_edge[2])
        else:
            break
    return edge_path


def _init_edge_transition_matrix(networkf, net_delimiter):
    """Initialize edge type transition matrix.
    :param networkf: Network file path
    :param net_delimiter: Delimiter for network file
    """
    edgetypes = set()
    with open(networkf, "r") as fin:
        lines = fin.readlines()
    for line in lines:
        edgetypes.add(line.split(net_delimiter)[2])
    type_count = len(edgetypes)  # Number of edge types

    matrix = np.ones((type_count, type_count))
    matrix = matrix / matrix.sum()
    return matrix


def _update_trans_matrix(walks, matrix):
    """Update edge type transition matrix.
    :param walks: List of edge walks
    :param matrix: Edge type transition matrix
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


def pearsonr_test(v1, v2):  # original metric: the larger the more similar
    """Calculate Pearson correlation coefficient.
    :param v1: List of values
    :param v2: List of values
    """
    result = stats.mstats.pearsonr(v1, v2)[0]
    return sigmoid(result)


def sigmoid(x):
    """Calculate sigmoid function.
    :param x: Input value
    """
    return 1 / (1 + math.exp(-x))


# Teleport guided random walk
def generate_DREAMwalk_paths(
    G, G_sim, trans_matrix, p, q, num_walks, walk_length, tp_factor, workers
):
    """Generate DREAMwalk paths.
    :param G: NetworkX graph object
    :param G_sim: NetworkX graph object
    :param trans_matrix: Edge type transition matrix
    :param p: Return parameter
    :param q: In-out parameter
    :param num_walks: Number of walks
    :param walk_length: Length of random walk
    :param tp_factor: Teleport factor
    :param workers: Number of workers
    """
    tot_walks = []
    nodes = list(G.nodes())

    #     print(f'# of walkers : {min(num_walks,workers)} for {num_walks} times')
    walks = parmap.map(
        _parmap_walks,
        range(num_walks),
        nodes,
        G,
        G_sim,
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


# parallel walks
def _parmap_walks(_, nodes, G, G_sim, trans_matrix, p, q, walk_length, tp_factor):
    """Generate DREAMwalk paths in parallel.
    :param _: Dummy variable
    :param nodes: List of nodes
    :param G: NetworkX graph object
    :param G_sim: NetworkX graph object
    :param trans_matrix: Edge type transition matrix
    :param p: Return parameter
    :param q: In-out parameter
    :param walk_length: Length of random walk
    :param tp_factor: Teleport factor
    """
    walks = []
    random.shuffle(nodes)
    for node in nodes:
        walks.append(_DREAMwalker(node, G, G_sim, trans_matrix, p, q, walk_length, tp_factor))
    return walks


def _DREAMwalker(start_node, G, G_sim, trans_matrix, p, q, walk_length, tp_factor):
    """Generate DREAMwalk path.
    :param start_node: Start node
    :param G: NetworkX graph object
    :param G_sim: NetworkX graph object
    :param trans_matrix: Edge type transition matrix
    :param p: Return parameter
    :param q: In-out parameter
    :param walk_length: Length of random walk
    :param tp_factor: Teleport factor
    """
    walk = [start_node]
    edge_walk = []

    # select first edge from any neighbors
    nbrs_list = []
    weights_list = []
    for nbr in sorted(G.neighbors(start_node)):
        for cur_edge in G[start_node][nbr].values():
            nbrs_list.append((nbr, cur_edge["type"]))
            weights_list.append(cur_edge["weight"])
    next_edge = random.choices(nbrs_list, weights=weights_list)[0]
    walk.append(next_edge[0])
    edge_walk.append(next_edge[1])

    while len(walk) < walk_length:
        prev = walk[-2]
        cur = walk[-1]
        cur_edge_type = edge_walk[-1]

        cur_nbrs = sorted(G.neighbors(cur))
        if len(cur_nbrs) > 0:
            # perform DREAMwalk path generation
            if (cur in G_sim.nodes()) & (np.random.rand() < tp_factor):
                next_node = _teleport_operation(cur, G_sim)
            else:
                prev = (prev, cur_edge_type)
                next_node, next_edge_type = _network_traverse(cur, prev, G, trans_matrix, p, q)
                edge_walk.append(next_edge_type)
            walk.append(next_node)
        else:  # dead end
            break  # if start node has 0 neighbour : dead end
    return walk


def _network_traverse(cur, prev, G, trans_matrix, p, q):
    """Traverse network to find next node and edge type.
    :param cur: Current node
    :param prev: Previous node and edge type
    :param G: NetworkX graph object
    :param trans_matrix: Edge type transition matrix
    :param p: Return parameter
    :param q: In-out parameter
    """
    prev_node = prev[0]
    cur_edge_type = prev[1]
    cur_nbrs = sorted(G.neighbors(cur))
    nbrs_list = []
    weights_list = []

    # search for reachable edges and their weights
    for nbr in cur_nbrs:
        nbr_edges = G[cur][nbr]
        for nbr_edge in nbr_edges.values():
            nbr_edge_type = nbr_edge["type"]
            nbr_edge_weight = nbr_edge["weight"]
            trans_weight = trans_matrix[cur_edge_type - 1][nbr_edge_type - 1]

            if G.has_edge(nbr, prev_node) or G.has_edge(prev_node, nbr):
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


def _teleport_operation(cur, G_sim):
    """Perform teleport operation.
    :param cur: Current node
    :param G_sim: NetworkX graph object
    """
    cur_nbrs = sorted(G_sim.neighbors(cur))
    random.shuffle(cur_nbrs)
    selected_nbrs = []
    distance_sum = 0
    for nbr in cur_nbrs:
        nbr_links = G_sim[cur][nbr]
        for i in nbr_links:
            nbr_link_weight = nbr_links[i]["weight"]
            distance_sum += nbr_link_weight
            selected_nbrs.append(nbr)
    if distance_sum == 0:
        return False
    rand = np.random.rand() * distance_sum
    threshold = 0

    for nbr in set(selected_nbrs):
        nbr_links = G_sim[cur][nbr]
        for i in nbr_links:
            nbr_link_weight = nbr_links[i]["weight"]
            threshold += nbr_link_weight
            if threshold >= rand:
                next = nbr
                break
    return next


def save_embedding_files(
    netf: str,
    sim_netf: str,
    outputf: str,
    nodetypef: str,
    tp_factor: float = 0.5,
    seed: int = 42,
    directed: bool = False,
    weighted: bool = True,
    em_max_iter: int = 5,
    num_walks: int = 100,
    walk_length: int = 10,
    workers: int = os.cpu_count() - 2,
    dimension: int = 128,
    window_size: int = 4,
    p: float = 1,
    q: float = 1,
    net_delimiter: str = "\t",
):

    set_seed(seed)
    G = read_graph(netf, weighted=weighted, directed=directed, delimiter=net_delimiter)
    if sim_netf:
        G_sim = read_graph(sim_netf, weighted=True, directed=False)
    else:
        G_sim = nx.empty_graph()
        raise ValueError("No input similarity file detected!")

    logging.info("Training edge type transition matrix...")
    trans_matrix = train_edgetype_transition_matrix(
        em_max_iter, G, netf, net_delimiter, walk_length, p, q
    )

    logging.info("Generating paths...")
    walks = generate_DREAMwalk_paths(
        G, G_sim, trans_matrix, p, q, num_walks, walk_length, tp_factor, workers
    )
    #     with open('tmp_walk_file.pkl','wb') as fw:
    #         pickle.dump(walks,fw)

    logging.info("Generating node embeddings...")
    use_hetSG = True if nodetypef != None else False
    embeddings = HeterogeneousSG(
        use_hetSG,
        walks,
        set(G.nodes()),
        nodetypef=nodetypef,
        embedding_size=dimension,
        window_length=window_size,
        workers=workers,
    )
    with open(outputf, "wb") as fw:
        pickle.dump(embeddings, fw)

    logging.info(f"Node embeddings saved: {outputf}")
