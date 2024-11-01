"""Codes taken from DreamWalk repository: https://github.com/eugenebang/DREAMwalk."""

import os
from typing import Dict

import numpy as np
import pandas as pd

from pyBiodatafuse.algorithms.DREAMwalk.constant import CURRENT_DIR


def _prep_hetsg_walks(walks: list, node2id: Dict[str, str], nodetypef: str) -> list:
    """Prepare walks for HetSG algorithm.

    :param walks: List of walks
    :param node2id: Dictionary of node to id mapping
    :param nodetypef: Node type file path
    :returns: Annotated walks
    """
    node2type = pd.read_csv(nodetypef, sep="\t").set_index("node")["type"].to_dict()

    # we need to annotate the nodes with prefixes for the HetSG to recognize the nodetype
    type2meta = {"drug": "d", "disease": "i", "gene": "g"}

    annot_walks = []
    for walk in walks:
        annot_walk = []
        for node in walk:
            nodeid = node2id[node]
            try:
                node = type2meta[node2type[node]] + nodeid
            except KeyError:
                node = "e" + nodeid
            annot_walk.append(node)
        annot_walks.append(" ".join(annot_walk))
    return annot_walks


def _prep_sg_walks(walks: list, node2id: Dict[str, str]) -> list:
    """Prepare walks for SG algorithm.

    :param walks: List of walks
    :param node2id: Dictionary of node to id mapping
    :returns: Annotated walks
    """
    # does not require node type prefixes
    annot_walks = []
    for walk in walks:
        annot_walk = [node2id[node] for node in walk]
        annot_walks.append(" ".join(annot_walk))
    return annot_walks


def run_heterogeneous_sg(
    use_hetsg: bool,
    walks: list,
    nodes: set,
    nodetypef: str,
    embedding_size: int,
    window_length: int,
    workers: int,
) -> Dict[str, np.ndarray]:
    """Generate embeddings using HetSG algorithm.

    :param use_hetsg: Use HetSG or SG
    :param walks: List of walks
    :param nodes: Set of nodes
    :param nodetypef: Node type file path
    :param embedding_size: Embedding size
    :param window_length: Window length
    :param workers: Number of workers
    :returns: Embeddings
    """
    node2id = {node: str(i) for i, node in enumerate(nodes)}
    id2node = {idx: node for node, idx in node2id.items()}

    if use_hetsg:
        annot_walks = _prep_hetsg_walks(walks, node2id, nodetypef)
    else:
        annot_walks = _prep_sg_walks(walks, node2id)

    tmp_walkf = "tmp_walkfile"
    with open(tmp_walkf, "w") as fw:
        fw.write("\n".join(annot_walks))

    # use hetsg.cpp file for embedding vector generation -> outputs .txt file
    outputf = "tmp_outputf"
    cpp_command = f"g++ {CURRENT_DIR}/HeterogeneousSG.cpp -o {CURRENT_DIR}/HetSG -lm -pthread -O3 -march=native -Wall -funroll-loops -Wno-unused-result"
    os.system(cpp_command)

    run_command = f"{CURRENT_DIR}/HetSG -train {tmp_walkf} -output {outputf} -pp {int(use_hetsg)} -min-count 0 -size {embedding_size} -iter 1 -window {window_length} -threads {workers}"
    os.system(run_command)

    with open(outputf + ".txt", "r") as fr:
        lines = fr.readlines()

    # remove tmp files
    remove_command = f"rm -f {outputf} {tmp_walkf} {outputf}.txt {CURRENT_DIR}/HetSG"
    os.system(remove_command)

    embeddings = {}
    for line in lines[1:]:
        line_info = line.strip().split(" ")
        if line_info[0] == "</s>":
            pass
        else:
            nodeid = line_info[0]
            if use_hetsg:
                embedding = np.array([float(i) for i in line_info[1:]], dtype=np.float32)
                embeddings[id2node[nodeid[1:]]] = embedding
            else:
                embedding = np.array([float(i) for i in line_info[1:]], dtype=np.float32)
                embeddings[id2node[nodeid]] = embedding

    return embeddings
