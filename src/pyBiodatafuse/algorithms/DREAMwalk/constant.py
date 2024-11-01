"""File paths for DREAMwalk algorithm."""

import os

DDA_DIRECTORY = "dda/"
NODE_TYPE_FILE = "nodetypes.tsv"
GRAPH_FILE = "graph.txt"

ATC_HIERARCHY_FILE = "atc_hierarchy.csv"
DRUG_HIERACHY_FILE = "drug_hierarchy.csv"

DISEASE_SIM_FILE = "disease_similarity.tsv"
DRUG_SIM_FILE = "drug_similarity.tsv"
SIM_FILE = "similarty_graph.txt"

SUBGRAPH_FILE = "subgraph_graph.gml"

EMBEDDING_FILE = "embedding_file.pkl"

CURRENT_DIR = os.path.dirname(os.path.realpath(__file__))
