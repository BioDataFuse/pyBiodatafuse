#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Codes to run Dreamwalk algorithm on BDF graph."""


import logging

import pyBiodatafuse.algorithms.DREAMwalk.generate_dis_sim as dis_gen
import pyBiodatafuse.algorithms.DREAMwalk.generate_files as gen
from pyBiodatafuse.algorithms.DREAMwalk.calculate_drug_scores import find_candidates
from pyBiodatafuse.algorithms.DREAMwalk.generate_embeddings import save_embedding_files
from pyBiodatafuse.algorithms.DREAMwalk.generate_similarity_net import save_sim_graph
from pyBiodatafuse.algorithms.DREAMwalk.predict_associations import predict_dda
from pyBiodatafuse.analyzer.summarize import BioGraph

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def preprocess_files():
    """Generating base files for the algorithm."""
    graph_obj = BioGraph(graph_path="PCS_graph.gml", graph_format="gml")
    gen.create_files(graph_obj=graph_obj, output_dir="dreamwalk_data")
    logger.info("Files generated successfully.")


def generate_sim():
    """Generate similarity graph."""
    hierf = "dreamwalk_data/hierarchy.csv"
    simf = "dreamwalk_data/similarty_graph_drugs.tsv"
    cutoff = 0.4

    graph_obj = BioGraph(graph_path="dreamwalk_data/subgraph_graph.gml", graph_format="gml")
    save_sim_graph(graph_obj=graph_obj, hierf=hierf, outputf=simf, cutoff=cutoff)
    logger.info("Similarity graph generated successfully.")


# # create similarity graph merging disease and drugs
# dis_sim = pd.read_csv("dis_sim.tsv", sep="\t", header=None)
# drug_sim = pd.read_csv("similarty_graph_drugs.tsv", sep="\t", header=None)

# drug_sim[2] = 2
# dis_sim[2] = 1

# similarty_graph = pd.concat([dis_sim, drug_sim], ignore_index=True)

# similarty_graph[4] = range(0, len(similarty_graph))
# similarty_graph.to_csv("similarty_graph.txt", sep="\t", index=False, header=False)

"""
# Generate node embeddings by teleport-guided randomwalk

networkf='graph.txt'
hierf='hierarchy.csv'
nodetypef='nodetypes.tsv'
embeddingf='embedding_file.pkl'
simf='similarty_graph.txt'

save_embedding_files(netf=networkf,sim_netf=simf, outputf=embeddingf,
                    nodetypef=nodetypef,tp_factor=0.3)


# Predict drug-disease association

pairf='dda_files/dda1.tsv'
modelf='results/clf1.pkl'
embeddingf='embedding_file.pkl'

predict_dda(embeddingf=embeddingf, pairf=pairf, modelf=modelf)


# Calculate drug scores
embeddingf='embedding_file.pkl'
model_folder= 'results'
query_disease= 'C00000'
kgfile='preprocessed_graph.csv'
find_candidates(kgfile, embeddingf, model_folder, query_disease)

results1=find_candidates(kgfile, embeddingf, model_folder, query_disease, candidates_count=400)
results1.to_csv('results.csv')
"""

if __name__ == "__main__":
    preprocess_files()
    generate_sim()
    # generate_embedding_files()
    # predict_dda()
    # find_candidates()
