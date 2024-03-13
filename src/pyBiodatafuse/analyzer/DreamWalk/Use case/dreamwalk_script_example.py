#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 20 10:32:58 2024

@author: alejandroadriaquelozano
"""


# Set path
import os 
from pathlib import Path
path = Path(__file__).resolve().parent
os.chdir(path)

import pandas as pd
import DREAMwalk.generate_files as gen
import DREAMwalk.generate_dis_sim as dis_gen
from DREAMwalk.generate_embeddings import save_embedding_files
from DREAMwalk.predict_associations import predict_dda
from DREAMwalk.generate_similarity_net import save_sim_graph
from DREAMwalk.calculate_drug_scores import find_candidates


# GENERSTE FILES
kg_data= pd.read_csv('post-COVID_graph.csv')


gen.generate_files(kg_data)

dis_gen.save_dis_sim('post-COVID_graph.csv','dis_sim.csv')


networkf='graph.txt'
hierf='hierarchy.csv'
simf='similarty_graph_drugs.tsv'
cutoff=0.4

save_sim_graph(networkf=networkf,hierf=hierf,outputf=simf,cutoff=cutoff)

#create similarity graph merging disease and drugs
dis_sim = pd.read_csv('dis_sim.tsv',sep='\t',header=None)
drug_sim = pd.read_csv('similarty_graph_drugs.tsv',sep='\t',header=None)

similarty_graph=pd.concat([dis_sim,drug_sim],ignore_index=True)

similarty_graph[4] = range(0, len(similarty_graph))
similarty_graph.to_csv('similarty_graph.txt', sep='\t', index=False,header=False)



# Generate node embeddings by teleport-guided randomwalk

nodetypef='nodetypes.tsv'
embeddingf='embedding_file.pkl'
simf='similarty_graph.txt'

save_embedding_files(netf=networkf,sim_netf=hierf, outputf=embeddingf,
                    nodetypef=nodetypef,tp_factor=0.3)


# Predict drug-disease association

pairf='dda_files/dda10.tsv'
modelf='results/clf10.pkl'
embeddingf='embedding_file.pkl'

predict_dda(embeddingf=embeddingf, pairf=pairf, modelf=modelf)

# Calculate drug scores, run only this section to just check the results (running the algorithm is quite time consuming, especially the embedding file)
model_folder= 'results'
query_disease= 'C0000000' #Post COVID 19 node ID
kgfile='preprocessed_graph.csv'
find_candidates(kgfile, embeddingf, model_folder, query_disease, candidates_count=20)

