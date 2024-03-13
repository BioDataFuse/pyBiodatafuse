import pandas as pd
import requests
import numpy as np
import json 
from sklearn.metrics import jaccard_score

# function to map _id to the corresponding identifier based on _labels
def map_id(row):
    if row['_labels'] == ':Disease':
        return row['diseaseID']
    else:
        return None
    
def jaccard_similarity(list1, list2):
    set1 = set(list1)
    set2 = set(list2)
    intersection = len(set1.intersection(set2))
    union = len(set1.union(set2))
    return intersection / union    

def save_dis_sim(kgfile, sim_file):
    kg_data = pd.read_csv(kgfile, dtype=str)

    # filter rows with ':Disease' in '_labels'
    disease_nodes = kg_data[kg_data['_labels'].isin([':Disease'])].copy()

    # create a new dataframe for neo4j id -> actual id
    id_map = pd.DataFrame({
        'id': disease_nodes['_id'],
        'name': disease_nodes['name'],
        'mapped_id': disease_nodes.apply(map_id, axis=1)
    })

    diseases = id_map["mapped_id"].unique()

    disease_dict = {}
    for i in range(0, len(diseases)):
        disease_dict[diseases[i]] = []
    
    
    ########################################## Retrieve disease-gene information locally (API not working)
    gene_ids = pd.read_csv('/Users/alejandroadriaquelozano/Documents/Internships/MacsBio/KG/disgenet/genes_disgenet.csv')
    diseases_ids = pd.read_csv('/Users/alejandroadriaquelozano/Documents/Internships/MacsBio/KG/disgenet/diseases_disgenet.csv')
    gene_diseases = pd.read_csv('/Users/alejandroadriaquelozano/Documents/Internships/MacsBio/KG/disgenet/gene-disease_disgenet.csv')
    
    df_all= pd.DataFrame()
    df_all['diseaseId']  =   diseases
    df_all = pd.merge(df_all,diseases_ids,on='diseaseId',how='left')
    df_all = pd.merge(df_all,gene_diseases,on='diseaseNID',how='left')
    df_all = pd.merge(df_all,gene_ids,on='geneNID',how='left')
    
    df_all=df_all[['diseaseId','geneId']]
    df_all['geneId']=df_all['geneId'].astype(str)
    
    grouped = df_all.groupby('diseaseId')['geneId'].agg(list).reset_index()
    
    disease_dict = grouped.set_index('diseaseId').to_dict()['geneId']
    ##########################################
    disease_similarity_matrix = {}

    print('Jaccard similarities are being calculated...')
    for disease1, genes1 in disease_dict.items():
        for disease2, genes2 in disease_dict.items():
            if disease1 != disease2:
                similarity = jaccard_similarity(genes1, genes2)
                disease_similarity_matrix.setdefault(disease1, {})[disease2] = similarity

    print('Similarities are being written...')
    colnames = ['a', 'b', 'c', 'd', 'e'] 
    sim_graph = pd.read_csv(sim_file, names=colnames, sep="\t")
    data = []
    count = len(sim_graph.index)    
    for key, value in disease_similarity_matrix.items():
        for key2, value2 in value.items():
            if value2 > 0.4:
                data.append([key, key2, 3, value2, count] )
                count = count + 1

    df = pd.DataFrame(data, columns = ['a', 'b', 'c', 'd', 'e']) 
    sim_graph =  pd.concat([sim_graph, df], ignore_index=True, sort=False)                

    sim_graph.to_csv(sim_file, sep="\t", index = False, header = False)
    
