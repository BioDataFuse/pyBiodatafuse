import argparse
import pandas as pd
import os
import pickle

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--knowledge_graph_file', type=str, required=True)
    parser.add_argument('--embeddingf', type=str, required=True)
    parser.add_argument('--model_folder', type=str, required=True)
    #parser.add_argument('--output_file', type=str, required=True)
    
    parser.add_argument('--query_disease', type=str)
    parser.add_argument('--candidates_count', type=int, default=10)
    
    args=parser.parse_args()
    args={'kgfile':args.knowledge_graph_file,
     'embeddingf':args.embeddingf,
     'model_folder':args.model_folder,
     'query_disease':args.query_disease,
     #'outputf':args.output_file,
     'candidates_count':args.candidates_count}
    
    return args

# function to map _id to the corresponding identifier based on _labels
def map_id(row):
    if row['_labels'] == ':Drug':
        return row['drugID']
    elif row['_labels'] == ':Disease':
        return row['diseaseID']
    else:
        return None
    
# function to process unlabeled drugs - returns info of candidate drugs
def process_drugs(drugs_to_process, query_disease, embeddingf, model_list, map_name_dict, candidates_count):
    with open(embeddingf,'rb') as fin:
        embedding_dict = pickle.load(fin)
    result_df = pd.DataFrame(columns=['drug', 'name', 'avg_prob'])
    for drug in drugs_to_process:
        prob_sum = 0
        for i in range(0,10):
            prob_sum += model_list[i].predict_proba([embedding_dict[drug]-embedding_dict[query_disease]])[:, 1]
        #model_files = os.listdir(model_folder)
        #for file in model_files:
        #    file_path = os.path.join(model_folder, file)
        #    with open(file_path, 'rb') as f:
        #        model = pickle.load(f)     
        #    prob_sum += model.predict_proba([embedding_dict[drug]-embedding_dict[query_disease]])[:, 1]
            #print(prob_sum)
        prob_avg = prob_sum / 10
        result_df.loc[len(result_df) + 1] = {'drug': drug, 'name': map_name_dict[drug], 'avg_prob': prob_avg}
    candidates = result_df.sort_values(by='avg_prob', ascending=False).head(candidates_count)

    return candidates

#def find_candidates(kgfile:str, model_folder:str, outputf:str, candidates_count:int=10):
def find_candidates(kgfile, embeddingf, model_folder, query_disease, candidates_count:int=10):
    kg_data = pd.read_csv(kgfile, dtype=str)

    # filter rows with ':Drug' or ':Disease' in '_labels'
    nodes_filtered = kg_data[kg_data['_labels'].isin([':Drug', ':Disease'])].copy()
    # create a new dataframe for neo4j id -> actual id
    id_map = pd.DataFrame({
        'id': nodes_filtered['_id'],
        'label': nodes_filtered['_labels'],
        'name': nodes_filtered['name'],
        'mapped_id': nodes_filtered.apply(map_id, axis=1)
    })

    # create a dictionary mapping id to mapped_id
    id_map_dict = dict(zip(id_map['id'], id_map['mapped_id']))

    # create a dictionary mapping mapped_id to name
    map_name_dict = dict(zip(id_map['mapped_id'], id_map['name']))

    # filter rows where _label is equal to "Drug" - find all drugs
    drug_array = id_map.loc[id_map['label'] == ':Drug', 'mapped_id'].values    

    # filter rows with 'INDICATES' in '_type'
    drug_disease_ints_filtered = kg_data[kg_data['_type'].isin(['INDICATES'])].copy()

    # find drugs that are assocaiated with the query disease
    associated_drugs = []
    for index, row in drug_disease_ints_filtered.iterrows():
        if id_map_dict.get(row['_end']) == query_disease:
            associated_drugs.append(id_map_dict.get(row['_start'])) 

    # find drugs that are not assocaiated with the query disease
    temp = set(associated_drugs)
    drugs_to_process = [x for x in drug_array if x not in temp]

    # create a model list from the model files
    model_files = os.listdir(model_folder)
    model_list = []
    for file in model_files:
        file_path = os.path.join(model_folder, file)

        with open(file_path, 'rb') as f:
            model = pickle.load(f)

        # Process content as needed
        model_list.append(model)           

    candidate_drugs = process_drugs(drugs_to_process, query_disease, embeddingf, model_list, map_name_dict, candidates_count)
    print(candidate_drugs)

if __name__ == '__main__':
    args=parse_args()
    find_candidates(**args)