#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 30 10:21:39 2024

@author: alejandroadriaquelozano
"""

# coding: utf-8

"""Python file for queriying the MINERVA platform (https://minerva.pages.uni.lu/doc/)."""

from typing import Tuple

import pandas as pd
import requests

from pyBiodatafuse.utils import collapse_data_sources, get_identifier_of_interest

# URL of MINERVA's API endpoint



def get_version_minerva() -> dict:
    base_url = "https://covid19map.elixir-luxembourg.org/minerva/api/" #
    """Get version of minerva API.

    :returns: a string containing the version information
    """
    response = requests.get(base_url + "/configuration/")
    conf_dict = response.json()

    return conf_dict["version"]

def list_projects()-> pd.DataFrame:
    """Get information about MINERVA projects.
    
    :returns: a DataFrame containing url, names, and IDs from the different projects in minerva plattform
    """
    base_url = "https://minerva-net.lcsb.uni.lu/api/"
    
    response = requests.get(base_url + "/machines/")
    projects = response.json()
    projects_ids= projects['pageContent']
    project_df=pd.DataFrame()
    for x in projects_ids:
        entry = {'url':x['rootUrl'],'id': x['id']}
        entry_df = pd.DataFrame([entry])
        project_df = pd.concat([project_df, entry_df], ignore_index=True)
        
    map_id_list=[]
    names_list=[]
    for x in project_df['id']:
        x=str(x)
        map_id= requests.get(base_url + "/machines/"+ x+"/projects/").json()['pageContent'][0]['projectId']
        name = requests.get(base_url + "/machines/"+ x+"/projects/").json()['pageContent'][0]['mapName']
        map_id_list.append(map_id)
        names_list.append(name)
    project_df['map_id']=map_id_list
    project_df['names']=names_list
    return project_df

def get_minerva_components(project_df ,map_name,get_elements=True, get_reactions=True) -> dict:
    """Get information about MINERVA componenets from a specific project.
    
    :param project_df: dataframe containing information about all projects contained in Minerva plattform, it is the output from the list_projects() function
    :param map_name: name of the map you want to retrieve the information from. At the moment the options are: 'Asthma Map' 'COVID19 Disease Map' 'Expobiome Map' 'Atlas of Inflammation Resolution' 'SYSCID map' 'Aging Map' 'Meniere's disease map' 'Parkinson's disease map' 'RA-Atlas'   
    :param get_elements: if get_elements = True, the elements of the model will appear as a dictionary in the output of the function
    :param get_reactions: if get_reactions = True, the reactions of the model will appear as a dictionary in the output of the function

    :returns: a Dictionary containing two other dictionaries (map_elements and map_reactions) and a list (models).
        - 'map_elements' contains a list for each of the pathways in the model. Those lists provide information about Compartment,Complex, Drug, Gene, Ion,Phenotype, Protein,RNA and Simple molecules involved in that pathway
        - 'map_reactions' contains a list for each of the pathways in the model. Those lists provide information about the reactions involed in that pathway.
        - 'models' is a list containing pathway-specific information for each of the pathways in the model
    """
    
    #Get url from the project specified
    condition = (project_df['names'] == map_name)
    row= project_df.index[condition].tolist()
    url= project_df.loc[row,'url'].to_string(index=False, header=False)
    project_id = project_df.loc[row,'map_id'].to_string(index=False, header=False)
    


    # Request project data using the extracted project ID
    response = requests.get(url + "/api/projects/" + project_id + "/models/")
    models = response.json() # pull down only models and then iterate over them to extract element of interest
    map_components = {"models": models}

    if get_elements:
        # Get elements of the chosen diagram
        model_elements = {}
        for model in models:
            model = str(model["idObject"])
            url_complete = (
                url
                + "api/projects/"
                + project_id
                + "/models/"
                + model
                + "/"
                + "bioEntities/elements/"
            )
            response_data = requests.get(url_complete)
            model_elements[model] = response_data.json()
        map_components["map_elements"] = model_elements

    if get_reactions:
        # Get reactions of the chosen diagram
        model_reactions = {}
        for model in models:
            model = str(model["idObject"])
            url_complete = (
                url
                + "api/projects/"
                + project_id
                + "/models/"
                + model
                + "/"
                + "bioEntities/reactions/"
            )
            response_data = requests.get(url_complete)
            model_reactions[model] = response_data.json()
        map_components["map_reactions"] = model_reactions

    return map_components


def get_gene_minerva_pathways(
    bridgedb_df: pd.DataFrame,
    map_components: pd.DataFrame,
    input_type: str = "Protein",
) -> Tuple[pd.DataFrame, dict]:
    """Get information about MINERVA pathways associated with a gene.

    :param bridgedb_df: BridgeDb output for creating the list of gene ids to query
    :param map_components: output of the function 'get_minerva_components'
    :param input_type: 'Compartment','Complex', 'Drug', 'Gene', 'Ion','Phenotype','Protein','RNA','Simple molecule'

    :returns: a DataFrame containing DataFrame containing the MINERVA output.
    """
   
    map_elements = map_components.get("map_elements", {})
    models = map_components.get("models", {})

    data_df = get_identifier_of_interest(bridgedb_df, "NCBI Gene")

    names = []
    for value in models:
        name = value["name"]
        names.append(name)

    row = 1
    combined_df = pd.DataFrame()
    for x in names:
        index_to_extract = row
        row = 1 + row

        list_at_index = list(map_elements.values())[index_to_extract - 1]
        common_keys = ["type", "references", "symbol", "name"]
        # Initialize empty lists to store values for each common key
        type = []
        refs = []
        symbol = []
        name = []

        # Iterate through the list of dicts
        for d in list_at_index:
            for key in common_keys:
                if key in d:
                    if key == "type":
                        type.append(d[key])
                    elif key == "references":
                        refs.append(d[key])
                    elif key == "symbol":
                        symbol.append(d[key])
                    elif key == "name":
                        name.append(d[key])

        data = pd.DataFrame()
        data["symbol"] = symbol
        data["pathwayLabel"] = x
        data["pathwayGeneCount"] = len(symbol) - symbol.count(None)
        data["pathwayId"] = models[index_to_extract - 1]["idObject"]
        data["refs"] = refs
        data["type"] = type

        combined_df = pd.concat([combined_df, data], ignore_index=True)
        combined_df = combined_df[combined_df["type"] == input_type]

    if "symbol" not in combined_df:
        return pd.DataFrame()
    else:
        # Add MINERVA output as a new column to BridgeDb file
        combined_df.rename(columns={"symbol": "identifier"}, inplace=True)
        combined_df["identifier"] = combined_df["identifier"].values.astype(str)

        selected_columns = ["pathwayId", "pathwayLabel", "pathwayGeneCount"]

        # Merge the two DataFrames based on 'geneid', 'gene_symbol', 'identifier', and 'target'
        merged_df = collapse_data_sources(
            data_df=data_df,
            source_namespace="NCBI Gene",
            target_df=combined_df,
            common_cols=["identifier"],
            target_specific_cols=selected_columns,
            col_name="MINERVA",
        )

        # Remove duplicates
        minerva_colum = dict(merged_df.MINERVA)
        new_minerva_colum = []
        for key in minerva_colum:
            current_list = minerva_colum[key]

            one_gene = []
            for item in current_list:
                item = str(item)
                one_gene.append(item)
            unique_list = list(set(one_gene))
            new_minerva_colum.append(unique_list)
            
        #
        import ast
        row=1
        for x in new_minerva_colum:
            if 'nan' in x[0]:
                new_minerva_colum[row-1]=[]
            else:
                string_dict= x[0]
                dic = ast.literal_eval(string_dict)
                dic_list = [dic]
                new_minerva_colum[row-1]= dic_list
            row= 1  + row
            
        merged_df["MINERVA"] = new_minerva_colum
        
        # remove 
        #new_minerva_colum_2 =[]
        #or x in new_minerva_colum:
          #  x= x[0].replace('"','')
           # new_minerva_colum_2.append(x)
            
        #merged_df["MINERVA"] = new_minerva_colum_2

            
    return merged_df




