#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Python file for queriying the MINERVA platform (https://minerva.pages.uni.lu/doc/)."""

import datetime
import logging
import warnings
from typing import Any, Dict, Optional, Tuple

import pandas as pd
import requests
from tqdm import tqdm

import pyBiodatafuse.constants as Cons
from pyBiodatafuse.utils import (
    check_columns_against_constants,
    collapse_data_sources,
    get_identifier_of_interest,
    give_annotator_warning,
)

logger = logging.getLogger(__name__)


def check_endpoint_minerva() -> bool:
    """Check the availability of the MINERVA API endpoint.

    :returns: True if the endpoint is available, False otherwise.
    """
    response = requests.get(f"{Cons.MINERVA_ENDPOINT}/machines/")

    # Check if API is down
    if response.status_code == 200:
        return True
    else:
        return False


def get_version_minerva(map_endpoint: str) -> dict:
    """Get version of minerva API.

    :param map_endpoint: MINERVA map API endpoint (eg. "https://covid19map.elixir-luxembourg.org/minerva/")
    :returns: a dictionary containing the version information
    """
    response = requests.get(map_endpoint + "api/configuration/")

    conf_dict = response.json()
    minerva_version = {"source_version": conf_dict["version"]}

    return minerva_version


def list_projects() -> pd.DataFrame:
    """Get information about MINERVA projects.

    :returns: a dataFrame containing url, names and IDs from the different projects in MINERVA plattform
    """
    base_endpoint = f"{Cons.MINERVA_ENDPOINT}/machines/"
    projects = requests.get(base_endpoint).json()
    projects_ids = projects["pageContent"]

    project_df = pd.DataFrame()

    for x in projects_ids:
        entry = {"url": x["rootUrl"], "id": x["id"]}
        entry_df = pd.DataFrame([entry])
        project_df = pd.concat([project_df, entry_df], ignore_index=True)

    map_id_list = []
    names_list = []
    for x in project_df["id"]:
        x = str(x)
        if len(requests.get(f"{base_endpoint}/{x}/projects/").json()["pageContent"]) != 0:
            map_id = requests.get(f"{base_endpoint}/{x}/projects/").json()["pageContent"][0][
                "projectId"
            ]
            name = requests.get(f"{base_endpoint}/{x}/projects/").json()["pageContent"][0][
                "mapName"
            ]
            map_id_list.append(map_id)
            names_list.append(name)
        else:
            project_df = project_df[
                project_df["id"] != int(x)
            ]  # If pageContent is not present, then delete this entry

    project_df["map_id"] = map_id_list
    project_df["names"] = names_list

    return project_df


def get_minerva_components(
    map_name: str,
    get_elements: Optional[bool] = True,
    get_reactions: Optional[bool] = True,
) -> Tuple[str, dict]:
    """Get information about MINERVA componenets from a specific project.

    :param map_name: MINERVA map name. The extensive list can be found at https://minerva-net.lcsb.uni.lu/table.html.
    :param get_elements: boolean to get elements of the chosen diagram
    :param get_reactions: boolean to get reactions of the chosen diagram
    :returns: a tuple of map endpoint and dictionary containing:
        - 'map_elements' contains a list for each of the pathways in the model.
        Those lists provide information about Compartment, Complex, Drug, Gene, Ion, Phenotype,
        Protein, RNA and Simple molecules involved in that pathway
        - 'map_reactions' contains a list for each of the pathways in the model.
        Those lists provide information about the reactions involed in that pathway.
        - 'models' is a list containing pathway-specific information for each of the pathways in the model.

    """
    # Get list of projects
    project_df = list_projects()

    # Get url from the project specified
    condition = project_df["names"] == map_name
    row = project_df.index[condition].tolist()
    map_url = project_df.loc[row, "url"].to_string(index=False, header=False)
    project_id = project_df.loc[row, "map_id"].to_string(index=False, header=False)

    # Request project data using the extracted project ID
    response = requests.get(map_url + "/api/projects/" + project_id + "/models/")

    models = (
        response.json()
    )  # pull down only models and then iterate over them to extract element of interest
    map_components = {"models": models}

    if get_elements:
        # Get elements of the chosen diagram
        model_elements = {}
        for model in models:
            model = str(model["idObject"])
            url_complete = (
                map_url
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
                map_url
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

    return map_url, map_components


def get_gene_pathways(
    bridgedb_df: pd.DataFrame,
    map_name: str,
    get_elements: Optional[bool] = True,
    get_reactions: Optional[bool] = True,
) -> Tuple[pd.DataFrame, dict]:
    """Get information about MINERVA pathways associated with a gene.

    :param bridgedb_df: BridgeDb output for creating the list of gene ids to query
    :param map_name: name of the map you want to retrieve the information from. The extensive list
        can be found at https://minerva-net.lcsb.uni.lu/table.html.
    :param get_elements: boolean to get elements of the chosen diagram.
    :param get_reactions: if get_reactions = boolean to get reactions of the chosen diagram.
    :returns: a tuple containing MINERVA outputs and dictionary of the MINERVA metadata.
    """
    # Check if the MINERVA API is available
    api_available = check_endpoint_minerva()
    if not api_available:
        warnings.warn(
            f"{Cons.MINERVA} API endpoint is not available. Unable to retrieve data.", stacklevel=2
        )
        return pd.DataFrame(), {}

    data_df = get_identifier_of_interest(bridgedb_df, Cons.MINERVA_GENE_INPUT_ID)

    # Record the start time
    start_time = datetime.datetime.now()

    logger.info("Getting minerva components")
    map_url, map_components = get_minerva_components(
        map_name=map_name, get_elements=get_elements, get_reactions=get_reactions
    )
    minerva_version = get_version_minerva(map_endpoint=map_url)

    map_elements = map_components.get("map_elements", {})
    models = map_components.get("models", {})

    names = []
    for value in models:
        name = value[Cons.NAME]
        names.append(name)

    intermediate_df = pd.DataFrame()

    for idx, pathway_name in tqdm(enumerate(names), total=len(names), desc="Processing pathways"):
        pathway_data = list(map_elements.values())[idx]

        # Initialize empty lists to store values for each common key
        entity_type = []
        refs = []
        symbol = []
        ensembl = []

        # Iterate through the list of dicts
        for data in pathway_data:
            for col in Cons.INTERESTED_INFO:
                if col not in data:
                    continue

                value = data[col]

                if col == Cons.ENTITY_TYPE:
                    entity_type.append(value)
                elif col == Cons.ENTITY_SYMBOL:
                    symbol.append(value)
                elif col == Cons.ENTITY_REFS:
                    refs.append(value)

                if isinstance(value, list):
                    ensembl_id = None
                    for p in value:
                        if p[Cons.ENTITY_TYPE].lower() == Cons.ENSEMBL.lower():
                            ensembl_id = p["resource"]
                    ensembl.append(ensembl_id)

        tmp_df = pd.DataFrame()
        tmp_df[Cons.ENTITY_SYMBOL] = symbol
        tmp_df[Cons.PATHWAY_LABEL] = pathway_name
        tmp_df[Cons.PATHWAY_GENE_COUNTS] = len(symbol) - symbol.count(None)
        tmp_df[Cons.PATHWAY_ID] = f"{Cons.MINERVA}:" + str(models[idx]["idObject"])
        tmp_df[Cons.ENTITY_REFS] = refs
        tmp_df[Cons.ENSEMBL] = ensembl
        tmp_df[Cons.ENTITY_TYPE] = entity_type
        pathway_subset = tmp_df[
            tmp_df[Cons.ENTITY_TYPE] == "Protein"
        ]  # to get gene linked information

        intermediate_df = pd.concat([intermediate_df, pathway_subset], ignore_index=True)

    # Record the end time
    end_time = datetime.datetime.now()

    """Metdata details"""
    # Get the current date and time
    current_date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    # Calculate the time elapsed
    time_elapsed = str(end_time - start_time)

    # Add the datasource, query, query time, and the date to metadata
    minerva_metadata: Dict[str, Any] = {
        "datasource": Cons.MINERVA,
        "metadata": minerva_version,
        "query": {
            "size": data_df[Cons.TARGET_COL].nunique(),
            "input_type": Cons.MINERVA_GENE_INPUT_ID,
            "MINERVA project": map_name,
            "time": time_elapsed,
            "date": current_date,
            "url": map_url,
        },
    }

    # Organize the annotation results as an array of dictionaries
    intermediate_df.rename(columns={Cons.ENSEMBL: Cons.TARGET_COL}, inplace=True)
    intermediate_df[Cons.TARGET_COL] = intermediate_df[Cons.TARGET_COL].values.astype(str)

    intermediate_df = intermediate_df.drop_duplicates(
        subset=[Cons.TARGET_COL, Cons.PATHWAY_ID, Cons.PATHWAY_LABEL, Cons.PATHWAY_GENE_COUNTS]
    )
    intermediate_df = intermediate_df[
        intermediate_df[Cons.TARGET_COL].isin(data_df[Cons.TARGET_COL])
    ]

    if intermediate_df.empty:
        warnings.warn(
            f"There is no annotation for your input list in {Cons.MINERVA}, project {map_name}.",
            stacklevel=2,
        )
        return pd.DataFrame(), minerva_metadata

    # Check if all keys in df match the keys in OUTPUT_DICT
    check_columns_against_constants(
        data_df=intermediate_df,
        output_dict=Cons.MINERVA_PATHWAY_OUTPUT_DICT,
        check_values_in=[Cons.MINERVA_PATHWAY_DEFAULT_ID],
    )

    # Merge the two DataFrames on the target column
    merged_df = collapse_data_sources(
        data_df=data_df,
        source_namespace=Cons.MINERVA_GENE_INPUT_ID,
        target_df=intermediate_df,
        common_cols=[Cons.TARGET_COL],
        target_specific_cols=list(Cons.MINERVA_PATHWAY_OUTPUT_DICT.keys()),
        col_name=Cons.MINERVA_PATHWAY_COL,
    )

    """Update metadata"""
    # Calculate the number of new nodes
    num_new_nodes = intermediate_df[Cons.PATHWAY_ID].nunique()
    # Calculate the number of new edges
    num_new_edges = intermediate_df.drop_duplicates(
        subset=[Cons.TARGET_COL, Cons.PATHWAY_ID]
    ).shape[0]

    # Check the intermediate_df
    if num_new_edges != len(intermediate_df):
        give_annotator_warning(Cons.MINERVA)

    # Add the number of new nodes and edges to metadata
    minerva_metadata[Cons.QUERY][Cons.NUM_NODES] = num_new_nodes
    minerva_metadata[Cons.QUERY][Cons.NUM_EDGES] = num_new_edges

    return merged_df, minerva_metadata
