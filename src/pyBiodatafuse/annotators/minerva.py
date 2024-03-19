#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Python file for queriying the MINERVA platform (https://minerva.pages.uni.lu/doc/)."""

import datetime
import warnings
from typing import Optional, Tuple

import pandas as pd
import requests

from pyBiodatafuse.constants import MINERVA, MINERVA_ENDPOINT, MINERVA_INPUT_ID, MINERVA_OUTPUT_DICT
from pyBiodatafuse.utils import collapse_data_sources, get_identifier_of_interest, check_columns_against_constants


def check_endpoint_minerva() -> bool:
    """Check the availability of the MINERVA API endpoint.

    :returns: True if the endpoint is available, False otherwise.
    """
    response = requests.get(f"{MINERVA_ENDPOINT}/machines/")

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
    base_endpoint = f"{MINERVA_ENDPOINT}/machines/"
    response = requests.get(base_endpoint).json()
    projects = response
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
        - 'models' is a list containing pathway-specific information for each of the pathways in the model
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


def get_gene_minerva_pathways(
    bridgedb_df: pd.DataFrame,
    map_name: str,
    input_type: Optional[str] = "Protein",
    get_elements: Optional[bool] = True,
    get_reactions: Optional[bool] = True,
) -> Tuple[pd.DataFrame, dict]:
    """Get information about MINERVA pathways associated with a gene.

    :param bridgedb_df: BridgeDb output for creating the list of gene ids to query
    :param map_name: name of the map you want to retrieve the information from. The extensive list
        can be found at https://minerva-net.lcsb.uni.lu/table.html.
    :param input_type: type of input gene. Default is "Protein"
    :param get_elements: boolean to get elements of the chosen diagram
    :param get_reactions: if get_reactions = boolean to get reactions of the chosen diagram
    :returns: a tuple containing MINERVA outputs and dictionary of the MINERVA metadata.
    """
    # Check if the MINERVA API is available
    api_available = check_endpoint_minerva()
    if not api_available:
        warnings.warn(
            f"{MINERVA} API endpoint is not available. Unable to retrieve data.", stacklevel=2
        )
        return pd.DataFrame(), {}

    assert input_type in [
        "Compartment",
        "Complex",
        "Drug",
        "Gene",
        "Ion",
        "Phenotype",
        "Protein",
        "RNA",
        "Simple molecule",
    ], "Incorrect Input_type provided. Please provide a valid input_type."

    # Record the start time
    start_time = datetime.datetime.now()

    map_url, map_components = get_minerva_components(
        map_name=map_name, get_elements=get_elements, get_reactions=get_reactions
    )
    map_elements = map_components.get("map_elements", {})
    models = map_components.get("models", {})

    data_df = get_identifier_of_interest(bridgedb_df, MINERVA_INPUT_ID)

    names = []
    for value in models:
        name = value["name"]
        names.append(name)

    row = 1
    intermediate_df = pd.DataFrame()
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
        data["pathway_label"] = x
        data["pathway_gene_count"] = len(symbol) - symbol.count(None)
        data["pathway_id"] = models[index_to_extract - 1]["idObject"]
        data["refs"] = refs
        data["type"] = type

        intermediate_df = pd.concat([intermediate_df, data], ignore_index=True)
        intermediate_df = intermediate_df[intermediate_df["type"] == input_type]

    # Record the end time
    end_time = datetime.datetime.now()

    if "symbol" not in intermediate_df:
        return pd.DataFrame(), {}

    # Organize the annotation results as an array of dictionaries
    # TODO the merge is based on the gene symbol, what if another id is being used as input
    intermediate_df.rename(columns={"symbol": "identifier"}, inplace=True)
    intermediate_df["identifier"] = intermediate_df["identifier"].values.astype(str)

    intermediate_df = intermediate_df.drop_duplicates(subset=["identifier", "pathway_id", "pathway_label", "pathway_gene_count"])

    # Check if all keys in df match the keys in OUTPUT_DICT
    check_columns_against_constants(
        data_df=intermediate_df,
        output_dict=MINERVA_OUTPUT_DICT,
        check_values_in=["pathway_id"],
    )

    # Merge the two DataFrames on the target column
    merged_df = collapse_data_sources(
        data_df=data_df,
        source_namespace=MINERVA_INPUT_ID,
        target_df=intermediate_df,
        common_cols=["identifier"],
        target_specific_cols=list(MINERVA_OUTPUT_DICT.keys()),
        col_name=MINERVA,
    )

    """Metdata details"""
    # Get the current date and time
    current_date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    # Calculate the time elapsed
    time_elapsed = str(end_time - start_time)
    # Add version to metadata file
    minerva_version = get_version_minerva(map_endpoint=map_url)
    # Add the datasource, query, query time, and the date to metadata
    minerva_metadata = {
        "datasource": MINERVA,
        "metadata": minerva_version,
        "query": {
            "size": data_df["target"].nunique(),
            "input_type": MINERVA_INPUT_ID,
            "MINERVA project": map_name,
            "time": time_elapsed,
            "date": current_date,
            "url": map_url,
        },
    }

    return merged_df, minerva_metadata
