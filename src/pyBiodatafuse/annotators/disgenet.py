# coding: utf-8

"""Python file for queriying DisGeNet database (https://www.disgenet.org/home/)."""

import datetime
import logging
import os
import warnings
from string import Template
from typing import Tuple

import pandas as pd
from SPARQLWrapper import JSON, SPARQLWrapper
from SPARQLWrapper.SPARQLExceptions import SPARQLWrapperException
import requests

from pyBiodatafuse.constants import (
    DISGENET,
    DISGENET_ENDPOINT,
    DISGENET_INPUT_ID,
    DISGENET_OUTPUT_DICT,
)
from pyBiodatafuse.utils import (
    check_columns_against_constants,
    collapse_data_sources,
    get_identifier_of_interest,
)

logger = logging.getLogger("disgenet")


def check_endpoint_disgenet(api_key: str) -> bool:
    """Check the availability of the DisGeNET API.

    :returns: True if the endpoint is available, False otherwise.
    """
    #Set HTTP headers
    HTTPheadersDict = {}
    HTTPheadersDict['Authorization'] = api_key
    HTTPheadersDict['accept'] = 'application/json'
    #Set the DisGeNET API
    s = requests.Session()
    #Get version
    response = s.get("https://api.disgenet.com/api/v1/public/version", headers=HTTPheadersDict)
    # Check if API is down
    if response.json()["status"] == "OK":
        return True
    else:
        return False


def get_version_disgenet(api_key: str) -> dict:
    """Get version of DisGeNET API.

    :returns: a dictionary containing the version information
    """
    #Set HTTP headers
    HTTPheadersDict = {}
    HTTPheadersDict['Authorization'] = api_key
    HTTPheadersDict['accept'] = 'application/json'
    #Set the DisGeNET API
    s = requests.Session()
    #Get version
    version_response = s.get("https://api.disgenet.com/api/v1/public/version", headers=HTTPheadersDict)
    disgenet_version = version_response.json()["payload"]

    return disgenet_version


def get_gene_disease(api_key: str, bridgedb_df: pd.DataFrame) -> Tuple[pd.DataFrame, dict]:
    """Query gene-disease associations from DisGeNET.

    :param api_key: DisGeNET API key (more details can be found at https://disgenet.com/plans)
    :param bridgedb_df: BridgeDb output for creating the list of gene ids to query.
    :returns: a DataFrame containing the DisGeNET output and dictionary of the DisGeNET metadata.
    """
    # Check if the DisGeNET API is available
    api_available = check_endpoint_disgenet(api_key)

    if not api_available:
        warnings.warn(
            f"{DISGENET} endpoint is not available. Unable to retrieve data.", stacklevel=2
        )
        return pd.DataFrame(), {}
    
    
    # Add the API key to the requests headers
    HTTPheadersDict = {}
    HTTPheadersDict['Authorization'] = API_KEY
    HTTPheadersDict['accept'] = 'application/json'
    #Set the DisGeNET API
    s = requests.Session()

    # Extract the "target" values and join them into a single string separated by commas
    data_df = get_identifier_of_interest(bridgedb_df, "NCBI Gene")
    disgenet_input = ",".join(data_df["target"])

    # Record the start time
    start_time = datetime.datetime.now()

    # Split the targets into chunks of 99 or fewer
    targets_list = disgenet_input.split(",")
    chunk_size = 99
    chunks = [targets_list[i : i + chunk_size] for i in range(0, len(targets_list), chunk_size)]

    disgenet_output = []

    if not params:
        params = {"source": "CURATED", "format": "json"}
    else:
        params["format"] = "json"
        params["source"] = "CURATED"

    for chunk in chunks:
        # Join the chunked targets into a comma-separated string
        chunked_input = ",".join(chunk)
        # Get all the diseases associated with genes for the current chunk
        gda_response = s.get(f"{api_host}/gda/gene/{chunked_input}", params=params)
        chunk_output = gda_response.json()
        disgenet_output.extend(chunk_output)

    # Record the end time
    end_time = datetime.datetime.now()

    # Convert disgenet_output to a DataFrame
    disgenet_df = pd.DataFrame(disgenet_output)
    if "geneid" not in disgenet_df:
        return pd.DataFrame()
    else:
        # Drop the uniprotid column
        disgenet_df.drop("uniprotid", axis=1, inplace=True)
        # Add DisGeNET output as a new column to BridgeDb file
        disgenet_df.rename(columns={"geneid": "target"}, inplace=True)
        disgenet_df["target"] = disgenet_df["target"].values.astype(str)

        selected_columns = [
            "gene_dsi",
            "gene_dpi",
            "gene_pli",
            "protein_class",
            "protein_class_name",
            "diseaseid",
            "disease_name",
            "disease_class",
            "disease_class_name",
            "disease_type",
            "disease_semantic_type",
            "score",
            "ei",
            "el",
            "year_initial",
            "year_final",
            "source",
        ]

        # Merge the two DataFrames based on 'geneid', 'gene_symbol', 'identifier', and 'target'
        merged_df = collapse_data_sources(
            data_df=data_df,
            source_namespace="NCBI Gene",
            target_df=disgenet_df,
            common_cols=["target"],
            target_specific_cols=selected_columns,
            col_name="DisGeNET",
        )

        """Metdata details"""
        # Get the current date and time
        current_date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        # Calculate the time elapsed
        time_elapsed = str(end_time - start_time)
        # Add version to metadata file
        disgenet_version = get_version_disgenet(api_host=api_host)
        # Add the datasource, query, query time, and the date to metadata
        disgenet_metadata = {
            "datasource": "DisGeNET",
            "metadata": {"source_version": disgenet_version},
            "query": {
                "size": len(disgenet_input.split(",")),
                "input_type": "NCBI Gene",
                "time": time_elapsed,
                "date": current_date,
                "url": gda_response.request.url,
            },
        }

        if s:
            s.close()

        return merged_df, disgenet_metadata
    # Extract the "target" values and join them into a single string separated by commas
    data_df = get_identifier_of_interest(bridgedb_df, DISGENET_INPUT_ID)
    gene_list = data_df["target"].tolist()
    gene_list = list(set(gene_list))

    query_gene_lists = []
    if len(gene_list) > 25:
        for i in range(0, len(gene_list), 25):
            tmp_list = gene_list[i : i + 25]
            query_gene_lists.append(" ".join(f'"{g}"' for g in tmp_list))

    else:
        query_gene_lists.append(" ".join(f'"{g}"' for g in gene_list))

    with open(os.path.dirname(__file__) + "/queries/disgenet-genes-disease.rq", "r") as fin:
        sparql_query = fin.read()

    # Record the start time
    disgenet_version = get_version_disgenet()
    start_time = datetime.datetime.now()

    sparql = SPARQLWrapper(DISGENET_ENDPOINT)
    sparql.setReturnFormat(JSON)

    query_count = 0

    intermediate_df = pd.DataFrame()

    for gene_list_str in query_gene_lists:
        query_count += 1
        sparql_query_template = Template(sparql_query)
        substit_dict = dict(gene_list=gene_list_str)
        sparql_query_template_sub = sparql_query_template.substitute(substit_dict)
        sparql.setQuery(sparql_query_template_sub)
        res = sparql.queryAndConvert()
        res = res["results"]["bindings"]
        df = pd.DataFrame(res)
        df = df.applymap(lambda x: x["value"])

        intermediate_df = pd.concat(
            [intermediate_df, df], ignore_index=True
        )  # this is also adding to the time

    # Record the end time
    end_time = datetime.datetime.now()

    # Organize the annotation results as an array of dictionaries
    if "gene_id" not in intermediate_df:
        return pd.DataFrame(), {"datasource": DISGENET, "metadata": disgenet_version}
    intermediate_df.drop_duplicates(inplace=True)
    intermediate_df["target"] = intermediate_df["gene_id"].apply(lambda x: x.split("/")[-1])
    intermediate_df["disease_id"] = intermediate_df["description"].apply(
        lambda x: "umls:" + x.split("umls:")[1].split("]")[0].strip()
    )
    intermediate_df["disease_name"] = intermediate_df["description"].apply(
        lambda x: x.split("[umls:")[0].strip()
    )
    intermediate_df["score"] = intermediate_df["disease_score"].astype(float)
    intermediate_df["evidence_source"] = intermediate_df["evidence_source"].apply(
        lambda x: x.split("/")[-1]
    )

    intermediate_df["target"] = intermediate_df["target"].values.astype(str)

    intermediate_df = intermediate_df[
        ["target", "disease_id", "disease_name", "score", "evidence_source"]
    ]

    # Check if all keys in df match the keys in OUTPUT_DICT
    check_columns_against_constants(
        data_df=intermediate_df,
        output_dict=DISGENET_OUTPUT_DICT,
        check_values_in=["disease_id"],
    )

    merged_df = collapse_data_sources(
        data_df=bridgedb_df,
        source_namespace=DISGENET_INPUT_ID,
        target_df=intermediate_df,
        common_cols=["target"],
        target_specific_cols=list(DISGENET_OUTPUT_DICT.keys()),
        col_name=DISGENET,
    )

    """Metdata details"""
    # Get the current date and time
    current_date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    # Calculate the time elapsed
    time_elapsed = str(end_time - start_time)
    # Add version, datasource, query, query time, and the date to metadata
    disgenet_metadata = {
        "datasource": DISGENET,
        "metadata": disgenet_version,
        "query": {
            "size": len(gene_list),
            "input_type": DISGENET_INPUT_ID,
            "time": time_elapsed,
            "date": current_date,
            "url": DISGENET_ENDPOINT,
        },
    }

    return merged_df, disgenet_metadata
