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


def check_endpoint_disgenet() -> bool:
    """Check the availability of the DisGeNET SPARQL endpoint.

    :returns: True if the endpoint is available, False otherwise.
    """
    with open(os.path.dirname(__file__) + "/queries/disgenet-metadata.rq", "r") as fin:
        sparql_query = fin.read()

    sparql = SPARQLWrapper(DISGENET_ENDPOINT)
    sparql.setReturnFormat(JSON)

    sparql.setQuery(sparql_query)

    try:
        sparql.queryAndConvert()
        return True
    except SPARQLWrapperException:
        return False


def get_version_disgenet() -> dict:
    """Get version of DisGeNET API.

    :returns: a dictionary containing the version information
    """
    with open(os.path.dirname(__file__) + "/queries/disgenet-metadata.rq", "r") as fin:
        sparql_query = fin.read()

    sparql = SPARQLWrapper(DISGENET_ENDPOINT)
    sparql.setReturnFormat(JSON)

    sparql.setQuery(sparql_query)

    version_response = sparql.queryAndConvert()
    bindings = version_response["results"]["bindings"]
    pattern = "RDF Distribution"
    disgenet_version = {"version": ""}  # Set default value

    for binding in bindings:
        title = binding["title"]["value"]
        if pattern in title:
            disgenet_version = {"source_version": title}

    return disgenet_version


def get_gene_disease(bridgedb_df: pd.DataFrame) -> Tuple[pd.DataFrame, dict]:
    """Query gene-disease associations from DisGeNET.

    :param bridgedb_df: BridgeDb output for creating the list of gene ids to query.
    :returns: a DataFrame containing the DisGeNET output and dictionary of the DisGeNET metadata.
    """
    # Check if the DisGeNET API is available
    api_available = check_endpoint_disgenet()

    if not api_available:
        warnings.warn(
            f"{DISGENET} SPARQL endpoint is not available. Unable to retrieve data.", stacklevel=2
        )
        return pd.DataFrame(), {}

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
