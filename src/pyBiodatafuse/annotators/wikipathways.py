#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Python file for queriying Wikipathways SPARQL endpoint ()."""

import datetime
import logging
import os
import warnings
from string import Template

import pandas as pd
from SPARQLWrapper import JSON, SPARQLWrapper
from SPARQLWrapper.SPARQLExceptions import SPARQLWrapperException

from pyBiodatafuse.constants import (
    WIKIPATHWAYS,
    WIKIPATHWAYS_ENDPOINT,
    WIKIPATHWAYS_INPUT_ID,
    WIKIPATHWAYS_OUTPUT_DICT,
)
from pyBiodatafuse.utils import (
    check_columns_against_constants,
    collapse_data_sources,
    get_identifier_of_interest,
)

logger = logging.getLogger("wikipathways")


def check_endpoint_wikipathways() -> bool:
    """Check the availability of the WikiPathways SPARQL endpoint.

    :returns: True if the endpoint is available, False otherwise.
    """
    with open(os.path.dirname(__file__) + "/queries/wikipathways-metadata.rq", "r") as fin:
        sparql_query = fin.read()

    sparql = SPARQLWrapper(WIKIPATHWAYS_ENDPOINT)
    sparql.setReturnFormat(JSON)

    sparql.setQuery(sparql_query)

    try:
        sparql.queryAndConvert()
        return True
    except SPARQLWrapperException:
        return False


def get_version_wikipathways() -> dict:
    """Get version of WikiPathways.

    :returns: a dictionary containing the version information
    """
    with open(os.path.dirname(__file__) + "/queries/wikipathways-metadata.rq", "r") as fin:
        sparql_query = fin.read()

    sparql = SPARQLWrapper(WIKIPATHWAYS_ENDPOINT)
    sparql.setReturnFormat(JSON)

    sparql.setQuery(sparql_query)
    res = sparql.queryAndConvert()

    wikipathways_version = {"source_version": res["results"]["bindings"][0]["title"]["value"]}

    return wikipathways_version


def get_gene_wikipathways(bridgedb_df: pd.DataFrame):
    """Query WikiPathways for pathways associated with genes.

    :param bridgedb_df: BridgeDb output for creating the list of gene ids to query
    :returns: a DataFrame containing the WikiPathways output and dictionary of the WikiPathways metadata.
    """
    # Check if the DisGeNET API is available
    api_available = check_endpoint_wikipathways()

    if not api_available:
        warnings.warn(
            f"{WIKIPATHWAYS} SPARQL endpoint is not available. Unable to retrieve data.",
            stacklevel=2,
        )
        return pd.DataFrame(), {}

    data_df = get_identifier_of_interest(bridgedb_df, WIKIPATHWAYS_INPUT_ID)

    wikipathways_version = get_version_wikipathways()
    gene_list = data_df["target"].tolist()
    gene_list = list(set(gene_list))

    query_gene_lists = []

    if len(gene_list) > 25:
        for i in range(0, len(gene_list), 25):
            tmp_list = gene_list[i : i + 25]
            query_gene_lists.append(" ".join(f'"{g}"' for g in tmp_list))

    else:
        query_gene_lists.append(" ".join(f'"{g}"' for g in gene_list))

    with open(os.path.dirname(__file__) + "/queries/wikipathways-genes-pathways.rq", "r") as fin:
        sparql_query = fin.read()

    # Record the start time
    start_time = datetime.datetime.now()

    sparql = SPARQLWrapper(WIKIPATHWAYS_ENDPOINT)
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

        intermediate_df = pd.concat([intermediate_df, df], ignore_index=True)

    # Record the end time
    end_time = datetime.datetime.now()

    if "gene_id" not in intermediate_df.columns:
        return pd.DataFrame(), {"datasource": WIKIPATHWAYS, "metadata": wikipathways_version}

    # Organize the annotation results as an array of dictionaries
    intermediate_df.rename(columns={"gene_id": "target"}, inplace=True)
    intermediate_df["pathway_gene_count"] = pd.to_numeric(intermediate_df["pathway_gene_count"])
    intermediate_df = intermediate_df.drop_duplicates()

    # Check if all keys in df match the keys in OUTPUT_DICT
    check_columns_against_constants(
        data_df=intermediate_df,
        output_dict=WIKIPATHWAYS_OUTPUT_DICT,
        check_values_in=["pathway_id"],
    )

    # Merge the two DataFrames on the target column
    merged_df = collapse_data_sources(
        data_df=data_df,
        source_namespace=WIKIPATHWAYS_INPUT_ID,
        target_df=intermediate_df,
        common_cols=["target"],
        target_specific_cols=list(WIKIPATHWAYS_OUTPUT_DICT.keys()),
        col_name=WIKIPATHWAYS,
    )

    # Metdata details
    # Get the current date and time
    current_date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    # Calculate the time elapsed
    time_elapsed = str(end_time - start_time)
    # Calculate the number of disease nodes
    num_new_nodes = intermediate_df["pathway_id"].nunique()
    # Calculate the number of edges between gene and disease nodes
    num_edges = len(intermediate_df)

    # Add the datasource, query, query time, and the date to metadata
    wikipathways_metadata = {
        "datasource": WIKIPATHWAYS,
        "metadata": wikipathways_version,
        "query": {
            "size": len(gene_list),
            "number_of_added_nodes": num_new_nodes,
            "number_of_added_edges": num_edges,
            "input_type": WIKIPATHWAYS_INPUT_ID,
            "time": time_elapsed,
            "date": current_date,
            "url": WIKIPATHWAYS_ENDPOINT,
        },
    }

    return merged_df, wikipathways_metadata
