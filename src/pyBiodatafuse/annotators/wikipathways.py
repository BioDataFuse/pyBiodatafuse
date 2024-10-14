#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Python file for queriying Wikipathways SPARQL endpoint ()."""

import datetime
import os
import time
import warnings
from string import Template
from typing import Any, Dict

import pandas as pd
from SPARQLWrapper import JSON, SPARQLWrapper
from SPARQLWrapper.SPARQLExceptions import SPARQLWrapperException
from tqdm import tqdm

from pyBiodatafuse.constants import (
    WIKIPATHWAYS,
    WIKIPATHWAYS_ENDPOINT,
    WIKIPATHWAYS_GENE_INPUT_ID,
    WIKIPATHWAYS_PATHWAYS_OUTPUT_DICT,
)
from pyBiodatafuse.utils import (
    check_columns_against_constants,
    collapse_data_sources,
    get_identifier_of_interest,
)


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
    # Check if the endpoint is available
    api_available = check_endpoint_wikipathways()

    if not api_available:
        warnings.warn(
            f"{WIKIPATHWAYS} SPARQL endpoint is not available. Unable to retrieve data.",
            stacklevel=2,
        )
        return pd.DataFrame(), {}

    data_df = get_identifier_of_interest(bridgedb_df, WIKIPATHWAYS_GENE_INPUT_ID)

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

    intermediate_df = pd.DataFrame()

    for gene_list_str in tqdm(query_gene_lists, desc="Querying WikiPathways"):
        sparql_query_template = Template(sparql_query)
        substit_dict = dict(gene_list=gene_list_str)
        sparql_query_template_sub = sparql_query_template.substitute(substit_dict)

        sparql.setQuery(sparql_query_template_sub)

        res = sparql.queryAndConvert()

        res = res["results"]["bindings"]

        df = pd.DataFrame(res)
        for col in df:
            df[col] = df[col].map(lambda x: x["value"], na_action="ignore")

        intermediate_df = pd.concat([intermediate_df, df], ignore_index=True)

    # Record the end time
    end_time = datetime.datetime.now()

    """Metdata details"""
    # Get the current date and time
    current_date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    # Calculate the time elapsed
    time_elapsed = str(end_time - start_time)

    # Add the datasource, query, query time, and the date to metadata
    wikipathways_metadata: Dict[str, Any] = {
        "datasource": WIKIPATHWAYS,
        "metadata": wikipathways_version,
        "query": {
            "size": len(gene_list),
            "input_type": WIKIPATHWAYS_GENE_INPUT_ID,
            "time": time_elapsed,
            "date": current_date,
            "url": WIKIPATHWAYS_ENDPOINT,
        },
    }

    if "gene_id" not in intermediate_df.columns:
        warnings.warn(
            f"There is no annotation for your input list in {WIKIPATHWAYS}.",
            stacklevel=2,
        )
        return pd.DataFrame(), wikipathways_metadata

    # Organize the annotation results as an array of dictionaries
    intermediate_df.rename(columns={"gene_id": "target"}, inplace=True)
    intermediate_df["pathway_gene_count"] = pd.to_numeric(intermediate_df["pathway_gene_count"])
    intermediate_df = intermediate_df.drop_duplicates()
    intermediate_df["pathway_id"] = intermediate_df["pathway_id"].apply(lambda x: f"WP:{x}")

    # Check if all keys in df match the keys in OUTPUT_DICT
    check_columns_against_constants(
        data_df=intermediate_df,
        output_dict=WIKIPATHWAYS_PATHWAYS_OUTPUT_DICT,
        check_values_in=["pathway_id"],
    )

    # Merge the two DataFrames on the target column
    merged_df = collapse_data_sources(
        data_df=data_df,
        source_namespace=WIKIPATHWAYS_GENE_INPUT_ID,
        target_df=intermediate_df,
        common_cols=["target"],
        target_specific_cols=list(WIKIPATHWAYS_PATHWAYS_OUTPUT_DICT.keys()),
        col_name=WIKIPATHWAYS,
    )

    """Update metadata"""
    # Calculate the number of new nodes
    num_new_nodes = intermediate_df["pathway_id"].nunique()
    # Calculate the number of new edges
    num_new_edges = intermediate_df.drop_duplicates(subset=["target", "pathway_id"]).shape[0]

    # Check the intermediate_df
    if num_new_edges != len(intermediate_df):
        warnings.warn(
            f"The intermediate_df in {WIKIPATHWAYS} annotatur should be checked, please create an issue on https://github.com/BioDataFuse/pyBiodatafuse/issues/.",
            stacklevel=2,
        )

    # Add the number of new nodes and edges to metadata
    wikipathways_metadata["query"]["number_of_added_nodes"] = num_new_nodes
    wikipathways_metadata["query"]["number_of_added_edges"] = num_new_edges

    return merged_df, wikipathways_metadata
