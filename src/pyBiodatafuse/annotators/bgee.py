#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Python file for queriying Bgee database (https://bgee.org)."""

import datetime
import os
import warnings
from string import Template
from typing import Any, Dict

import pandas as pd
from SPARQLWrapper import JSON, SPARQLWrapper
from SPARQLWrapper.SPARQLExceptions import SPARQLWrapperException
from tqdm import tqdm

import pyBiodatafuse.constants as Cons
from pyBiodatafuse.utils import (
    check_columns_against_constants,
    collapse_data_sources,
    get_identifier_of_interest,
    give_annotator_warning,
)


def check_sparql_endpoint_bgee() -> bool:
    """Check the availability of the Bgee SPARQL endpoint.

    :returns: True if the endpoint is available, False otherwise.
    """
    query_file_path = os.path.join(
        os.path.dirname(__file__), "queries", "bgee-get-last-modified.rq"
    )
    with open(query_file_path, "r", encoding="utf-8") as fin:
        sparql_query = fin.read()

    sparql = SPARQLWrapper(Cons.BGEE_ENDPOINT)
    sparql.setReturnFormat(JSON)

    sparql.setQuery(sparql_query)

    try:
        sparql.queryAndConvert()
        return True
    except SPARQLWrapperException:
        return False


def get_version_bgee() -> dict:
    """Get version of Bgee RDF data from its SPARQL endpoint.

    # not sure if a version per-se can be retrieved, but the endpoint supports
    # http://purl.org/dc/terms/modified
    :returns: a dictionary containing the last modified date information
    """
    with open(
        os.path.join(os.path.dirname(__file__), "queries", "bgee-get-last-modified.rq"),
        "r",
        encoding="utf-8",
    ) as fin:
        sparql_query = fin.read()

    sparql = SPARQLWrapper(Cons.BGEE_ENDPOINT)
    sparql.setReturnFormat(JSON)

    sparql.setQuery(sparql_query)
    res = sparql.queryAndConvert()
    bgee_version = {"source_version": res["results"]["bindings"][0]["date_modified"]["value"]}

    return bgee_version


def get_gene_expression(bridgedb_df: pd.DataFrame):
    """Query gene-tissue expression information from Bgee with SPARQL.

    :param bridgedb_df: BridgeDb output for creating the list of gene ids to query
    :returns: a DataFrame containing the Bgee output and dictionary of the Bgee metadata.
    """
    # Check if the Bgee SPARQL endpoint is available
    api_available = check_sparql_endpoint_bgee()

    if not api_available:
        warnings.warn(
            f"{Cons.BGEE} SPARQL endpoint is not available. Unable to retrieve data.", stacklevel=2
        )
        return pd.DataFrame(), {}

    # Extract the "target" values and join them into a single string separated by commas
    data_df = get_identifier_of_interest(bridgedb_df, Cons.BGEE_GENE_INPUT_ID)
    gene_list = data_df[Cons.TARGET_COL].tolist()
    gene_list = list(set(gene_list))

    query_gene_lists = []
    if len(gene_list) > 25:
        for i in range(0, len(gene_list), 25):
            tmp_list = gene_list[i : i + 25]
            query_gene_lists.append(" ".join(f'"{g}"' for g in tmp_list))

    else:
        query_gene_lists.append(" ".join(f'"{g}"' for g in gene_list))

    anatomical_entities_list = [
        anatomical_entity.strip()
        for anatomical_entity in Cons.ANATOMICAL_ENTITIES_LIST
        if anatomical_entity.strip() != ""
    ]

    query_file_path = os.path.join(
        os.path.dirname(__file__), "queries", "bgee-genes-tissues-expression-level.rq"
    )
    with open(query_file_path, "r", encoding="utf-8") as fin:
        sparql_query = fin.read()

    # Add version to metadata file
    bgee_version = get_version_bgee()

    # Record the start time
    start_time = datetime.datetime.now()

    sparql = SPARQLWrapper(Cons.BGEE_ENDPOINT)
    sparql.setReturnFormat(JSON)

    query_count = 0

    intermediate_df = pd.DataFrame()

    for gene_list_str in tqdm(query_gene_lists, desc="Querying Bgee"):
        query_count += 1

        gene_ids = gene_list_str.split(" ")

        for gene_id in gene_ids:
            sparql_query_template = Template(sparql_query)

            for anatomical_entity in anatomical_entities_list:
                # for the query text, need to put each name in between quotes
                substit_dict = dict(gene_list=gene_id, anat_entities_list=f'"{anatomical_entity}"')
                sparql_query_template_sub = sparql_query_template.substitute(substit_dict)

                sparql.setQuery(sparql_query_template_sub)
                res = sparql.queryAndConvert()

                df = pd.DataFrame(res["results"]["bindings"])

                for col in df:
                    df[col] = df[col].map(lambda x: x["value"], na_action="ignore")

                if df.empty:
                    continue

                df.drop_duplicates(subset=[Cons.ANATOMICAL_ID], inplace=True)
                intermediate_df = pd.concat([intermediate_df, df], ignore_index=True)

    # Record the end time
    end_time = datetime.datetime.now()

    """Metadata details"""
    # Get the current date and time
    current_date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    # Calculate the time elapsed
    time_elapsed = str(end_time - start_time)

    # Add the datasource, query, query time, and the date to metadata
    bgee_metadata: Dict[str, Any] = {
        "datasource": Cons.BGEE,
        "metadata": bgee_version,
        "query": {
            "size": len(gene_list),
            "input_type": Cons.BGEE_GENE_INPUT_ID,
            "time": time_elapsed,
            "date": current_date,
            "url": Cons.BGEE_ENDPOINT,
        },
    }

    if Cons.ANATOMICAL_ID not in intermediate_df:
        warnings.warn(
            f"There is no annotation for your input list in {Cons.BGEE}.",
            stacklevel=2,
        )
        return pd.DataFrame(), bgee_metadata

    # Organize the annotation results as an array of dictionaries
    intermediate_df.rename(columns={"ensembl_id": Cons.TARGET_COL}, inplace=True)
    intermediate_df[Cons.ANATOMICAL_ID] = intermediate_df[Cons.ANATOMICAL_ID].apply(
        lambda x: x.split("/")[-1]
    )
    intermediate_df[Cons.ANATOMICAL_ID] = intermediate_df[Cons.ANATOMICAL_ID].str.replace("_", ":")
    intermediate_df[Cons.CONFIDENCE_ID] = intermediate_df[Cons.CONFIDENCE_ID].apply(
        lambda x: x.split("/")[-1]
    )
    intermediate_df[Cons.CONFIDENCE_ID] = intermediate_df[Cons.CONFIDENCE_ID].str.replace("_", ":")
    intermediate_df[Cons.DEVELOPMENTAL_ID] = intermediate_df[Cons.DEVELOPMENTAL_ID].apply(
        lambda x: x.split("/")[-1]
    )
    intermediate_df[Cons.DEVELOPMENTAL_ID] = intermediate_df[Cons.DEVELOPMENTAL_ID].str.replace(
        "_", ":"
    )
    intermediate_df[Cons.EXPRESSION_LEVEL] = pd.to_numeric(intermediate_df[Cons.EXPRESSION_LEVEL])

    # Check if all keys in df match the keys in OUTPUT_DICT
    check_columns_against_constants(
        data_df=intermediate_df,
        output_dict=Cons.BGEE_GENE_EXPRESSION_OUTPUT_DICT,
        check_values_in=Cons.BGEE_VALUE_CHECK_LIST,
    )

    # Merge the two DataFrames on the target column
    merged_df = collapse_data_sources(
        data_df=data_df,
        source_namespace=Cons.BGEE_GENE_INPUT_ID,
        target_df=intermediate_df,
        common_cols=[Cons.TARGET_COL],
        target_specific_cols=list(Cons.BGEE_GENE_EXPRESSION_OUTPUT_DICT.keys()),
        col_name=Cons.BGEE_GENE_EXPRESSION_LEVELS_COL,
    )

    """Update metadata"""
    # Calculate the number of new nodes
    num_new_nodes = intermediate_df[Cons.ANATOMICAL_ID].nunique()
    # Calculate the number of new edges
    num_new_edges = intermediate_df.drop_duplicates(
        subset=[Cons.ANATOMICAL_ID, Cons.TARGET_COL]
    ).shape[0]

    # Check the intermediate_df
    if num_new_edges != len(intermediate_df):
        give_annotator_warning(Cons.BGEE)

    # Add the number of new nodes and edges to metadata
    bgee_metadata[Cons.QUERY][Cons.NUM_NODES] = num_new_nodes
    bgee_metadata[Cons.QUERY][Cons.NUM_EDGES] = num_new_edges

    return merged_df, bgee_metadata
