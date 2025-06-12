#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Python file for queriying Wikipathways SPARQL endpoint ()."""

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


def check_endpoint_wikipathways() -> bool:
    """Check the availability of the WikiPathways SPARQL endpoint.

    :returns: True if the endpoint is available, False otherwise.
    """
    with open(os.path.dirname(__file__) + "/queries/wikipathways-metadata.rq", "r") as fin:
        sparql_query = fin.read()

    sparql = SPARQLWrapper(Cons.WIKIPATHWAYS_ENDPOINT)
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

    sparql = SPARQLWrapper(Cons.WIKIPATHWAYS_ENDPOINT)
    sparql.setReturnFormat(JSON)

    sparql.setQuery(sparql_query)
    res = sparql.queryAndConvert()

    wikipathways_version = {"source_version": res["results"]["bindings"][0]["title"]["value"]}

    return wikipathways_version


def get_gene_wikipathways(
    bridgedb_df: pd.DataFrame, query_interactions: bool = False, organism: str = "Homo sapiens"
) -> pd.DataFrame:
    """Query WikiPathways for pathways associated with genes.

    :param bridgedb_df: BridgeDb output for creating the list of gene ids to query
    :param query_interactions: Set whether to retrieve gene part_of pathways relationships (False) or all molecular interactions (True).
    :param organism: The organism to query. Default is "Homo sapiens".
    :returns: a DataFrame containing the WikiPathways output and dictionary of the WikiPathways metadata.
    """
    # Check if the endpoint is available
    api_available = check_endpoint_wikipathways()

    if not api_available:
        warnings.warn(
            f"{Cons.WIKIPATHWAYS} SPARQL endpoint is not available. Unable to retrieve data.",
            stacklevel=2,
        )
        return pd.DataFrame(), {}

    data_df = get_identifier_of_interest(bridgedb_df, Cons.WIKIPATHWAYS_GENE_INPUT_ID)

    wikipathways_version = get_version_wikipathways()
    gene_list = data_df[Cons.TARGET_COL].tolist()
    gene_list = list(set(gene_list))

    query_gene_lists = []
    output_dict: Dict[Any, Any]
    if len(gene_list) > 25:
        for i in range(0, len(gene_list), 25):
            tmp_list = gene_list[i : i + 25]
            query_gene_lists.append(" ".join(f'"{g}"' for g in tmp_list))

    else:
        if query_interactions:
            gene_list = [f"<https://identifiers.org/ncbigene/{g}>" for g in gene_list]
        query_gene_lists.append(" ".join(g for g in gene_list))
    col_name = ""

    if query_interactions:
        file = os.path.join("queries", "wikipathways-mims.rq")
        output_dict = Cons.WIKIPATHWAYS_MOLECULAR_GENE_OUTPUT_DICT
        col_name = Cons.WIKIPATHWAYS_MOLECULAR_COL
    else:
        file = os.path.join("queries", "wikipathways-genes-pathways.rq")
        output_dict = Cons.WIKIPATHWAYS_PATHWAYS_OUTPUT_DICT
        col_name = Cons.WIKIPATHWAYS_PATHWAY_COL

    with open(os.path.join(os.path.dirname(__file__), file), "r", encoding="utf-8") as fin:
        sparql_query = fin.read()

    # Record the start time
    start_time = datetime.datetime.now()

    sparql = SPARQLWrapper(Cons.WIKIPATHWAYS_ENDPOINT)
    sparql.setReturnFormat(JSON)

    intermediate_df = pd.DataFrame()

    for gene_list_str in tqdm(query_gene_lists, desc="Querying WikiPathways"):
        sparql_query_template = Template(sparql_query)
        substit_dict = dict()
        if query_interactions:
            substit_dict = dict(gene_list=gene_list_str, organism_name=f'"{organism}"')

        if not query_interactions:
            substit_dict = dict(gene_list=gene_list_str)

        sparql_query_template_sub = sparql_query_template.substitute(substit_dict)
        sparql.setQuery(sparql_query_template_sub)

        result = sparql.queryAndConvert()

        res = result["results"]["bindings"]  # get data
        df = pd.DataFrame(res)
        for col in df:
            df[col] = df[col].map(lambda x: x["value"], na_action="ignore")

        # Retrieve the expected columns from the SPARQL query results' "vars" (deal with empty optionals)
        expected_columns = result["head"]["vars"]

        # Ensure all expected columns are present in intermediate_df
        for col in expected_columns:
            if col not in intermediate_df.columns:
                intermediate_df[col] = None  # Add missing columns with default value None

        # Concatenate the new results into the intermediate DataFrame
        intermediate_df = pd.concat([intermediate_df, df], ignore_index=True)
    # Record the end time
    end_time = datetime.datetime.now()

    """Metadata details"""
    # Get the current date and time
    current_date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    # Calculate the time elapsed
    time_elapsed = str(end_time - start_time)

    # Add the datasource, query, query time, and the date to metadata
    wikipathways_metadata: Dict[str, Any] = {
        "datasource": Cons.WIKIPATHWAYS,
        "metadata": wikipathways_version,
        "query": {
            "size": len(gene_list),
            "input_type": Cons.WIKIPATHWAYS_GENE_INPUT_ID,
            "time": time_elapsed,
            "date": current_date,
            "url": Cons.WIKIPATHWAYS_ENDPOINT,
        },
    }

    if Cons.WIKIPATHWAYS_GENE_ID not in intermediate_df.columns:
        warnings.warn(
            f"There is no annotation for your input list in {Cons.WIKIPATHWAYS}.",
            stacklevel=2,
        )
        return pd.DataFrame(), wikipathways_metadata

    # Fix identifiers
    for col in Cons.WIKIPATHWAY_ID_CLEANER_DICT:
        if col not in intermediate_df.columns:
            continue

        intermediate_df[col] = (
            intermediate_df[col].str.removeprefix(Cons.WIKIPATHWAY_ID_CLEANER_DICT[col]).fillna("")
        )

    # Organize the annotation results as an array of dictionaries
    intermediate_df.rename(
        columns={
            Cons.WIKIPATHWAYS_GENE_ID: Cons.TARGET_COL,
            "pathway_gene_count": Cons.PATHWAY_GENE_COUNTS,
        },
        inplace=True,
    )

    if Cons.PATHWAY_GENE_COUNTS in intermediate_df.columns:
        intermediate_df[Cons.PATHWAY_GENE_COUNTS] = pd.to_numeric(
            intermediate_df[Cons.PATHWAY_GENE_COUNTS], errors="coerce"
        )

    # Add namespaces
    for col, namespace in Cons.WIKIPATHWAY_NAMESPACE_DICT.items():
        if col not in intermediate_df.columns:
            continue
        intermediate_df[col] = intermediate_df[col].apply(
            lambda x, ns=namespace: f"{ns}:{x}" if x else x
        )

    if query_interactions:
        intermediate_df[Cons.PATHWAY_ID] = (
            intermediate_df[Cons.PATHWAY_ID]
            .str.removeprefix("https://identifiers.org/wikipathways/")
            .str.replace(r"(WP\d+)_.*", r"\1", regex=True)
            .str.replace("WP", f"{Cons.WIKIPATHWAY}:", regex=False)
            .str.strip()
        )
    else:
        intermediate_df[Cons.PATHWAY_ID] = intermediate_df[Cons.PATHWAY_ID].apply(
            lambda x: x.replace("WP", f"{Cons.WIKIPATHWAY}:WP") if "WP" in x else x
        )

    # Check if all keys in df match the keys in OUTPUT_DICT
    check_columns_against_constants(
        data_df=intermediate_df,
        output_dict=output_dict,
        check_values_in=[Cons.WIKIPATHWAY],
    )

    # Merge the two DataFrames on the target column
    merged_df = collapse_data_sources(
        data_df=data_df,
        source_namespace=Cons.WIKIPATHWAYS_GENE_INPUT_ID,
        target_df=intermediate_df,
        common_cols=[Cons.TARGET_COL],
        target_specific_cols=list(output_dict.keys()),
        col_name=col_name,
    )

    """Update metadata"""
    # Calculate the number of new nodes
    num_new_nodes = intermediate_df[Cons.PATHWAY_ID].nunique()
    # Calculate the number of new edges
    num_new_edges = intermediate_df.drop_duplicates(
        subset=[Cons.TARGET_COL, Cons.PATHWAY_ID]
    ).shape[0]
    if not query_interactions:
        # Check the intermediate_df
        if num_new_edges != len(intermediate_df):
            give_annotator_warning(Cons.WIKIPATHWAYS)

    # Add the number of new nodes and edges to metadata
    wikipathways_metadata[Cons.QUERY][Cons.NUM_NODES] = num_new_nodes
    wikipathways_metadata[Cons.QUERY][Cons.NUM_EDGES] = num_new_edges

    return merged_df, wikipathways_metadata
