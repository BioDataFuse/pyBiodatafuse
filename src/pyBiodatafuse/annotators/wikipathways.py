#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Python file for queriying Wikipathways SPARQL endpoint ()."""

import datetime
import os
import re
import time
import warnings
from string import Template
from typing import Any, Dict, Mapping, Type

import pandas as pd
from SPARQLWrapper import JSON, SPARQLWrapper
from SPARQLWrapper.SPARQLExceptions import SPARQLWrapperException
from tqdm import tqdm

from pyBiodatafuse.constants import (
    WIKIPATHWAYS,
    WIKIPATHWAYS_ENDPOINT,
    WIKIPATHWAYS_GENE_INPUT_ID,
    WIKIPATHWAYS_MOLECULAR_COL,
    WIKIPATHWAYS_MOLECULAR_GENE_OUTPUT_DICT,
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


def get_gene_wikipathways(
    bridgedb_df: pd.DataFrame, query_interactions: bool = False
) -> pd.DataFrame:
    """Query WikiPathways for pathways associated with genes.

    :param bridgedb_df: BridgeDb output for creating the list of gene ids to query
    :param query_interactions: Set whether to retrieve gene part_of pathways relationships (False) or all molecular interactions (True).
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
        output_dict = WIKIPATHWAYS_MOLECULAR_GENE_OUTPUT_DICT
        col_name = WIKIPATHWAYS_MOLECULAR_COL
    else:
        file = os.path.join("queries", "wikipathways-genes-pathways.rq")
        output_dict = WIKIPATHWAYS_PATHWAYS_OUTPUT_DICT
        col_name = WIKIPATHWAYS
    with open(os.path.join(os.path.dirname(__file__), file), "r", encoding="utf-8") as fin:
        sparql_query = fin.read()

    # Record the start time
    start_time = datetime.datetime.now()

    sparql = SPARQLWrapper(WIKIPATHWAYS_ENDPOINT)
    sparql.setReturnFormat(JSON)

    intermediate_df = pd.DataFrame()

    for gene_list_str in tqdm(query_gene_lists, desc="Querying WikiPathways"):
        sparql_query_template = Template(sparql_query)
        substit_dict = dict()
        if query_interactions:
            substit_dict = dict(
                gene_list=gene_list_str, organism_name='"Homo sapiens"'
            )  # TODO allow setting organism
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
    # Fix identifiers
    intermediate_df["gene_id"] = (
        intermediate_df["gene_id"].str.removeprefix("https://identifiers.org/ncbigene/").fillna("")
    )
    if query_interactions:
        intermediate_df["targetGene"] = (
            intermediate_df["targetGene"]
            .str.removeprefix("https://identifiers.org/ncbigene/")
            .fillna("")
        )
        intermediate_df["targetMetabolite"] = (
            intermediate_df["targetMetabolite"]
            .str.removeprefix("http://rdf.ncbi.nlm.nih.gov/pubchem/compound/")
            .fillna("")
        )
        intermediate_df["targetProtein"] = (
            intermediate_df["targetProtein"]
            .str.removeprefix("https://identifiers.org/uniprot/")
            .fillna("")
        )
        intermediate_df["mimtype"] = (
            intermediate_df["mimtype"]
            .str.removeprefix("http://vocabularies.wikipathways.org/wp#")
            .fillna("")
        )
    else:
        intermediate_df["pathway_gene_count"] = pd.to_numeric(
            intermediate_df["pathway_gene_count"], errors="coerce"
        )
    # Organize the annotation results as an array of dictionaries
    intermediate_df.rename(columns={"gene_id": "target"}, inplace=True)
    if not query_interactions:
        intermediate_df["pathway_gene_count"] = pd.to_numeric(intermediate_df["pathway_gene_count"])
    if query_interactions:
        intermediate_df["pathway_id"] = (
            intermediate_df["pathway_id"]
            .str.removeprefix("https://identifiers.org/wikipathways/")
            .str.replace(r"(WP\d+)_.*", r"\1", regex=True)
            .str.replace("WP", "WP:", regex=False)
            .str.strip()
        )
    else:
        intermediate_df["pathway_id"] = intermediate_df["pathway_id"].apply(
            lambda x: x.replace("WP", "WP:WP") if "WP" in x else x
        )

    # Check if all keys in df match the keys in OUTPUT_DICT
    check_columns_against_constants(
        data_df=intermediate_df,
        output_dict=output_dict,
        check_values_in=["pathway_id"],
    )

    # Merge the two DataFrames on the target column
    merged_df = collapse_data_sources(
        data_df=data_df,
        source_namespace=WIKIPATHWAYS_GENE_INPUT_ID,
        target_df=intermediate_df,
        common_cols=["target"],
        target_specific_cols=list(output_dict.keys()),
        col_name=col_name,
    )

    """Update metadata"""
    # Calculate the number of new nodes
    num_new_nodes = intermediate_df["pathway_id"].nunique()
    # Calculate the number of new edges
    num_new_edges = intermediate_df.drop_duplicates(subset=["target", "pathway_id"]).shape[0]
    if not query_interactions:
        # Check the intermediate_df
        if num_new_edges != len(intermediate_df):
            warnings.warn(
                f"The intermediate_df in {WIKIPATHWAYS} annotator should be checked, please create an issue on https://github.com/BioDataFuse/pyBiodatafuse/issues/.",
                stacklevel=2,
            )

    # Add the number of new nodes and edges to metadata
    wikipathways_metadata["query"]["number_of_added_nodes"] = num_new_nodes
    wikipathways_metadata["query"]["number_of_added_edges"] = num_new_edges

    return merged_df, wikipathways_metadata
