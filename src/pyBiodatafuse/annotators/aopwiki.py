#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Python file for querying the AOP Wiki RDF SPARQL endpoint.

This module provides functionality to query the AOP Wiki RDF SPARQL endpoint for
Adverse Outcome Pathways (AOPs) associated with genes and compounds.
"""

import datetime
import logging
import os
import warnings
from string import Template
from typing import Any, Dict, Tuple

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

# Pre-requisite:
QUERY_LIMIT = 25
QUERY_COMPOUND = os.path.join(os.path.dirname(__file__), "queries", "aopwiki-compound.rq")
QUERY_GENE = os.path.join(os.path.dirname(__file__), "queries", "aopwiki-gene.rq")

logger = logging.getLogger(__name__)


def read_sparql_file(file_path: str) -> str:
    """Read a SPARQL query file.

    :param file_path: the path to the SPARQL query file
    :returns: the content of the SPARQL query file
    """
    with open(file_path, "r") as fin:
        sparql_query = fin.read()

    return sparql_query


def check_endpoint_aopwiki() -> bool:
    """Check the availability of the AOP Wiki SPARQL endpoint.

    :returns: True if the endpoint is available, False otherwise.
    """
    try:
        sparql = SPARQLWrapper(Cons.AOPWIKI_ENDPOINT)
        sparql.setReturnFormat(JSON)
        sparql.setQuery("SELECT * WHERE { ?s ?p ?o } LIMIT 1")
        sparql.queryAndConvert()
        return True
    except SPARQLWrapperException:
        return False


def get_aops_gene(bridgedb_df: pd.DataFrame) -> Tuple[pd.DataFrame, dict]:
    """Query for AOPs associated with genes from AOP Wiki RDF.

    :param bridgedb_df: BridgeDb output for creating the list of gene ids to query
    :returns: a DataFrame containing the AOP Wiki RDF output and dictionary of the AOP Wiki RDF metadata
    """
    # Check if the endpoint is available
    if not check_endpoint_aopwiki():
        warnings.warn(
            f"{Cons.AOPWIKIRDF} SPARQL endpoint is not available. Unable to retrieve data.",
            stacklevel=2,
        )
        return pd.DataFrame(), {}

    # Record the start time
    start_time = datetime.datetime.now()

    # Step 1: Identifier mapping and harmonization
    data_df = get_identifier_of_interest(bridgedb_df, Cons.AOPWIKI_GENE_INPUT_ID)
    gene_list = data_df[Cons.TARGET_COL].unique().tolist()

    # Step 2: Prepare target list and batch queries
    query_batches = [
        " ".join(f'"{target}"' for target in gene_list[i : i + QUERY_LIMIT])
        for i in range(0, len(gene_list), QUERY_LIMIT)
    ]

    # Step 3: Run SPARQL queries
    sparql = SPARQLWrapper(Cons.AOPWIKI_ENDPOINT)
    sparql.setReturnFormat(JSON)

    intermediate_df = pd.DataFrame()

    for batch in tqdm(query_batches, desc=f"Querying {Cons.AOPWIKIRDF} for genes"):
        # Prepare the substitution dictionary
        substit_dict = {
            "genes": str(['"' + target.replace('"', "") + '"' for target in batch.split(" ")])
            .replace("[", "")
            .replace("]", "")
            .replace("'", "")
            .replace(",", "")
        }

        # Load and substitute the query template
        with open(QUERY_GENE, "r") as f:
            query = Template(f.read()).substitute(substit_dict)

        # Execute the query and process results
        sparql.setQuery(query)
        res = sparql.queryAndConvert()
        res_df = pd.DataFrame(
            [
                {k: (v["value"] if "value" in v else "") for k, v in item.items()}
                for item in res["results"]["bindings"]
            ]
        )

        # Retrieve the expected columns from the SPARQL query results' "vars"
        expected_columns = res["head"]["vars"]

        # Ensure all expected columns are present in intermediate_df
        for col in expected_columns:
            if col not in intermediate_df.columns:
                intermediate_df[col] = None

        # Concatenate the new results into the intermediate DataFrame
        intermediate_df = pd.concat([intermediate_df, res_df], ignore_index=True)

    # Record the end time
    end_time = datetime.datetime.now()

    # Step 4: Check if the query returned any results
    if intermediate_df.empty:
        warnings.warn(
            f"There is no annotation for your input list in {Cons.AOPWIKIRDF}.",
            stacklevel=2,
        )
        return pd.DataFrame(), {}

    # Step 5: Clean and process the results
    input_col = Cons.AOPWIKI_GENE_INPUT_ID
    output_dict = Cons.AOPWIKI_GENE_OUTPUT_DICT

    for key in output_dict.keys():
        if key in intermediate_df.columns:
            intermediate_df[key] = intermediate_df[key].apply(
                lambda x: x.split("/")[-1] if isinstance(x, str) and "http" in x else x
            )

    intermediate_df.rename(columns={input_col: Cons.TARGET_COL}, inplace=True)
    intermediate_df = intermediate_df.drop_duplicates()

    # Step 6: Generate metadata
    current_date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    time_elapsed = str(end_time - start_time)

    metadata_dict: Dict[str, Any] = {
        "datasource": Cons.AOPWIKIRDF,
        "query": {
            "size": len(gene_list),
            "input_type": Cons.AOPWIKI_GENE_INPUT_ID,
            "time": time_elapsed,
            "date": current_date,
            "url": Cons.AOPWIKI_ENDPOINT,
        },
    }

    # Check if all keys in df match the keys in OUTPUT_DICT
    check_columns_against_constants(
        data_df=intermediate_df,
        output_dict=output_dict,
        check_values_in=[],  # no prefix for AOP Wiki RDF
    )

    # Step 7: Integrate into main dataframe
    merged_df = collapse_data_sources(
        data_df=data_df,
        source_namespace=Cons.AOPWIKI_GENE_INPUT_ID,
        target_df=intermediate_df,
        common_cols=[Cons.TARGET_COL],
        target_specific_cols=list(output_dict.keys()),
        col_name=Cons.AOPWIKI_GENE_COL,
    )

    # Calculate the number of new nodes and edges
    num_new_nodes = intermediate_df[Cons.TARGET_COL].nunique()
    num_new_edges = intermediate_df.drop_duplicates(subset=[Cons.TARGET_COL]).shape[0]

    # Check the intermediate_df
    if num_new_edges != len(intermediate_df):
        give_annotator_warning(Cons.AOPWIKIRDF)

    # Add the number of new nodes and edges to metadata
    metadata_dict[Cons.QUERY][Cons.NUM_NODES] = num_new_nodes
    metadata_dict[Cons.QUERY][Cons.NUM_EDGES] = num_new_edges

    return merged_df, metadata_dict


def get_aops_compound(bridgedb_df: pd.DataFrame) -> Tuple[pd.DataFrame, dict]:
    """Query for AOPs associated with compounds from AOP Wiki RDF.

    :param bridgedb_df: BridgeDb output for creating the list of compound ids to query
    :returns: a DataFrame containing the AOP Wiki RDF output and dictionary of the AOP Wiki RDF metadata
    """
    # Check if the endpoint is available
    if not check_endpoint_aopwiki():
        warnings.warn(
            f"{Cons.AOPWIKIRDF} SPARQL endpoint is not available. Unable to retrieve data.",
            stacklevel=2,
        )
        return pd.DataFrame(), {}

    # Record the start time
    start_time = datetime.datetime.now()

    # Step 1: Identifier mapping and harmonization
    data_df = get_identifier_of_interest(bridgedb_df, Cons.AOPWIKI_COMPOUND_INPUT_ID)
    compound_list = data_df[Cons.TARGET_COL].unique().tolist()

    # Step 2: Prepare target list and batch queries
    query_batches = [
        " ".join(f'"{target}"' for target in compound_list[i : i + QUERY_LIMIT])
        for i in range(0, len(compound_list), QUERY_LIMIT)
    ]

    # Step 3: Run SPARQL queries
    sparql = SPARQLWrapper(Cons.AOPWIKI_ENDPOINT)
    sparql.setReturnFormat(JSON)

    intermediate_df = pd.DataFrame()

    for batch in tqdm(query_batches, desc=f"Querying {Cons.AOPWIKIRDF} for compounds"):
        # Prepare the substitution dictionary
        substit_dict = {
            "compounds": str(
                [
                    "<https://identifiers.org/pubchem.compound/" + target.replace('"', "") + ">"
                    for target in batch.split(" ")
                ]
            )
            .replace("[", "")
            .replace("]", "")
            .replace("'", "")
            .replace(",", "")
        }

        # Load and substitute the query template
        with open(QUERY_COMPOUND, "r") as f:
            query = Template(f.read()).substitute(substit_dict)

        # Execute the query and process results
        sparql.setQuery(query)
        res = sparql.queryAndConvert()
        res_df = pd.DataFrame(
            [
                {k: (v["value"] if "value" in v else "") for k, v in item.items()}
                for item in res["results"]["bindings"]
            ]
        )

        # Retrieve the expected columns from the SPARQL query results' "vars"
        expected_columns = res["head"]["vars"]

        # Ensure all expected columns are present in intermediate_df
        for col in expected_columns:
            if col not in intermediate_df.columns:
                intermediate_df[col] = None

        # Concatenate the new results into the intermediate DataFrame
        intermediate_df = pd.concat([intermediate_df, res_df], ignore_index=True)

    # Record the end time
    end_time = datetime.datetime.now()

    # Step 4: Check if the query returned any results
    if intermediate_df.empty:
        warnings.warn(
            f"There is no annotation for your input list in {Cons.AOPWIKI_COMPOUND_COL}.",
            stacklevel=2,
        )
        return pd.DataFrame(), {}

    # Step 5: Clean and process the results
    input_col = "pubchem_compound"
    output_dict = Cons.AOPWIKI_COMPOUND_OUTPUT_DICT

    # Clean URLs in output columns
    for key in output_dict.keys():
        if key in intermediate_df.columns:
            intermediate_df[key] = intermediate_df[key].apply(
                lambda x: x.split("/")[-1] if isinstance(x, str) and "http" in x else x
            )

    # Clean the compound ID column
    if input_col in intermediate_df.columns:
        intermediate_df[input_col] = intermediate_df[input_col].apply(lambda x: x.split("/")[-1])

    intermediate_df.rename(columns={input_col: Cons.TARGET_COL}, inplace=True)
    intermediate_df = intermediate_df.drop_duplicates()

    # Step 6: Generate metadata
    current_date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    time_elapsed = str(end_time - start_time)

    metadata_dict: Dict[str, Any] = {
        "datasource": Cons.AOPWIKIRDF,
        "query": {
            "size": len(compound_list),
            "input_type": Cons.AOPWIKI_COMPOUND_INPUT_ID,
            "time": time_elapsed,
            "date": current_date,
            "url": Cons.AOPWIKI_ENDPOINT,
        },
    }

    # Check if all keys in df match the keys in OUTPUT_DICT
    check_columns_against_constants(
        data_df=intermediate_df,
        output_dict=output_dict,
        check_values_in=[],  # no prefix for AOP Wiki RDF,
    )

    # Step 7: Integrate into main dataframe
    merged_df = collapse_data_sources(
        data_df=data_df,
        source_namespace=Cons.AOPWIKI_COMPOUND_INPUT_ID,
        target_df=intermediate_df,
        common_cols=[Cons.TARGET_COL],
        target_specific_cols=list(output_dict.keys()),
        col_name=Cons.AOPWIKI_COMPOUND_COL,
    )

    # Calculate the number of new nodes and edges
    num_new_nodes = intermediate_df[Cons.TARGET_COL].nunique()
    num_new_edges = intermediate_df.drop_duplicates(subset=[Cons.TARGET_COL]).shape[0]

    # Check the intermediate_df
    if num_new_edges != len(intermediate_df):
        give_annotator_warning(Cons.AOPWIKI_COMPOUND_COL)

    # Add the number of new nodes and edges to metadata
    metadata_dict[Cons.QUERY][Cons.NUM_NODES] = num_new_nodes
    metadata_dict[Cons.QUERY][Cons.NUM_EDGES] = num_new_edges

    return merged_df, metadata_dict


def get_aops(
    bridgedb_df: pd.DataFrame,
) -> Tuple[pd.DataFrame, dict]:
    """Query for AOPs associated with genes or compounds.

    :param bridgedb_df: BridgeDb output for creating the list of gene/compound ids to query
    :raises ValueError: if the input identifiers are not recognized or if they are not admitted gene or compound identifiers
    :returns: a DataFrame containing the AOP Wiki RDF output and dictionary of the AOP Wiki RDF metadata
    """
    # Find the matching type based on which input_identifier we find
    if (
        Cons.AOPWIKI_GENE_INPUT_ID in bridgedb_df["identifier.source"].values
        or Cons.AOPWIKI_GENE_INPUT_ID in bridgedb_df["target.source"].values
    ):
        return get_aops_gene(bridgedb_df)
    elif (
        Cons.AOPWIKI_COMPOUND_INPUT_ID in bridgedb_df["identifier.source"].values
        or Cons.AOPWIKI_COMPOUND_INPUT_ID in bridgedb_df["target.source"].values
    ):
        return get_aops_compound(bridgedb_df)
    else:
        raise ValueError(
            f"Input identifiers must be either '{Cons.AOPWIKI_GENE_INPUT_ID}' or '{Cons.AOPWIKI_COMPOUND_INPUT_ID}'"
        )
