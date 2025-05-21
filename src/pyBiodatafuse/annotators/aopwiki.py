#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Python file for queriying the AOP Wiki RDF SPARQL endpoint ()."""

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

from pyBiodatafuse.constants import (
    AOPWIKI_COMPOUND_COL,
    AOPWIKI_COMPOUND_INPUT_ID,
    AOPWIKI_COMPOUND_OUTPUT_DICT,
    AOPWIKI_ENDPOINT,
    AOPWIKI_GENE_INPUT_ID,
    AOPWIKI_GENE_OUTPUT_DICT,
    AOPWIKIRDF,
    OPENTARGETS_GO_COL,
)
from pyBiodatafuse.utils import (
    check_columns_against_constants,
    collapse_data_sources,
    get_identifier_of_interest,
)

# Pre-requisite:
# VERSION_QUERY_FILE = os.path.join(os.path.dirname(__file__), "queries", "aopwiki-metadata.rq")
DATABASE_SPARQL_DICT = {"aopwiki": AOPWIKI_ENDPOINT}
DATABASE_QUERY_IDENTIFER_GENE = AOPWIKI_GENE_INPUT_ID
DATABASE_QUERY_IDENTIFER_COMPOUND = AOPWIKI_COMPOUND_INPUT_ID
QUERY_LIMIT = 25
QUERY_COMPOUND = os.path.join(os.path.dirname(__file__), "queries", "aopwiki-compound.rq")
QUERY_GENE = os.path.join(os.path.dirname(__file__), "queries", "aopwiki-gene.rq")
# QUERY_PROCESS = os.path.join(os.path.dirname(__file__), "queries", "aopwiki-get-by-biological-process.rq.rq")
# PROCESS_INPUT_COL = OPENTARGETS_GO_COL

# Unique inputs and outputs:
INPUT_OPTIONS = [
    "gene",
    "compound",
]  # "biological_process"]


def read_sparql_file(file_path: str) -> str:
    """Read a SPARQL query file.

    :param file_path: the path to the SPARQL query file
    :returns: the content of the SPARQL query file
    """
    with open(file_path, "r") as fin:
        sparql_query = fin.read()

    return sparql_query


# def check_endpoint(db: str) -> bool:
#    """Check the availability of the a SPARQL endpoint.
#
#    :param db: the database to query
#    :returns: True if the endpoint is available, False otherwise.
#    """
#    sparql_query = read_sparql_file(VERSION_QUERY_FILE)
#
#    sparql = SPARQLWrapper(DATABASE_SPARQL_DICT[db])
#
#    sparql.setReturnFormat(JSON)
#    sparql.setQuery(sparql_query)
#
#    try:
#        sparql.queryAndConvert()
#        return True
#    except SPARQLWrapperException:
#        return False


#  def get_version(db: str) -> dict:
#      """Get version of RDF graph.
#
#      :param db: the database to query
#      :returns: a dictionary containing the version information
#      """
#      sparql_query = read_sparql_file(VERSION_QUERY_FILE)
#
#      sparql = SPARQLWrapper(DATABASE_SPARQL_DICT[db])
#      sparql.setReturnFormat(JSON)
#
#      sparql.setQuery(sparql_query)
#      res = sparql.queryAndConvert()
#
#      version = {"source_version": str(res["results"]["bindings"][0]["date"]["value"])}
#
#      return version


def get_aops(
    bridgedb_df: pd.DataFrame, db: str, input_type: str, input_identifier: str
) -> Tuple[pd.DataFrame, dict]:
    """Query for AOPs associated with genes or compounds.

    :param bridgedb_df: BridgeDb output for creating the list of gene/compound ids to query
    :type bridgedb_df: pd.DataFrame
    :param db: the database to query
    :type db: str
    :param input_type: the input type used to query the database (e.g., "gene" or "compound")
    :type input_type: str
    :param input_identifier: the input identifier used to query the database
    :type input_identifier: str
    :returns: a DataFrame containing the AOP Wiki RDF output and a dictionary of the AOP Wiki RDF metadata
    :rtype: Tuple[pd.DataFrame, dict]
    :raises ValueError: If the database (`db`) is not valid or the input type (`input_type`) is not valid
    """
    # Validate inputs
    if db not in DATABASE_SPARQL_DICT:
        raise ValueError(f"{db} is not a valid database.")
    if input_type not in INPUT_OPTIONS:
        raise ValueError(f"{input_type} is not a valid input.")

    # Check if the endpoint is available
    # if not check_endpoint(db=db):
    #    warnings.warn(
    #        f"{db} SPARQL endpoint is not available. Unable to retrieve data.", stacklevel=2
    #    )
    #    return pd.DataFrame(), {}

    # Step 1: Identifier mapping and harmonization
    data_df = get_identifier_of_interest(bridgedb_df, input_identifier)
    # version = get_version(db=db)  # Get the version of the RDF graph

    # Step 2: Prepare target list and batch queries
    target_list = data_df["target"].unique().tolist()
    query_batches = [
        " ".join(f'"{target}"' for target in target_list[i : i + QUERY_LIMIT])
        for i in range(0, len(target_list), QUERY_LIMIT)
    ]

    # Step 3: Run SPARQL queries
    sparql = SPARQLWrapper(DATABASE_SPARQL_DICT[db])
    sparql.setReturnFormat(JSON)

    intermediate_df = pd.DataFrame()
    start_time = datetime.datetime.now()
    col = ""
    for batch in tqdm(query_batches, desc=f"Querying {db} for {input_type}"):
        # Prepare the substitution dictionary
        if input_type == "gene":
            col = AOPWIKIRDF
            substit_dict = {
                "genes": str(
                    [
                        "<https://identifiers.org/ensembl/" + target.replace('"', "") + ">"
                        for target in batch.split(" ")
                    ]
                )
                .replace("[", "")
                .replace("]", "")
                .replace("'", "")
                .replace(",", "")
            }
            query_file = QUERY_GENE
        else:  # input_type == "compound"
            col = AOPWIKI_COMPOUND_COL
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
            query_file = QUERY_COMPOUND

        # Load and substitute the query template
        with open(query_file, "r") as f:
            query = Template(f.read()).substitute(substit_dict)
        # Execute the query and process results
        sparql.setQuery(query)
        res = sparql.queryAndConvert()
        res_df = pd.DataFrame(
            [
                {k: (v["value"] if "value" in v else "") for k, v in item.items()}
                for item in res["results"]["bindings"]
            ]
        )  #
        # Retrieve the expected columns from the SPARQL query results' "vars"
        expected_columns = res["head"]["vars"]

        # Ensure all expected columns are present in intermediate_df
        for col in expected_columns:
            if col not in intermediate_df.columns:
                intermediate_df[col] = None  # Add missing columns with default value None

        # Concatenate the new results into the intermediate DataFrame
        intermediate_df = pd.concat([intermediate_df, res_df], ignore_index=True)

    end_time = datetime.datetime.now()
    # Step 4: Check if the query returned any results
    if intermediate_df.empty:
        warnings.warn(f"There are no results for your input list in {db}.", stacklevel=2)
        return pd.DataFrame(), {}
    # if AOPWIKI_GENE_INPUT_ID or AOPWIKI_COMPOUND_INPUT_ID not in intermediate_df.columns:
    #    warnings.warn(f"There is no annotation for your input list in {db}.", stacklevel=2)
    #    return pd.DataFrame(), {}

    # Step 5: Clean and process the results
    source_namespace = AOPWIKI_GENE_INPUT_ID
    if input_type == "gene":
        input_col = AOPWIKI_GENE_INPUT_ID
        output_dict = AOPWIKI_GENE_OUTPUT_DICT
    else:
        input_col = "pubchem_compound"
        output_dict = AOPWIKI_COMPOUND_OUTPUT_DICT
        source_namespace = "PubChem Compound"
        intermediate_df[input_col] = intermediate_df[input_col].apply(lambda x: x.split("/")[-1])
    for key in output_dict.keys():
        intermediate_df[key] = intermediate_df[key].apply(
            lambda x: x.split("/")[-1] if isinstance(x, str) and "http" in x else x
        )
    intermediate_df.rename(columns={input_col: "target"}, inplace=True)
    col = "target"
    intermediate_df = intermediate_df.drop_duplicates()
    # Strip all text before the last "/" in the 'target' column
    # Step 6: Generate metadata
    metadata_dict = {
        "datasource": db,
        # "metadata": version,
        "query": {
            "size": len(target_list),
            "input_type": input_type,
            "time": str(end_time - start_time),
            "date": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "url": DATABASE_SPARQL_DICT[db],
            "number_of_added_nodes": intermediate_df[col].nunique(),
            "number_of_added_edges": intermediate_df.drop_duplicates(subset=["target", col]).shape[
                0
            ],
        },
    }
    # Step 7: Integrate into main dataframe
    merged_df = collapse_data_sources(
        data_df=data_df,
        source_namespace=source_namespace,
        target_df=intermediate_df,
        common_cols=["target"],
        target_specific_cols=list(output_dict.keys()),
        col_name=AOPWIKIRDF,
    )

    return merged_df, metadata_dict
