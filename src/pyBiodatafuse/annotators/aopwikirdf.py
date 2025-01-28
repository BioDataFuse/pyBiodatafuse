#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Python file for queriying the AOP Wiki RDF SPARQL endpoint ()."""

import datetime
import os
import warnings
from string import Template
from typing import Any, Dict, Tuple

import pandas as pd
from SPARQLWrapper import JSON, SPARQLWrapper
from SPARQLWrapper.SPARQLExceptions import SPARQLWrapperException
from tqdm import tqdm

from pyBiodatafuse.utils import (
    check_columns_against_constants,
    collapse_data_sources,
    get_identifier_of_interest,
)
from pyBiodatafuse.constants import (
    OPENTARGETS_GO_COL,
    AOPWIKI_ENDPOINT,
    AOPWIKI_GENE_INPUT_ID,
    AOPWIKI_GENE_OUTPUT_DICT,
    AOPWIKI_COMPOUND_OUTPUT_DICT
)
# Pre-requisite:
VERSION_QUERY_FILE = os.path.join(os.path.dirname(__file__), "queries", "aopwiki-metadata.rq")
DATABASE_SPARQL_DICT = {"aopwiki": AOPWIKI_ENDPOINT}
DATABASE_QUERY_IDENTIFER_GENE = AOPWIKI_GENE_INPUT_ID
QUERY_LIMIT = 25
QUERY_COMPOUND = os.path.join(os.path.dirname(__file__), "queries", "aopwiki-get-by-compound.rq.rq")
QUERY_GENE = os.path.join(os.path.dirname(__file__), "queries", "aopwiki-get-by-gene.rq")
#QUERY_PROCESS = os.path.join(os.path.dirname(__file__), "queries", "aopwiki-get-by-biological-process.rq.rq")
GENE_INPUT_COL = "Ensembl"
COMPOUND_INPUT = list()
#PROCESS_INPUT_COL = OPENTARGETS_GO_COL
NEW_DATA_COL_GENE = "aop"

# Unique inputs and outputs:
INPUT_OPTIONS = ["gene", "compound",]  #"biological_process"]


def read_sparql_file(file_path: str) -> str:
    """Read a SPARQL query file.

    :param file_path: the path to the SPARQL query file
    :returns: the content of the SPARQL query file
    """
    with open(file_path, "r") as fin:
        sparql_query = fin.read()

    return sparql_query


def check_endpoint(db: str) -> bool:
    """Check the availability of the a SPARQL endpoint.

    :param db: the database to query
    :returns: True if the endpoint is available, False otherwise.
    """
    sparql_query = read_sparql_file(VERSION_QUERY_FILE)

    sparql = SPARQLWrapper(DATABASE_SPARQL_DICT[db])

    sparql.setReturnFormat(JSON)
    sparql.setQuery(sparql_query)

    try:
        sparql.queryAndConvert()
        return True
    except SPARQLWrapperException:
        return False


def get_version(db: str) -> dict:
    """Get version of RDF graph.

    :param db: the database to query
    :returns: a dictionary containing the version information
    """
    sparql_query = read_sparql_file(VERSION_QUERY_FILE)

    sparql = SPARQLWrapper(DATABASE_SPARQL_DICT[db])
    sparql.setReturnFormat(JSON)

    sparql.setQuery(sparql_query)
    res = sparql.queryAndConvert()

    version = {"source_version": str(res["results"]["bindings"][0]["date"]["value"])}

    return version


def get_aops(
    bridgedb_df: pd.DataFrame, db: str, input_type: str, input_identifier: str
) -> Tuple[pd.DataFrame, dict]:
    """Query for AOPs associated with genes or compounds.

    :param bridgedb_df: BridgeDb output for creating the list of gene/compound ids to query
    :param db: the database to query
    :param input_type: the input type used to query the database (e.g., "gene" or "compound")
    :param input_identifier: the input identifier used to query the database
    :returns: a DataFrame containing the AOP Wiki RDF output and dictionary of the AOP Wiki RDF metadata.
    """
    # Validate inputs
    if db not in DATABASE_SPARQL_DICT:
        raise ValueError(f"{db} is not a valid database.")
    if input_type not in INPUT_OPTIONS:
        raise ValueError(f"{input_type} is not a valid input.")

    # Check if the endpoint is available
    if not check_endpoint(db=db):
        warnings.warn(f"{db} SPARQL endpoint is not available. Unable to retrieve data.", stacklevel=2)
        return pd.DataFrame(), {}

    # Step 1: Identifier mapping and harmonization
    data_df = get_identifier_of_interest(bridgedb_df, input_identifier)
    version = get_version(db=db)  # Get the version of the RDF graph

    # Step 2: Prepare target list and batch queries
    target_list = data_df["target"].unique().tolist()
    query_batches = [
        " ".join(f'"{target}"' for target in target_list[i:i + QUERY_LIMIT])
        for i in range(0, len(target_list), QUERY_LIMIT)
    ]

    # Step 3: Run SPARQL queries
    sparql = SPARQLWrapper(DATABASE_SPARQL_DICT[db])
    sparql.setReturnFormat(JSON)

    intermediate_df = pd.DataFrame()
    start_time = datetime.datetime.now()

    for batch in tqdm(query_batches, desc=f"Querying {db} for {input_type}"):
        # Prepare the substitution dictionary
        if input_type == "gene":
            substit_dict = {
                'genes': ' '.join(f'<https://identifiers.org/ensembl/{target.replace('"', '')}>' for target in batch.split())
            }
            query_file = QUERY_GENE
        else:  # input_type == "compound"
            substit_dict = {
                'compounds': ' '.join(f'<https://identifiers.org/pubchem.compound/{target.replace('"', '')}>' for target in batch.split())
            }
            query_file = QUERY_COMPOUND

        # Load and substitute the query template
        with open(query_file, 'r') as f:
            query = Template(f.read()).substitute(substit_dict)

        # Execute the query and process results
        sparql.setQuery(query)
        res = sparql.queryAndConvert()
        res_df = pd.DataFrame([{k: v["value"] for k, v in item.items()} for item in res["results"]["bindings"]])
        intermediate_df = pd.concat([intermediate_df, res_df], ignore_index=True)

    end_time = datetime.datetime.now()

    # Step 4: Check if the query returned any results
    if GENE_INPUT_COL not in intermediate_df.columns:
        warnings.warn(f"There is no annotation for your input list in {db}.", stacklevel=2)
        return pd.DataFrame(), {}

    # Step 5: Clean and process the results
    intermediate_df.rename(columns={GENE_INPUT_COL: "target"}, inplace=True)
    intermediate_df = intermediate_df.drop_duplicates()

    # Step 6: Generate metadata
    metadata_dict = {
        "datasource": db,
        "metadata": version,
        "query": {
            "size": len(target_list),
            "input_type": input_type,
            "time": str(end_time - start_time),
            "date": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "url": DATABASE_SPARQL_DICT[db],
            "number_of_added_nodes": intermediate_df[NEW_DATA_COL_GENE].nunique(),
            "number_of_added_edges": intermediate_df.drop_duplicates(subset=["target", NEW_DATA_COL_GENE]).shape[0],
        },
    }

    # Step 7: Integrate into main dataframe
    merged_df = collapse_data_sources(
        data_df=data_df,
        source_namespace=GENE_INPUT_COL,
        target_df=intermediate_df,
        common_cols=["target"],
        target_specific_cols=list(AOPWIKI_GENE_OUTPUT_DICT.keys() if input_type == "gene" else AOPWIKI_COMPOUND_OUTPUT_DICT.keys()),
        col_name=db,
    )

    return merged_df, metadata_dict