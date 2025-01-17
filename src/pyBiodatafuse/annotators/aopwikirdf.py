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


def get_aops_from_gene(
    bridgedb_df: pd.DataFrame, db: str, input: str, input_identifier: str
) -> Tuple[pd.DataFrame, dict]:
    """Query for AOPs associated with compounds.

    :param bridgedb_df: BridgeDb output for creating the list of gene ids to query
    :param db: the database to query
    :param input: the input type used to query the database
    :param input_identifier: the input identifier used to query the database
    :returns: a DataFrame containing the AOP Wiki RDF output and dictionary of the AOP Wiki RDF metadata.
    """
    assert db in DATABASE_SPARQL_DICT.keys(), f"{db} is not a valid database."

    assert input in INPUT_OPTIONS, f"{input} is not a valid input."

    # Check if the endpoint is available
    api_available = check_endpoint(db=db)

    if not api_available:
        warnings.warn(
            f"{db} SPARQL endpoint is not available. Unable to retrieve data.",
            stacklevel=2,
        )
        return pd.DataFrame(), {}

    # Step 1: Identifier mapping and harmonization
    data_df = get_identifier_of_interest(bridgedb_df, input_identifier)
    version = get_version(db=db)  # Get the version of the RDF graph
    if input == "gene":
        # Step 2: Aggregating query to avoid the query limit
        gene_list = list(data_df["target"].unique())
        query_gene_lists = []
        if len(gene_list) > QUERY_LIMIT:
            for i in range(0, len(gene_list), QUERY_LIMIT):
                tmp_list = gene_list[i : i + QUERY_LIMIT]
                query_gene_lists.append(" ".join(f'"{g}"' for g in tmp_list))
        else:
            query_gene_lists.append(" ".join(f'"{g}"' for g in gene_list))
    
    # Step 3: Running the SPARQL query

    start_time = datetime.datetime.now()  # Record the start time

    sparql = SPARQLWrapper(DATABASE_SPARQL_DICT[db])
    sparql.setReturnFormat(JSON)

    intermediate_df = pd.DataFrame()
    for gene_id in tqdm(query_gene_lists, desc="Querying aopwiki"):
        # Prepare the substitution dictionary
        substit_dict = {'genes': str(['<https://identifiers.org/ensembl/'+i.replace('"','') + '>' for i in gene_id.split(" ")]).replace("[", "").replace("]", "").replace("'", "").replace(",", "")}  #TODO fix
        with open(QUERY_GENE, 'r') as f:
            query = Template(f.read())
        sparql_query = query.substitute(substit_dict)
        
        
        # Set and execute the query
        sparql.setQuery(sparql_query)
        res = sparql.queryAndConvert()
        # Process the results
        
        res = res["results"]["bindings"]
        df = pd.DataFrame(res)
        for col in df:
            df[col] = df[col].map(lambda x: x["value"], na_action="ignore")
        intermediate_df = pd.concat([intermediate_df, df], ignore_index=True)
        
    end_time = datetime.datetime.now()  # Record the end time

    # Step 4: Checking if the query returned any results
    if GENE_INPUT_COL not in intermediate_df.columns:
        warnings.warn(
            f"There is no annotation for your input list in {db}.",
            stacklevel=2,
        )
        return pd.DataFrame(), {}

    # Step 5: Dtype conversion and renaming columns
    intermediate_df.rename(columns={GENE_INPUT_COL: "target"}, inplace=True)
    # Examples:
    # intermediate_df["pathway_gene_count"] = pd.to_numeric(intermediate_df["pathway_gene_count"])
    # intermediate_df["pathway_id"] = intermediate_df["pathway_id"].apply(lambda x: f"WP:{x}")
    intermediate_df = intermediate_df.drop_duplicates()

    # Step 6: Generating a metadata dictionary
    current_date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    time_elapsed = str(end_time - start_time)  # Calculate the time elapsed
    metadata_dict: Dict[str, Any] = {
        "datasource": db,
        "metadata": version,
        "query": {
            "size": len(gene_list),
            "input_type": DATABASE_QUERY_IDENTIFER_GENE,
            "time": time_elapsed,
            "date": current_date,
            "url": DATABASE_SPARQL_DICT[db],
        },
    }

    # Step 7: Cataloging the outputs from the resource
    output_dict = AOPWIKI_GENE_OUTPUT_DICT

    # Step 8: Quality check for dtypes
    #check_columns_against_constants(
    #    data_df=intermediate_df,
    #    output_dict=output_dict,
    #    check_values_in=output_dict.keys())

    # Step 9: Adding node and edge statistics to metadata
    num_new_nodes = intermediate_df[NEW_DATA_COL_GENE].nunique()  # Calculate the number of new nodes
    num_new_edges = intermediate_df.drop_duplicates(subset=["target", NEW_DATA_COL_GENE]).shape[
        0
    ]  # Calculate the number of new edges

    if num_new_edges != len(intermediate_df):
        warnings.warn(
            f"The intermediate_df in {db} annotator should be checked, please create an issue on https://github.com/BioDataFuse/pyBiodatafuse/issues/.",
            stacklevel=2,
        )
    print(num_new_nodes, " nodes")
    metadata_dict["query"]["number_of_added_nodes"] = num_new_nodes
    metadata_dict["query"]["number_of_added_edges"] = num_new_edges
    # Step 10: Integrating into main dataframe
    merged_df = collapse_data_sources(
        data_df=data_df,
        source_namespace=GENE_INPUT_COL,
        target_df=intermediate_df,
        common_cols=["target"],
        target_specific_cols=list(output_dict.keys()),
        col_name=db,
    )

    return merged_df, metadata_dict

