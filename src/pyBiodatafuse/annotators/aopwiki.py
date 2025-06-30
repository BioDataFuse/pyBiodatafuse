#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Python file for queriying the AOP Wiki RDF SPARQL endpoint."""

import datetime
import os
import warnings
from string import Template
from typing import Any, Dict, List, Tuple, Union

import pandas as pd
from SPARQLWrapper import JSON, SPARQLWrapper
from tqdm import tqdm

import pyBiodatafuse.constants as Cons
from pyBiodatafuse.utils import (
    collapse_data_sources,
    get_identifier_of_interest,
    give_annotator_warning,
)


def read_sparql_file(file_path: str) -> str:
    """Read a SPARQL query file.

    :param file_path: the path to the SPARQL query file
    :returns: the content of the SPARQL query file
    """
    with open(file_path, "r") as fin:
        sparql_query = fin.read()

    return sparql_query


def check_endpoint() -> bool:
    """Check the availability of the a SPARQL endpoint.

    :returns: True if the endpoint is available, False otherwise.
    """
    #    sparql_query = read_sparql_file(VERSION_QUERY_FILE)
    #    sparql = SPARQLWrapper(DATABASE_SPARQL_DICT[db])
    #    sparql.setReturnFormat(JSON)
    #    sparql.setQuery(sparql_query)
    #    try:
    #        sparql.queryAndConvert()
    #        return True
    #    except SPARQLWrapperException:
    #        return False
    return True


def get_version() -> Dict[str, str]:
    """Get version of RDF graph.

    :returns: a dictionary containing the version information
    """
    # VERSION_QUERY_FILE = os.path.join(os.path.dirname(__file__), "queries", "aopwiki-metadata.rq")
    #  sparql_query = read_sparql_file(VERSION_QUERY_FILE)
    #  sparql = SPARQLWrapper(DATABASE_SPARQL_DICT[db])
    #  sparql.setReturnFormat(JSON)
    #  sparql.setQuery(sparql_query)
    #  res = sparql.queryAndConvert()
    #  version = {"source_version": str(res["results"]["bindings"][0]["date"]["value"])}
    #  return version
    return {"source_version": "Not available"}


def get_aops(
    bridgedb_df: pd.DataFrame, input_type: str, input_identifier: str
) -> Tuple[pd.DataFrame, Dict[str, Any]]:
    """Query for AOPs associated with genes or compounds.

    :param bridgedb_df: BridgeDb output for creating the list of gene/compound ids to query
    :type bridgedb_df: pd.DataFrame
    :param input_type: the input type used to query the database (e.g., "gene" or "compound")
    :type input_type: str
    :param input_identifier: the input identifier used to query the database
    :type input_identifier: str
    :returns: a DataFrame containing the AOP Wiki RDF output and a dictionary of the AOP Wiki RDF metadata
    :rtype: Tuple[pd.DataFrame, dict]
    :raises ValueError: If the input type (`input_type`) is not valid
    """
    # Validate inputs
    if input_type not in Cons.AOPWIKI_INPUT_OPTIONS:
        raise ValueError(f"{input_type} is not a valid input.")

    # Check if the endpoint is available
    if not check_endpoint():
        warnings.warn(
            f"{Cons.AOPWIKIRDF} SPARQL endpoint is not available. Unable to retrieve data.",
            stacklevel=2,
        )
        return pd.DataFrame(), {}

    # Step 1: Identifier mapping and harmonization
    data_df = get_identifier_of_interest(bridgedb_df, input_identifier)
    version = get_version()  # Get the version of the RDF graph

    # Step 2: Prepare target list and batch queries
    target_list = data_df[Cons.TARGET_COL].unique().tolist()
    query_batches = [
        " ".join(f'"{target}"' for target in target_list[i : i + 25])
        for i in range(0, len(target_list), 25)
    ]

    # Step 3: Run SPARQL queries
    sparql = SPARQLWrapper(Cons.AOPWIKI_ENDPOINT)
    sparql.setReturnFormat(JSON)
    sparql.setOnlyConneg(True)

    intermediate_df = pd.DataFrame()
    start_time = datetime.datetime.now()

    # Prepare the substitution dictionary
    if input_type == "gene":
        data_params: Dict[str, Any] = Cons.AOPWIKI_GENE_PARAMS
    else:
        data_params = Cons.AOPWIKI_COMPOUND_PARAMS

    for batch in tqdm(query_batches, desc=f"Querying {Cons.AOPWIKIRDF} for {input_type}"):
        uri = data_params["uri"]
        entity_type = data_params["type_of_input"]

        query_file = f"{os.path.dirname(__file__)}/queries/{data_params['query_file_name']}"
        with open(query_file, "r") as fin:
            sparql_query = fin.read()

        sparql_query_template = Template(sparql_query)
        substit_dict: Dict[str, str] = {
            entity_type: uri + target.replace('"', "") + ">" for target in batch.split(" ")
        }
        substit_dict[entity_type] = (
            str(substit_dict[entity_type])
            .replace("[", "")
            .replace("]", "")
            .replace("'", "")
            .replace(",", "")
        )
        sparql_query_template_sub = sparql_query_template.substitute(substit_dict)

        # Execute the query and process results
        sparql.setQuery(sparql_query_template_sub)
        res = sparql.queryAndConvert()

        # Type check for SPARQL response
        if not isinstance(res, dict) or "results" not in res or "bindings" not in res["results"]:
            continue

        res_df = pd.DataFrame(
            [
                {
                    k: (v["value"] if isinstance(v, dict) and "value" in v else "")
                    for k, v in item.items()
                }
                for item in res["results"]["bindings"]
            ]
        )

        # Retrieve the expected columns from the SPARQL query results' "vars"
        if "head" in res and "vars" in res["head"]:
            expected_columns = res["head"]["vars"]
        else:
            expected_columns = []

        # Ensure all expected columns are present in intermediate_df
        for expected_col in expected_columns:
            if expected_col not in intermediate_df.columns:
                intermediate_df[expected_col] = None  # Add missing columns with default value None

        # Concatenate the new results into the intermediate DataFrame
        intermediate_df = pd.concat([intermediate_df, res_df], ignore_index=True)

    end_time = datetime.datetime.now()

    # Step 4: Check if the query returned any results
    if intermediate_df.empty:
        warnings.warn(
            f"There are no results for your input list in {Cons.AOPWIKIRDF}.", stacklevel=2
        )
        return pd.DataFrame(), {}

    if (
        Cons.AOPWIKI_GENE_INPUT_ID not in intermediate_df.columns
        and Cons.AOPWIKI_COMPOUND_INPUT_ID not in intermediate_df.columns
    ):
        give_annotator_warning(Cons.AOPWIKIRDF)

    # Step 5: Clean and process the results
    source_namespace = data_params["column"]
    input_col = data_params["input_col"]
    output_dict = data_params["output_dict"]

    intermediate_df[input_col] = intermediate_df[input_col].apply(lambda x: x.split("/")[-1])

    for key in output_dict.keys():
        intermediate_df[key] = intermediate_df[key].apply(
            lambda x: x.split("/")[-1] if isinstance(x, str) and "http" in x else x
        )
    intermediate_df.rename(columns={input_col: Cons.TARGET_COL}, inplace=True)
    intermediate_df = intermediate_df.drop_duplicates()

    # Step 6: Generate metadata
    metadata_dict = {
        "datasource": Cons.AOPWIKIRDF,
        "metadata": version,
        "query": {
            "size": len(target_list),
            "input_type": input_type,
            "time": str(end_time - start_time),
            "date": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "url": Cons.AOPWIKI_ENDPOINT,
            "number_of_added_nodes": intermediate_df[Cons.TARGET_COL].nunique(),
            "number_of_added_edges": intermediate_df.drop_duplicates(
                subset=[Cons.TARGET_COL]
            ).shape[0],
        },
    }
    # Step 7: Integrate into main dataframe
    merged_df = collapse_data_sources(
        data_df=data_df,
        source_namespace=source_namespace,
        target_df=intermediate_df,
        common_cols=[Cons.TARGET_COL],
        target_specific_cols=list(output_dict.keys()),
        col_name=Cons.AOPWIKI_COL,
    )

    return merged_df, metadata_dict
