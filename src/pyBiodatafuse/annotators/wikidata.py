# -*- coding: utf-8 -*-


"""Python file for querying the Wikidata database (https://www.wikidata.org/)."""

import datetime
import os
import warnings
from string import Template
from typing import Any, Dict

import pandas as pd
from SPARQLWrapper import JSON, SPARQLWrapper
from SPARQLWrapper.SPARQLExceptions import SPARQLWrapperException

import pyBiodatafuse.constants as Cons
from pyBiodatafuse.utils import (
    collapse_data_sources,
    get_identifier_of_interest,
    give_annotator_warning,
)


def check_endpoint_wikidata() -> bool:
    """Check the availability of the Wikidata SPARQL endpoint.

    :returns: True if the endpoint is available, False otherwise.
    """
    with open(os.path.dirname(__file__) + "/queries/wikidata-test.rq", "r") as fin:
        sparql_query = fin.read()

    sparql = SPARQLWrapper(Cons.WIKIDATA_ENDPOINT)
    sparql.setReturnFormat(JSON)

    sparql.setQuery(sparql_query)

    try:
        sparql.queryAndConvert()
        return True
    except SPARQLWrapperException:
        return False


# TODO: Fix this information to be fetched from the server
def get_version_wikidata() -> dict:
    """Get version of Wikidata content.

    :returns: a dictionary containing the (data) version information
    """
    now = str(datetime.datetime.now())

    metadata = {
        "metadata": {
            "data_version": {
                "dataVersion": {
                    "year": now[0:4],
                    "month": now[5:7],
                }
            },
        },
    }

    return metadata


# TODO: Remove this functionabilty
def get_gene_cellular_component(bridgedb_df: pd.DataFrame):
    """Get cellcular component information and Wikidata identifiers for a gene's encoded protein.

    :param bridgedb_df: BridgeDb output for creating the list of gene ids to query
    :returns: a DataFrame containing the Wikidata output and dictionary of the query metadata.
    """
    # Check if the Wikidata API is available
    api_available = check_endpoint_wikidata()

    if not api_available:
        warnings.warn(
            f"{Cons.WIKIDATA} SPARQL endpoint is not available. Unable to retrieve data.",
            stacklevel=2,
        )
        return pd.DataFrame(), {}

    # Record the start time
    start_time = datetime.datetime.now()

    # Add version to metadata file
    wikidata_version = get_version_wikidata()

    data_df = get_identifier_of_interest(bridgedb_df, Cons.NCBI_GENE)
    gene_list = data_df[Cons.TARGET_COL].tolist()
    gene_list = list(set(gene_list))

    query_gene_lists = []
    if len(gene_list) > 25:
        for i in range(0, len(gene_list), 25):
            tmp_list = gene_list[i : i + 25]
            query_gene_lists.append(" ".join(f'"{g}"' for g in tmp_list))

    else:
        query_gene_lists.append(" ".join(f'"{g}"' for g in gene_list))

    with open(
        os.path.dirname(__file__) + "/queries/wikidata-genes-cellularComponent.rq", "r"
    ) as fin:
        sparql_query = fin.read()

    sparql = SPARQLWrapper(Cons.WIKIDATA_ENDPOINT)
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

        df = pd.DataFrame(res["results"]["bindings"])
        for col in df:
            df[col] = df[col].map(lambda x: x["value"], na_action="ignore")

        intermediate_df = pd.concat([intermediate_df, df], ignore_index=True)

    # Record the end time and build metadata
    end_time = datetime.datetime.now()
    current_date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    time_elapsed = str(end_time - start_time)

    # Add the datasource, query, query time, and the date to metadata
    wikidata_metadata: Dict[str, Any] = {
        "datasource": Cons.WIKIDATA,
        "metadata": {"source_version": wikidata_version},
        "query": {
            "size": len(gene_list),
            "time": time_elapsed,
            "date": current_date,
            "url": Cons.WIKIDATA_ENDPOINT,
        },
    }

    if "cellularComp" not in intermediate_df.columns:
        return pd.DataFrame(), wikidata_metadata

    # Organize the annotation results as an array of dictionaries
    intermediate_df = intermediate_df.rename(
        columns={
            "cellularComp": Cons.WIKIDATA_ID_COL,
            "cellularCompLabel": Cons.WIKIDATA_LABEL_COL,
            "go": Cons.GO_ID,
        }
    )

    intermediate_df = intermediate_df.rename(columns={"geneId": Cons.TARGET_COL})

    # Merge the two DataFrames on the target column
    merged_df = collapse_data_sources(
        data_df=data_df,
        source_namespace=Cons.NCBI_GENE,
        target_df=intermediate_df,
        common_cols=[Cons.TARGET_COL],
        target_specific_cols=list(Cons.WIKIDATA_OUTPUT_DICT.keys()),
        col_name=Cons.WIKIDATA_CC_COL,
    )

    # Calculate the number of new nodes
    num_new_nodes = intermediate_df[Cons.WIKIDATA_ID_COL].nunique()
    # Calculate the number of new edges
    num_new_edges = intermediate_df.drop_duplicates(
        subset=[Cons.TARGET_COL, Cons.WIKIDATA_ID_COL]
    ).shape[0]

    # Check the intermediate_df
    if num_new_edges != len(intermediate_df):
        give_annotator_warning(Cons.WIKIDATA)

    # Add the number of new nodes and edges to metadata
    wikidata_metadata[Cons.QUERY][Cons.NUM_NODES] = num_new_nodes
    wikidata_metadata[Cons.QUERY][Cons.NUM_EDGES] = num_new_edges

    return merged_df, wikidata_metadata
