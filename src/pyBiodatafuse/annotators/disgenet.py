# coding: utf-8

"""Python file for queriying DisGeNet database (https://www.disgenet.org/home/)."""

import datetime
import logging
import os
import warnings
from string import Template
from typing import Tuple

import pandas as pd
from SPARQLWrapper import JSON, SPARQLWrapper
from SPARQLWrapper.SPARQLExceptions import SPARQLWrapperException

from pyBiodatafuse.utils import collapse_data_sources, get_identifier_of_interest

logger = logging.getLogger("disgenet")


def test_endpoint_disgenet(endpoint: str) -> bool:
    """Test the availability of the DisGeNET SPARQL endpoint.

    :param endpoint: DisGeNET SAPRQL endpoint ("http://rdf.disgenet.org/sparql/")
    :returns: True if the endpoint is available, False otherwise.
    """
    with open(os.path.dirname(__file__) + "/queries/disgenet-metadata.rq", "r") as fin:
        sparql_query = fin.read()

    sparql = SPARQLWrapper(endpoint)
    sparql.setReturnFormat(JSON)

    sparql.setQuery(sparql_query)

    try:
        sparql.queryAndConvert()
        return True
    except SPARQLWrapperException:
        return False


def get_version_disgenet(endpoint: str) -> dict:
    """Get version of DisGeNET API.

    :param endpoint: DisGeNET SAPRQL endpoint ("http://rdf.disgenet.org/sparql/")
    :returns: a dictionary containing the version information
    """
    with open(os.path.dirname(__file__) + "/queries/disgenet-metadata.rq", "r") as fin:
        sparql_query = fin.read()

    sparql = SPARQLWrapper(endpoint)
    sparql.setReturnFormat(JSON)

    sparql.setQuery(sparql_query)

    version_response = sparql.queryAndConvert()
    bindings = version_response["results"]["bindings"]
    pattern = "RDF Distribution"
    disgenet_version = {"disgenet_version": ""}  # Set default value

    for binding in bindings:
        title = binding["title"]["value"]
        if pattern in title:
            disgenet_version = {"disgenet_version": title}

    return disgenet_version


def get_gene_disease(
    bridgedb_df: pd.DataFrame, endpoint: str = "http://rdf.disgenet.org/sparql/"
) -> Tuple[pd.DataFrame, dict]:
    """Query gene-disease associations from DisGeNET.

    :param bridgedb_df: BridgeDb output for creating the list of gene ids to query.
    :param endpoint: DisGeNET SAPRQL endpoint ("http://rdf.disgenet.org/sparql/").
    :returns: a DataFrame containing the DisGeNET output and dictionary of the DisGeNET metadata.
    """
    # Check if the DisGeNET API is available
    api_available = test_endpoint_disgenet(endpoint=endpoint)
    if not api_available:
        warnings.warn(
            "DisGeNET SPARQL endpoint is not available. Unable to retrieve data.", stacklevel=2
        )
        return pd.DataFrame(), {}

    # Extract the "target" values and join them into a single string separated by commas
    data_df = get_identifier_of_interest(bridgedb_df, "NCBI Gene")
    hgnc_gene_list = data_df["target"].tolist()
    hgnc_gene_list = list(set(hgnc_gene_list))
    query_gene_lists = []

    if len(hgnc_gene_list) > 25:
        for i in range(0, len(hgnc_gene_list), 25):
            tmp_list = hgnc_gene_list[i : i + 25]
            query_gene_lists.append(" ".join(f'"{g}"' for g in tmp_list))

    else:
        query_gene_lists.append(" ".join(f'"{g}"' for g in hgnc_gene_list))

    with open(os.path.dirname(__file__) + "/queries/disgenet-genes-disease.rq", "r") as fin:
        sparql_query = fin.read()

    # Record the start time
    start_time = datetime.datetime.now()

    sparql = SPARQLWrapper(endpoint)
    sparql.setReturnFormat(JSON)

    query_count = 0

    results_df_list = list()

    for gene_list_str in query_gene_lists:
        query_count += 1
        sparql_query_template = Template(sparql_query)
        substit_dict = dict(gene_list=gene_list_str)
        sparql_query_template_sub = sparql_query_template.substitute(substit_dict)
        sparql.setQuery(sparql_query_template_sub)
        res = sparql.queryAndConvert()
        res = res["results"]["bindings"]
        df = pd.DataFrame(res)
        df = df.applymap(lambda x: x["value"])

        results_df_list.append(df)

    # Record the end time
    end_time = datetime.datetime.now()
    # Organize the annotation results as an array of dictionaries
    disgenet_df = pd.concat(results_df_list)
    if "gene_id" not in disgenet_df:
        return pd.DataFrame()
    else:
        disgenet_df.drop_duplicates(inplace=True)
        disgenet_df["target"] = disgenet_df["gene_id"].apply(lambda x: x.split("/")[-1])
        disgenet_df["disease_id"] = disgenet_df["description"].apply(
            lambda x: x.split("[")[1].split("]")[0]
        )
        disgenet_df["disease_label"] = disgenet_df["description"].apply(lambda x: x.split(" [")[0])
        disgenet_df["score"] = disgenet_df["disease_score"].astype(float)
        disgenet_df["source"] = disgenet_df["source"].apply(lambda x: x.split("/")[-1])

        disgenet_df["target"] = disgenet_df["target"].values.astype(str)
        disgenet_df = disgenet_df[["target", "disease_id", "disease_label", "score", "source"]]

        selected_columns = ["disease_id", "disease_label", "score", "source"]

        merged_df = collapse_data_sources(
            data_df=bridgedb_df,
            source_namespace="NCBI Gene",
            target_df=disgenet_df,
            common_cols=["target"],
            target_specific_cols=selected_columns,
            col_name="DisGeNET",
        )

        """Metdata details"""
        # Get the current date and time
        current_date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        # Calculate the time elapsed
        time_elapsed = str(end_time - start_time)
        # Add version to metadata file
        disgenet_version = get_version_disgenet(endpoint=endpoint)
        # Add the datasource, query, query time, and the date to metadata
        disgenet_metadata = {
            "datasource": "DisGeNET",
            "metadata": {"source_version": disgenet_version},
            "query": {
                "size": len(hgnc_gene_list),
                "input_type": "NCBI Gene",
                "time": time_elapsed,
                "date": current_date,
                "url": "http://rdf.disgenet.org/sparql/",
            },
        }

        return merged_df, disgenet_metadata
