#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Python file for queriying Bgee database (https://bgee.org) via its SPARQL endpoint (https://www.bgee.org/sparql/)."""

import datetime
import os
from string import Template

import pandas as pd
from SPARQLWrapper import JSON, SPARQLWrapper

from pyBiodatafuse.utils import collapse_data_sources, get_identifier_of_interest


def get_version_bgee() -> dict:
    """Get version of Bgee RDF data from its SPARQL endpoint.

    # not sure if a version per-se can be retrieved, but the endpoint supports
    # http://purl.org/dc/terms/modified
    :returns: a dictionary containing the last modified date information
    """
    with open(os.path.dirname(__file__) + "/queries/bgee-get-last-modified.rq", "r") as fin:
        sparql_query = fin.read()

    sparql = SPARQLWrapper("https://www.bgee.org/sparql/")
    sparql.setReturnFormat(JSON)

    sparql.setQuery(sparql_query)
    res = sparql.queryAndConvert()

    bgee_version = {"bgee_version": res["results"]["bindings"][0]["date_modified"]["value"]}

    return bgee_version


def get_gene_expression(bridgedb_df: pd.DataFrame, anatomical_entities: pd.DataFrame):
    """Query gene-tissue expression information from Bgee.

    :param bridgedb_df: BridgeDb output for creating the list of gene ids to query
    :param anatomical_entities: a dataframe containing the names of Anatomical entities of interest
    :returns: a DataFrame containing the Bgee output and dictionary of the Bgee metadata.
    """
    # Record the start time
    start_time = datetime.datetime.now()

    data_df = get_identifier_of_interest(bridgedb_df, "Ensembl")
    gene_list = data_df["target"].tolist()
    gene_list = list(set(gene_list))

    query_gene_lists = []
    if len(gene_list) > 25:
        for i in range(0, len(gene_list), 25):
            tmp_list = gene_list[i : i + 25]
            query_gene_lists.append(" ".join(f'"{g}"' for g in tmp_list))

    else:
        query_gene_lists.append(" ".join(f'"{g}"' for g in gene_list))

    anat_entities_list = anatomical_entities["AnatomicalEntityNames"].tolist()
    anat_entities_list = list(set(anat_entities_list))

    query_anat_entities_lists = []
    if len(anat_entities_list) > 25:
        for i in range(0, len(anat_entities_list), 25):
            tmp_list = anat_entities_list[i : i + 25]
            query_anat_entities_lists.append(" ".join(f'"{g}"' for g in tmp_list))

    else:
        query_anat_entities_lists.append(" ".join(f'"{g}"' for g in anat_entities_list))

    with open(os.path.dirname(__file__) + "/queries/bgee-genes-tissues-expression.rq", "r") as fin:
        sparql_query = fin.read()

    sparql = SPARQLWrapper("https://www.bgee.org/sparql/")
    sparql.setReturnFormat(JSON)

    query_count = 0

    results_df_list = list()

    for gene_list_str in query_gene_lists:
        for query_anat_entities_str in query_anat_entities_lists:
            query_count += 1

            sparql_query_template = Template(sparql_query)
            substit_dict = dict(gene_list=gene_list_str, anat_entities_list=query_anat_entities_str)
            sparql_query_template_sub = sparql_query_template.substitute(substit_dict)

            sparql.setQuery(sparql_query_template_sub)
            res = sparql.queryAndConvert()

            df = pd.DataFrame(res["results"]["bindings"])
            df = df.applymap(lambda x: x["value"])

            results_df_list.append(df)

    # Organize the annotation results as an array of dictionaries
    intermediate_df = pd.concat(results_df_list)

    intermediate_df.rename(columns={"ensembl_id": "target"}, inplace=True)

    # Record the end time
    end_time = datetime.datetime.now()

    # Metadata details
    # Get the current date and time
    current_date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    # Calculate the time elapsed
    time_elapsed = str(end_time - start_time)

    # Add version to metadata file
    bgee_version = get_version_bgee()

    # Add the datasource, query, query time, and the date to metadata
    bgee_metadata = {
        "datasource": "Bgee",
        "metadata": {"source_version": bgee_version},
        "query": {
            "size": len(gene_list),
            "time": time_elapsed,
            "date": current_date,
            "url": "https://www.bgee.org/sparql/",
        },
    }

    # Merge the two DataFrames on the target column
    merged_df = collapse_data_sources(
        data_df=data_df,
        source_namespace="Ensembl",
        target_df=intermediate_df,
        common_cols=["target"],
        target_specific_cols=[
            "anatomical_entity_id",
            "anatomical_entity_name",
            "expression_level",
            "confidence_level",
        ],
        col_name="Bgee",
    )

    return merged_df, bgee_metadata
