#!/usr/bin/env python
# -*- coding: utf-8 -*-

import datetime
import os
from string import Template

import pandas as pd
from SPARQLWrapper import JSON, SPARQLWrapper

from pyBiodatafuse.utils import collapse_data_sources, get_identifier_of_interest


def get_version_bgee() -> dict:
    # not sure if a version per-se can be retrieved, but the endpoint supports
    # http://purl.org/dc/terms/modified

    with open(os.path.dirname(__file__) + "/queries/bgee-get-last-modified.rq", "r") as fin:
        sparql_query = fin.read()

    sparql = SPARQLWrapper("https://www.bgee.org/sparql/")
    sparql.setReturnFormat(JSON)

    sparql.setQuery(sparql_query)
    res = sparql.queryAndConvert()

    bgee_version = {"bgee_version": res["results"]["bindings"][0]["dateModified"]["value"]}

    return bgee_version


def get_gene_literature(bridgedb_df: pd.DataFrame):
    # Record the start time
    start_time = datetime.datetime.now()
    print(start_time)

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

    with open(os.path.dirname(__file__) + "/queries/bgee-genes-tissues-expression.rq", "r") as fin:
        sparql_query = fin.read()

    sparql = SPARQLWrapper("https://www.bgee.org/sparql/")
    sparql.setReturnFormat(JSON)

    query_count = 0

    results_df_list = list()

    for gene_list_str in query_gene_lists:
        query_count += 1

        sparql_query_template = Template(sparql_query)
        substit_dict = dict(gene_list=gene_list_str)
        sparql_query_template_sub = sparql_query_template.substitute(substit_dict)
        print(sparql_query_template_sub)

        sparql.setQuery(sparql_query_template_sub)
        res = sparql.queryAndConvert()

        df = pd.DataFrame(res["results"]["bindings"])
        df = df.applymap(lambda x: x["value"])

        results_df_list.append(df)

    # Organize the annotation results as an array of dictionaries
    intermediate_df = pd.concat(results_df_list)
    print(intermediate_df)

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

    return data_df, bgee_metadata