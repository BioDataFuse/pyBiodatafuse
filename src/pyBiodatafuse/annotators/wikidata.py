#!/usr/bin/env python
# -*- coding: utf-8 -*-

import datetime
import os
from string import Template

import pandas as pd
from SPARQLWrapper import JSON, SPARQLWrapper

from pyBiodatafuse.utils import collapse_data_sources, get_identifier_of_interest


def get_version_wikidata() -> dict:

    wikidata_version = {"wikidata_version": "0"}

    return wikidata_version


def get_gene_literature(bridgedb_df: pd.DataFrame):
    # Record the start time
    start_time = datetime.datetime.now()

    data_df = get_identifier_of_interest(bridgedb_df, "NCBI Gene")
    gene_list = data_df["target"].tolist()
    gene_list = list(set(gene_list))

    query_gene_lists = []
    if len(gene_list) > 25:
        for i in range(0, len(gene_list), 25):
            tmp_list = gene_list[i : i + 25]
            query_gene_lists.append(" ".join(f'"{g}"' for g in tmp_list))

    else:
        query_gene_lists.append(" ".join(f'"{g}"' for g in gene_list))

    with open(os.path.dirname(__file__) + "/queries/wikidata-genes-literature.rq", "r") as fin:
        sparql_query = fin.read()

    sparql = SPARQLWrapper("https://query.wikidata.org/sparql")
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

        df = pd.DataFrame(res["results"]["bindings"])
        df = df.applymap(lambda x: x["value"])

        results_df_list.append(df)

    # Organize the annotation results as an array of dictionaries
    intermediate_df = pd.concat(results_df_list)
    intermediate_df = intermediate_df.rename(columns={"article": "wikidata_id"})
    intermediate_df = intermediate_df.rename(columns={"geneId": "target"})
    # the next line does some magic
    # before:
    #             target    pubmed       gene wikidata_id
    #         0     1103  10861222  Q14863671   Q22254344
    #         1    85365  11278427  Q18048007   Q24291011
    # after (grouped by geneId):
    #            target  Wikidata_publication
    #         0    1103  [{'pubmed': '10861222', 'wikidata_id': 'Q22254...
    #         4   85365  [{'pubmed': '11278427', 'wikidata_id': 'Q24291...
    intermediate_df = intermediate_df.groupby("target").apply(lambda x: x[["pubmed", "wikidata_id"]].to_dict(orient="records")).reset_index(name="Wikidata_publication")

    # Merge the two DataFrames on the target column
    merged_df = collapse_data_sources(
        data_df=data_df,
        source_namespace="NCBI Gene",
        target_df=intermediate_df,
        common_cols=["target"],
        target_specific_cols=["Wikidata_publication"],
        col_name="Wikidata_publication")

    # Record the end time
    end_time = datetime.datetime.now()

    # Metdata details
    # Get the current date and time
    current_date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    # Calculate the time elapsed
    time_elapsed = str(end_time - start_time)

    # Add version to metadata file
    wikidata_version = get_version_wikidata()

    # Add the datasource, query, query time, and the date to metadata
    wikidata_metadata = {
        "datasource": "Wikidata",
        "metadata": {"source_version": wikidata_version},
        "query": {
            "size": len(gene_list),
            "time": time_elapsed,
            "date": current_date,
            "url": "https://sparql.wikidata.org/sparql",
        },
    }

    return merged_df, wikidata_metadata

def foo():
    intermediate_df.rename(
        columns={"geneId": "target", "geneCount": "pathwayGeneCount"}, inplace=True
    )

    # Merge the two DataFrames on the target column
    merged_df = collapse_data_sources(
        data_df=data_df,
        source_namespace="NCBI Gene",
        target_df=intermediate_df,
        common_cols=["target"],
        target_specific_cols=["pathwayId", "pathwayLabel", "pathwayGeneCount"],
        col_name="WikiPathways",
    )

