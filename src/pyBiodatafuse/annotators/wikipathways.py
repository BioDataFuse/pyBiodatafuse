#!/usr/bin/env python
# -*- coding: utf-8 -*-

import datetime
import os
from string import Template

import pandas as pd
from SPARQLWrapper import JSON, SPARQLWrapper

from pyBiodatafuse.utils import collapse_data_sources, get_identifier_of_interest


def annotateGenesWithWikipathwaysPathways(bridgedb_df: pd.DataFrame):
    """Query WikiPathways for pathways associated with genes.

    @param bridgedb_df: BridgeDb output for creating the list of gene ids to query
    """

    # Record the start time
    start_time = datetime.datetime.now()

    data_df = get_identifier_of_interest(bridgedb_df, "NCBI Gene")
    hgnc_gene_list = data_df["target"].tolist()

    hgnc_gene_list = list(set(hgnc_gene_list))

    query_gene_lists = []

    if len(hgnc_gene_list) > 25:
        print(
            "The number of gene IDs is larger than 25. Therfore, multiple queries will be issued, each with a 25 ID."
        )

        for i in range(0, len(hgnc_gene_list), 25):
            tmp_list = hgnc_gene_list[i : i + 25]
            query_gene_lists.append(" ".join(f'"{g}"' for g in tmp_list))

    else:
        query_gene_lists.append(" ".join(f'"{g}"' for g in hgnc_gene_list))

    with open(os.path.dirname(__file__) + "/queries/wikipathways-genes-pathways.rq", "r") as fin:
        sparql_query = fin.read()

    sparql = SPARQLWrapper("https://sparql.wikipathways.org/sparql")
    sparql.setReturnFormat(JSON)

    query_count = 0

    results_df_list = list()

    for gene_list in query_gene_lists:
        query_count += 1

        sparql_query_template = Template(sparql_query)
        sparql_query_template = sparql_query_template.substitute(gene_list=gene_list)

        sparql.setQuery(sparql_query_template)

        try:
            res = sparql.queryAndConvert()

            df = pd.DataFrame(res["results"]["bindings"])
            df = df.applymap(lambda x: x["value"])

            results_df_list.append(df)

        except Exception as e:
            print(e)

    # Record the end time
    end_time = datetime.datetime.now()

    # Organize the annotation results as an array of dictionaries
    intermediate_df = pd.concat(results_df_list)

    # intermediate_df = intermediate_df.groupby('geneId').apply(lambda x: x.to_dict(orient='r')).rename('WikiPathways')
    # intermediate_df.drop(['pathwayId', 'pathwayLabel', 'geneCount'], axis=1, inplace=True)

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

    # Metdata details
    # Get the current date and time
    current_date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    # Calculate the time elapsed
    time_elapsed = str(end_time - start_time)
    # Add version to metadata file

    with open(os.path.dirname(__file__) + "/queries/wikipathways-metadata.rq", "r") as fin:
        sparql_query = fin.read()

    sparql.setQuery(sparql_query)

    wikipathways_version = ""

    try:
        res = sparql.queryAndConvert()

        wikipathways_version = res["results"]["bindings"][0]["title"]["value"]

    except Exception as e:
        print("SPARQL query exception: " + str(e))

    # Add the datasource, query, query time, and the date to metadata
    wikipathways_metadata = {
        "datasource": "WikiPathways",
        "metadata": {"source_version": str(wikipathways_version)},
        "query": {
            "size": len(hgnc_gene_list),
            "time": time_elapsed,
            "date": current_date,
            "url": "https://sparql.wikipathways.org/sparql",
        },
    }

    return merged_df, wikipathways_metadata
