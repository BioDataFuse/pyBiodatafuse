# coding: utf-8

"""Python file for querying the WikiPathways database (https://www.wikipathways.org/)."""

import datetime
import logging
import os
from string import Template

import pandas as pd
from SPARQLWrapper import JSON, SPARQLWrapper

from pyBiodatafuse.utils import create_or_append_to_metadata, get_identifier_of_interest

logger = logging.getLogger(__name__)


def get_gene_pathway(bridgedb_df: pd.DataFrame) -> pd.DataFrame:
    """Query WikiPathways for pathways associated with genes.

    :param bridgedb_df: BridgeDb output for creating the list of gene ids to query
    :return: DataFrame with WikiPathways pathways associated with genes
    :raises ValueError: if the query fails
    """
    # Record the start time
    start_time = datetime.datetime.now()

    data_df = get_identifier_of_interest(bridgedb_df, "HGNC")
    hgnc_gene_list = data_df["target"].tolist()

    hgnc_gene_list = list(set(hgnc_gene_list))

    query_gene_lists = []

    if len(hgnc_gene_list) > 25:
        logger.info(
            "The number of gene HGNC symbols is larger than 25. \
            Therfore, multiple queries will be issued, each with a 25 symbols."
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

        sparql_query_template = Template(sparql_query).substitute(gene_list=gene_list)
        sparql.setQuery(sparql_query_template)

        try:
            res = sparql.queryAndConvert()

            df = pd.DataFrame(res["results"]["bindings"])
            df = df.applymap(lambda x: x["value"])

            results_df_list.append(df)

        except Exception as e:
            raise ValueError(e)

    # Record the end time
    end_time = datetime.datetime.now()

    intermediate_df = pd.concat(results_df_list)
    intermediate_df = intermediate_df.groupby("geneId").agg(
        WikiPathways=("pathwayId", lambda x: ",".join(x))
    )

    # Merge the two DataFrames on the target column
    # TODO: use collapse_data_sources function for the consistency
    output_df = pd.merge(data_df, intermediate_df, how="left", left_on="target", right_on="geneId")

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
        raise ValueError("SPARQL query exception: " + str(e))

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

    create_or_append_to_metadata(
        wikipathways_metadata
    )  # Call the function from the metadata module

    return output_df
