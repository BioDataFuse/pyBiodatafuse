# coding: utf-8

"""Python file for extracting literature data.

The module contains special functions that could be server expensive and can only be performed for smaller datasets.
"""

import os
import warnings
from string import Template
from typing import Dict, Set

import pandas as pd
from SPARQLWrapper import JSON, SPARQLWrapper

from pyBiodatafuse.annotators.wikidata import check_endpoint_wikidata
from pyBiodatafuse.constants import WIKIDATA, WIKIDATA_ENDPOINT, WIKIDATA_GENE_INPUT_ID
from pyBiodatafuse.utils import get_identifier_of_interest


def get_wikidata_gene_literature(bridgedb_df: pd.DataFrame) -> Dict[str, Set[str]]:
    """Get PubMed articles linked to a gene or its encoded protein.

    :param bridgedb_df: BridgeDb output for creating the list of gene ids to query
    :returns: a dictionary with the NCBI gene id as the key and the PMIDs as the value.
    """
    # Check if the Wikidata API is available
    api_available = check_endpoint_wikidata()

    if not api_available:
        warnings.warn(
            f"{WIKIDATA} SPARQL endpoint is not available. Unable to retrieve data.", stacklevel=2
        )
        return {}

    data_df = get_identifier_of_interest(bridgedb_df, WIKIDATA_GENE_INPUT_ID)
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

    sparql = SPARQLWrapper(WIKIDATA_ENDPOINT)
    sparql.setReturnFormat(JSON)

    gene_pmid_dict = {}  # type: Dict[str, Set[str]]

    for gene_list_str in query_gene_lists:
        sparql_query_template = Template(sparql_query)
        substit_dict = dict(gene_list=gene_list_str)
        sparql_query_template_sub = sparql_query_template.substitute(substit_dict)
        sparql.setQuery(sparql_query_template_sub)
        res = sparql.queryAndConvert()

        for gene_result in res["results"]["bindings"]:
            gene_idx = gene_result["geneId"]["value"]

            if gene_idx not in gene_pmid_dict:
                gene_pmid_dict[gene_idx] = set()

            gene_pmid_dict[gene_idx].add("PMID:" + gene_result["pubmed"]["value"])

    return gene_pmid_dict
