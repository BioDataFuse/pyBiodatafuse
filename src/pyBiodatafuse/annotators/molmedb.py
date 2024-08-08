#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Python file for queriying MolMeDB (https://molmedb.upol.cz/detail/intro)."""

import datetime
import os
import warnings
from string import Template
from typing import Tuple

import numpy as np
import pandas as pd
from SPARQLWrapper import JSON, SPARQLWrapper

from pyBiodatafuse.constants import (
    MOLMEDB,
    MOLMEDB_COMPOUND_INPUT_ID,
    MOLMEDB_COMPOUND_OUTPUT_DICT,
    MOLMEDB_ENDPOINT,
    MOLMEDB_GENE_INPUT_ID,
    MOLMEDB_GENE_OUTPUT_DICT,
    MOLMEDB_INHIBITED_COL,
    MOLMEDB_INHIBITOR_COL,
)
from pyBiodatafuse.utils import (
    check_columns_against_constants,
    collapse_data_sources,
    get_identifier_of_interest,
)


def check_endpoint_molmedb() -> bool:
    """Check the availability of the MolmeDB endpoint.

    :returns: True if the endpoint is available, False otherwise.
    """
    query_string = """SELECT * WHERE {
        <https://identifiers.org/molmedb/MM00040> <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> ?t
        }
        """
    sparql = SPARQLWrapper(MOLMEDB_ENDPOINT)
    sparql.setOnlyConneg(True)
    sparql.setQuery(query_string)
    try:
        sparql.query()
        return True
    except BaseException:
        return False


# TODO - Add metadata function. Currently, no metadata is returned from IDSM servers


def get_gene_compound_inhibitor(bridgedb_df: pd.DataFrame) -> Tuple[pd.DataFrame, dict]:
    """Query MolMeDB for inhibitors of transporters encoded by genes.

    :param bridgedb_df: BridgeDb output for creating the list of gene ids to query
    :returns: a DataFrame containing the MolMeDB output and dictionary of the MolMeDB metadata.
    """
    # Check if the MolMeDB endpoint is available
    api_available = check_endpoint_molmedb()
    if not api_available:
        warnings.warn("MolMeDB endpoint is not available. Unable to retrieve data.", stacklevel=2)
        return pd.DataFrame(), {}

    # Record the start time
    start_time = datetime.datetime.now()

    data_df = get_identifier_of_interest(bridgedb_df, MOLMEDB_GENE_INPUT_ID)
    molmedb_transporter_list = data_df["target"].tolist()

    molmedb_transporter_list = list(set(molmedb_transporter_list))

    query_transporter_list = []

    if len(molmedb_transporter_list) > 25:
        for i in range(0, len(molmedb_transporter_list), 25):
            tmp_list = molmedb_transporter_list[i : i + 25]
            query_transporter_list.append(" ".join(f'"{g}"' for g in tmp_list))

    else:
        query_transporter_list.append(" ".join(f'"{g}"' for g in molmedb_transporter_list))

    with open(
        os.path.dirname(__file__) + "/queries/molmedb-transporters-inhibitors.rq", "r"
    ) as fin:
        sparql_query = fin.read()

    sparql = SPARQLWrapper(MOLMEDB_ENDPOINT)
    sparql.setReturnFormat(JSON)
    sparql.setOnlyConneg(True)

    query_count = 0

    intermediate_df = pd.DataFrame()

    for transporter_list_str in query_transporter_list:
        query_count += 1

        sparql_query_template = Template(sparql_query)
        substit_dict = dict(transporter_list=transporter_list_str)
        sparql_query_template_sub = sparql_query_template.substitute(substit_dict)

        sparql.setQuery(sparql_query_template_sub)

        res = sparql.queryAndConvert()

        df = pd.DataFrame(res["results"]["bindings"])
        df = df.applymap(lambda x: x["value"], na_action="ignore")

        intermediate_df = pd.concat([intermediate_df, df], ignore_index=True)  # adds to the time

    # Record the end time
    end_time = datetime.datetime.now()

    if "transporterID" not in intermediate_df.columns:
        return pd.DataFrame(), {"datasource": MOLMEDB, "metadata": ""}

    # Organize the annotation results as an array of dictionaries
    intermediate_df.rename(
        columns={
            "transporterID": "target",
            "pubchem_compound_id": "compound_cid",
            "label": "compound_name",
        },
        inplace=True,
    )
    intermediate_df["source_doi"] = intermediate_df["source_doi"].map(
        lambda x: "doi:" + x, na_action="ignore"
    )

    # Check if all keys in df match the keys in OUTPUT_DICT
    check_columns_against_constants(
        data_df=intermediate_df,
        output_dict=MOLMEDB_GENE_OUTPUT_DICT,
        check_values_in=["molmedb_id", "source_doi", "drugbank_id"],
    )

    # Merge the two DataFrames on the target column
    merged_df = collapse_data_sources(
        data_df=data_df,
        source_namespace=MOLMEDB_GENE_INPUT_ID,
        target_df=intermediate_df,
        common_cols=["target"],
        target_specific_cols=list(MOLMEDB_GENE_OUTPUT_DICT.keys()),
        col_name=MOLMEDB_INHIBITOR_COL,
    )

    # if mappings exist but SPARQL returns empty response
    if (not merged_df.empty) and merged_df[MOLMEDB_INHIBITOR_COL][0] is None:
        merged_df.drop_duplicates(subset=["identifier", MOLMEDB_INHIBITOR_COL], inplace=True)

    elif not merged_df.empty:
        res_keys = merged_df[MOLMEDB_INHIBITOR_COL][0][0].keys()
        # remove duplicate identifier and response row
        merged_df[MOLMEDB_INHIBITOR_COL] = merged_df[MOLMEDB_INHIBITOR_COL].map(
            lambda x: tuple(frozenset(d.items()) for d in x), na_action="ignore"
        )
        merged_df.drop_duplicates(subset=["identifier", MOLMEDB_INHIBITOR_COL], inplace=True)
        merged_df[MOLMEDB_INHIBITOR_COL] = merged_df[MOLMEDB_INHIBITOR_COL].map(
            lambda res_tup: list(dict((x, y) for x, y in res) for res in res_tup),
            na_action="ignore",
        )

        # drop rows with duplicate identifiers with empty response
        identifiers = merged_df["identifier"].unique()
        for identifier in identifiers:
            if merged_df.loc[merged_df["identifier"] == identifier].shape[0] > 1:
                mask = merged_df[MOLMEDB_INHIBITOR_COL].apply(
                    lambda lst: all(
                        [
                            all([isinstance(val, float) and np.isnan(val) for val in dct.values()])
                            for dct in lst
                        ]
                    )
                )
                mask2 = merged_df["identifier"].apply(lambda x, id=identifier: x == id)
                merged_df.drop(merged_df[mask & mask2].index, inplace=True)

        # set default order to response dictionaries to keep output consistency
        merged_df[MOLMEDB_INHIBITOR_COL] = merged_df[MOLMEDB_INHIBITOR_COL].apply(
            lambda res: list(dict((k, r[k]) for k in res_keys) for r in res)
        )
        # set numerical identifiers to int to kepp output consistency
        merged_df[MOLMEDB_INHIBITOR_COL] = merged_df[MOLMEDB_INHIBITOR_COL].apply(
            lambda res: int_response_value_types(res, ["compound_cid", "source_pmid"])
        )
    merged_df.reset_index(drop=True, inplace=True)

    # Metadata details
    # Get the current date and time
    current_date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    # Calculate the time elapsed
    time_elapsed = str(end_time - start_time)
    # Calculate the number of new nodes
    num_new_nodes = intermediate_df["molmedb_id"].nunique()
    # Calculate the number of new edges
    num_edges = len(intermediate_df[["target", "molmedb_id"]].drop_duplicates())

    # Add the datasource, query, query time, and the date to metadata
    molmedb_metadata = {
        "datasource": MOLMEDB,
        "query": {
            "size": len(molmedb_transporter_list),
            "input_type": MOLMEDB_GENE_INPUT_ID,
            "number_of_added_nodes": num_new_nodes,
            "number_of_added_edges": num_edges,
            "time": time_elapsed,
            "date": current_date,
            "url": MOLMEDB_ENDPOINT,
        },
    }

    return merged_df, molmedb_metadata


def get_compound_gene_inhibitor(bridgedb_df: pd.DataFrame) -> Tuple[pd.DataFrame, dict]:
    """Query MolMeDB for transporters inhibited by molecule.

    :param bridgedb_df: BridgeDb output for creating the list of gene ids to query.
    :returns: a DataFrame containing the MolMeDB output and dictionary of the MolMeDB metadata.
    """
    # Check if the MolMeDB endpoint is available
    api_available = check_endpoint_molmedb()
    if not api_available:
        warnings.warn("MolMeDB endpoint is not available. Unable to retrieve data.", stacklevel=2)
        return pd.DataFrame(), {}

    # Record the start time
    start_time = datetime.datetime.now()

    data_df = get_identifier_of_interest(bridgedb_df, MOLMEDB_COMPOUND_INPUT_ID)
    inhibitor_list_str = data_df["target"].tolist()

    inhibitor_list_str = list(set(inhibitor_list_str))

    query_inhibitor_list = []

    if len(inhibitor_list_str) > 25:
        for i in range(0, len(inhibitor_list_str), 25):
            tmp_list = inhibitor_list_str[i : i + 25]
            query_inhibitor_list.append(" ".join(f'"{g}"' for g in tmp_list))

    else:
        query_inhibitor_list.append(" ".join(f'"{g}"' for g in inhibitor_list_str))

    with open(
        os.path.dirname(__file__) + "/queries/molmedb-transporters-inhibited-by-molecule.rq", "r"
    ) as fin:
        sparql_query = fin.read()

    sparql = SPARQLWrapper(MOLMEDB_ENDPOINT)
    sparql.setReturnFormat(JSON)
    sparql.setOnlyConneg(True)

    query_count = 0

    intermediate_df = pd.DataFrame()

    for inhibitor_list_str in query_inhibitor_list:
        query_count += 1

        sparql_query_template = Template(sparql_query)
        substit_dict = dict(inhibitor_list=inhibitor_list_str)
        sparql_query_template_sub = sparql_query_template.substitute(substit_dict)

        sparql.setQuery(sparql_query_template_sub)

        res = sparql.queryAndConvert()

        df = pd.DataFrame(res["results"]["bindings"])
        df = df.applymap(lambda x: x["value"], na_action="ignore")

        intermediate_df = pd.concat([intermediate_df, df], ignore_index=True)  # adds to the time

    # Record the end time
    end_time = datetime.datetime.now()

    if "inhibitorInChIKey" not in intermediate_df.columns:
        return pd.DataFrame(), {"datasource": MOLMEDB, "metadata": ""}

    # Organize the annotation results as an array of dictionaries
    intermediate_df.rename(
        columns={"inhibitorInChIKey": "target", "hgcn_id": "hgnc_symbol"}, inplace=True
    )
    intermediate_df["source_doi"] = intermediate_df["source_doi"].map(
        lambda x: "doi:" + x, na_action="ignore"
    )

    # Check if all keys in df match the keys in OUTPUT_DICT
    check_columns_against_constants(
        data_df=intermediate_df,
        output_dict=MOLMEDB_COMPOUND_OUTPUT_DICT,
        check_values_in=["UNIPROT_TREMBL_ID"],
    )

    # Merge the two DataFrames on the target column
    merged_df = collapse_data_sources(
        data_df=data_df,
        source_namespace=MOLMEDB_COMPOUND_INPUT_ID,
        target_df=intermediate_df,
        common_cols=["target"],
        target_specific_cols=list(MOLMEDB_COMPOUND_OUTPUT_DICT.keys()),
        col_name=MOLMEDB_INHIBITED_COL,
    )

    # Metdata details
    # Get the current date and time
    current_date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    # Calculate the time elapsed
    time_elapsed = str(end_time - start_time)
    # Calculate the number of new nodes
    num_new_nodes = intermediate_df["uniprot_trembl_id"].nunique()
    # Calculate the number of new edges
    num_edges = len(intermediate_df)

    # Add the datasource, query, query time, and the date to metadata
    molmedb_metadata = {
        "datasource": MOLMEDB,
        "query": {
            "size": len(inhibitor_list_str),
            "input_type": MOLMEDB_COMPOUND_INPUT_ID,
            "number_of_added_nodes": num_new_nodes,
            "number_of_added_edges": num_edges,
            "time": time_elapsed,
            "date": current_date,
            "url": MOLMEDB_ENDPOINT,
        },
    }

    return merged_df, molmedb_metadata


def int_response_value_types(resp_list: list, key_list: list):
    """Change values in response dictionaries to int to stay consistent with other Annotators.

    :param: resp_list: list of response dictionaries.
    :param: key_list: list of keys to change to int.
    :returns: resp_list with int values in response dictionaries on keys in key_list.
    """
    for r in resp_list:
        for k in key_list:
            try:
                r[k] = int(r[k])
            except ValueError:
                continue
    return resp_list
