#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Python file for queriying MolMeDB (https://molmedb.upol.cz/detail/intro)."""

import datetime
import os
import warnings
from math import isnan
from string import Template
from typing import Any, Dict, Tuple

import pandas as pd
from SPARQLWrapper import JSON, SPARQLWrapper

import pyBiodatafuse.constants as Cons
from pyBiodatafuse.utils import (
    check_columns_against_constants,
    collapse_data_sources,
    get_identifier_of_interest,
    give_annotator_warning,
)

pd.set_option("mode.chained_assignment", None)


def check_endpoint_molmedb() -> bool:
    """Check the availability of the MolmeDB endpoint.

    :returns: True if the endpoint is available, False otherwise.
    """
    query_string = """SELECT * WHERE {
        <https://identifiers.org/molmedb/MM00040> <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> ?t
        }
        """
    sparql = SPARQLWrapper(Cons.MOLMEDB_ENDPOINT)
    sparql.setOnlyConneg(True)
    sparql.setQuery(query_string)
    try:
        sparql.query()
        return True
    except BaseException:
        return False


# TODO - Add metadata function. Currently, no metadata is returned from IDSM servers
# def check_version_molmedb():
#     """Check the current version of the MolMeDB database."""
#     pass


def get_gene_compound_inhibitor(bridgedb_df: pd.DataFrame) -> Tuple[pd.DataFrame, dict]:
    """Query MolMeDB for inhibitors of transporters encoded by genes in input.

    :param bridgedb_df: BridgeDb output for creating the list of gene ids to query
    :returns: a DataFrame containing the MolMeDB output and dictionary of the MolMeDB metadata.
    """
    # Check if the MolMeDB endpoint is available
    api_available = check_endpoint_molmedb()
    if not api_available:
        warnings.warn(
            f"{Cons.MOLMEDB} endpoint is not available. Unable to retrieve data.", stacklevel=2
        )
        return pd.DataFrame(), {}

    # Record the start time
    start_time = datetime.datetime.now()

    data_df = get_identifier_of_interest(bridgedb_df, Cons.MOLMEDB_PROTEIN_INPUT_ID)
    molmedb_transporter_list = data_df[Cons.TARGET_COL].tolist()

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

    sparql = SPARQLWrapper(Cons.MOLMEDB_ENDPOINT)
    sparql.setReturnFormat(JSON)
    sparql.setOnlyConneg(True)

    intermediate_df = pd.DataFrame()

    for transporter_list_str in query_transporter_list:
        sparql_query_template = Template(sparql_query)
        substit_dict = dict(transporter_list=transporter_list_str)
        sparql_query_template_sub = sparql_query_template.substitute(substit_dict)

        sparql.setQuery(sparql_query_template_sub)

        res = sparql.queryAndConvert()

        df = pd.DataFrame(res["results"]["bindings"])
        for col in df:
            df[col] = df[col].map(lambda x: x["value"], na_action="ignore")

        if df.empty:
            continue
        # Merging the source_pmid values for each unique compound-gene pair
        cols = [col for col in df.columns.to_list() if col != Cons.SOURCE_PMID]

        df2 = df.groupby(cols)[Cons.SOURCE_PMID].apply(lambda x: ", ".join(x)).reset_index()

        intermediate_df = pd.concat([intermediate_df, df2], ignore_index=True)  # adds to the time

    # Record the end time
    end_time = datetime.datetime.now()

    if intermediate_df.empty:
        warnings.warn(
            f"There is no annotation for your input list in {Cons.MOLMEDB}.",
            stacklevel=2,
        )
        return intermediate_df, {}

    """Metdata details"""
    # Get the current date and time
    current_date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    # Calculate the time elapsed
    time_elapsed = str(end_time - start_time)

    # Add the datasource, query, query time, and the date to metadata
    molmedb_metadata: Dict[str, Any] = {
        "datasource": Cons.MOLMEDB,
        "query": {
            "size": len(molmedb_transporter_list),
            "input_type": Cons.MOLMEDB_PROTEIN_INPUT_ID,
            "time": time_elapsed,
            "date": current_date,
            "url": Cons.MOLMEDB_ENDPOINT,
        },
    }

    # Organize the annotation results as an array of dictionaries
    intermediate_df.rename(
        columns=Cons.MAP_COMPOUND_COL_NAMES,
        inplace=True,
    )

    # Check if all keys in df match the keys in OUTPUT_DICT
    check_columns_against_constants(
        data_df=intermediate_df,
        output_dict=Cons.MOLMEDB_PROTEIN_COMPOUND_OUTPUT_DICT,
        check_values_in=[Cons.MOLMEDB_COMPOUND_DEFAULT_ID],
    )

    # Merge the two DataFrames on the target column
    merged_df = collapse_data_sources(
        data_df=data_df,
        source_namespace=Cons.MOLMEDB_PROTEIN_INPUT_ID,
        target_df=intermediate_df,
        common_cols=[Cons.TARGET_COL],
        target_specific_cols=list(Cons.MOLMEDB_PROTEIN_COMPOUND_OUTPUT_DICT.keys()),
        col_name=Cons.MOLMEDB_PROTEIN_COMPOUND_COL,
    )

    # Ensuring all the dictionaries in the MolMeDB_transporter_inhibitor column are same for multiple gene isoforms
    main_df = []

    for source in merged_df[Cons.IDENTIFIER_COL].unique():
        mm = merged_df[merged_df[Cons.IDENTIFIER_COL] == source]
        if len(mm) < 2:
            main_df.append(mm)
            continue

        molmedb_output = list(mm[Cons.MOLMEDB_PROTEIN_COMPOUND_COL].values)
        unique_output = get_unique_dicts(molmedb_output)
        mm[Cons.MOLMEDB_PROTEIN_COMPOUND_COL] = ([unique_output]) * len(mm)
        main_df.append(mm)

    main_df = pd.concat(main_df)

    """Update metadata"""
    # Calculate the number of new nodes
    num_new_nodes = intermediate_df[Cons.MOLMEDB_ID].nunique()
    # Calculate the number of new edges
    num_new_edges = intermediate_df.drop_duplicates(
        subset=[Cons.TARGET_COL, Cons.MOLMEDB_ID]
    ).shape[0]

    # Check the intermediate_df
    if num_new_edges != len(intermediate_df):
        warnings.warn(
            f"The intermediate_df in {Cons.MOLMEDB} annotator should be checked, please create an issue on https://github.com/BioDataFuse/pyBiodatafuse/issues/.",
            stacklevel=2,
        )

    # Add the number of new nodes and edges to metadata
    molmedb_metadata[Cons.QUERY][Cons.NUM_NODES] = num_new_nodes
    molmedb_metadata[Cons.QUERY][Cons.NUM_EDGES] = num_new_edges

    return main_df, molmedb_metadata


def get_compound_gene_inhibitor(bridgedb_df: pd.DataFrame) -> Tuple[pd.DataFrame, dict]:
    """Query MolMeDB for transporters inhibited by molecule.

    :param bridgedb_df: BridgeDb output for creating the list of gene ids to query.
    :returns: a DataFrame containing the MolMeDB output and dictionary of the MolMeDB metadata.
    """
    # Check if the MolMeDB endpoint is available
    api_available = check_endpoint_molmedb()
    if not api_available:
        warnings.warn(
            f"{Cons.MOLMEDB} endpoint is not available. Unable to retrieve data.", stacklevel=2
        )
        return pd.DataFrame(), {}

    # Record the start time
    start_time = datetime.datetime.now()

    data_df = get_identifier_of_interest(bridgedb_df, Cons.MOLMEDB_COMPOUND_INPUT_ID)
    inhibitor_list_str = data_df[Cons.TARGET_COL].tolist()

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

    sparql = SPARQLWrapper(Cons.MOLMEDB_ENDPOINT)
    sparql.setReturnFormat(JSON)
    sparql.setOnlyConneg(True)

    intermediate_df = pd.DataFrame()

    for inhibitor_list_str in query_inhibitor_list:
        sparql_query_template = Template(sparql_query)
        substit_dict = dict(inhibitor_list=inhibitor_list_str)
        sparql_query_template_sub = sparql_query_template.substitute(substit_dict)

        sparql.setQuery(sparql_query_template_sub)

        res = sparql.queryAndConvert()

        df = pd.DataFrame(res["results"]["bindings"])

        for col in df:
            df[col] = df[col].map(lambda x: x["value"], na_action="ignore")

        # Merging the source_pmid values for each unique compound-gene pair
        df2 = (
            df.groupby(
                [
                    Cons.MOLMEDB_INHIBITOR_INCHIKEY,
                    Cons.MOLMEDB_UNIPROT_TREMBL_ID,
                    Cons.MOLMEDB_HGNC_ID,
                ]
            )[Cons.SOURCE_PMID]
            .apply(lambda x: ", ".join(x))
            .reset_index()
        )
        df2[Cons.MOLMEDB_UNIPROT_TREMBL_ID] = df2[Cons.MOLMEDB_UNIPROT_TREMBL_ID].apply(
            lambda x: f"{Cons.UNIPROT_TREMBL}:{x}"
        )

        intermediate_df = pd.concat([intermediate_df, df2], ignore_index=True)  # adds to the time

    # Record the end time
    end_time = datetime.datetime.now()

    if intermediate_df.empty:
        warnings.warn(
            f"There is no annotation for your input list in {Cons.MOLMEDB}.",
            stacklevel=2,
        )
        return pd.DataFrame(), {}

    # Organize the annotation results as an array of dictionaries
    intermediate_df.rename(columns=Cons.MAP_GENE_COL_NAMES, inplace=True)

    """Metdata details"""
    # Get the current date and time
    current_date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    # Calculate the time elapsed
    time_elapsed = str(end_time - start_time)

    # Add the datasource, query, query time, and the date to metadata
    molmedb_metadata: Dict[str, Any] = {
        "datasource": Cons.MOLMEDB,
        "query": {
            "size": len(inhibitor_list_str),
            "input_type": Cons.MOLMEDB_COMPOUND_INPUT_ID,
            "time": time_elapsed,
            "date": current_date,
            "url": Cons.MOLMEDB_ENDPOINT,
        },
    }

    # Check if all keys in df match the keys in OUTPUT_DICT
    check_columns_against_constants(
        data_df=intermediate_df,
        output_dict=Cons.MOLMEDB_COMPOUND_PROTEIN_OUTPUT_DICT,
        check_values_in=[Cons.UNIPROT_TREMBL],
    )

    # Merge the two DataFrames on the target column
    merged_df = collapse_data_sources(
        data_df=data_df,
        source_namespace=Cons.MOLMEDB_COMPOUND_INPUT_ID,
        target_df=intermediate_df,
        common_cols=[Cons.TARGET_COL],
        target_specific_cols=list(Cons.MOLMEDB_COMPOUND_PROTEIN_OUTPUT_DICT.keys()),
        col_name=Cons.MOLMEDB_COMPOUND_PROTEIN_COL,
    )

    """Update metadata"""
    # Calculate the number of new nodes
    num_new_nodes = intermediate_df[Cons.MOLMEDB_UNIPROT_ID].nunique()
    # Calculate the number of new edges
    num_new_edges = intermediate_df.drop_duplicates(
        subset=[Cons.TARGET_COL, Cons.MOLMEDB_UNIPROT_ID]
    ).shape[0]

    # Check the intermediate_df
    if num_new_edges != len(intermediate_df):
        give_annotator_warning(Cons.MOLMEDB)

    # Add the number of new nodes and edges to metadata
    molmedb_metadata[Cons.QUERY][Cons.NUM_NODES] = num_new_nodes  #
    molmedb_metadata[Cons.QUERY][Cons.NUM_EDGES] = num_new_edges  # noqa

    return merged_df, molmedb_metadata


def get_unique_dicts(list_of_list_of_dicts: list) -> list:
    """Return list of unique dictionaries.

    :param list_of_list_of_dicts: list of list of dictionaries
    :returns: list of unique dictionaries
    """
    seen = set()
    unique_dicts = []

    for list_of_dicts in list_of_list_of_dicts:
        for d in list_of_dicts:
            # Convert dictionary to frozenset of its items (which is hashable)
            dict_items = frozenset(d.items())
            if not isinstance(d[Cons.MOLMEDB_INCHIKEY], str) and isnan(d[Cons.MOLMEDB_INCHIKEY]):
                empty = d
                continue
            if dict_items not in seen:
                seen.add(dict_items)
                unique_dicts.append(d)

    if len(unique_dicts) == 0:
        return [empty]
    return unique_dicts
