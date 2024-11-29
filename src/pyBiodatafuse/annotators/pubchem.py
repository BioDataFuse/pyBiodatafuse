#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Python file for queriying PubChem (https://pubchem.ncbi.nlm.nih.gov/)."""

import datetime
import os
import warnings
from string import Template
from typing import Any, Dict, Tuple

import pandas as pd
from SPARQLWrapper import JSON, SPARQLWrapper
from tqdm import tqdm

from pyBiodatafuse.constants import (
    PUBCHEM,
    PUBCHEM_COMPOUND_ASSAYS_COL,
    PUBCHEM_COMPOUND_INPUT_ID,
    PUBCHEM_COMPOUND_OUTPUT_DICT,
    PUBCHEM_ENDPOINT,
)
from pyBiodatafuse.utils import (
    check_columns_against_constants,
    collapse_data_sources,
    get_identifier_of_interest,
)


def check_endpoint_pubchem() -> bool:
    """Check the availability of the IDSM endpoint.

    :returns: True if the endpoint is available, False otherwise.
    """
    query_string = """SELECT * WHERE {
        <http://rdf.ncbi.nlm.nih.gov/pubchem/taxonomy> ?p ?o
        }
        LIMIT 1
        """
    sparql = SPARQLWrapper(PUBCHEM_ENDPOINT)
    sparql.setOnlyConneg(True)
    sparql.setQuery(query_string)
    try:
        sparql.query()
        return True
    except BaseException:
        return False


# TODO - Add metadata function. Currently, no metadata is returned from IDSM servers


def get_protein_compound_screened(bridgedb_df: pd.DataFrame) -> Tuple[pd.DataFrame, dict]:
    """Query PubChem for molecules screened on proteins as targets.

    :param bridgedb_df: BridgeDb output for creating the list of gene ids to query.
    :returns: a DataFrame containing the PubChem output and dictionary of the PubChem metadata.
    """
    # Check if the IDSM endpoint is available
    api_available = check_endpoint_pubchem()
    if not api_available:
        warnings.warn(
            f"{PUBCHEM} endpoint is not available. Unable to retrieve data.", stacklevel=2
        )
        return pd.DataFrame(), {}

    # Record the start time
    start_time = datetime.datetime.now()

    data_df = get_identifier_of_interest(bridgedb_df, PUBCHEM_COMPOUND_INPUT_ID)
    protein_list_str = data_df["target"].tolist()
    for i in range(len(protein_list_str)):
        protein_list_str[i] = "<http://purl.uniprot.org/uniprot/" + protein_list_str[i] + ">"

    protein_list_str = list(set(protein_list_str))

    query_protein_list = []

    if len(protein_list_str) > 25:
        for i in range(0, len(protein_list_str), 25):
            tmp_list = protein_list_str[i : i + 25]
            query_protein_list.append(" ".join(f"{g}" for g in tmp_list))
    else:
        query_protein_list.append(" ".join(f"{g}" for g in protein_list_str))

    with open(
        os.path.dirname(__file__) + "/queries/pubchem-proteins-screend-molecule.rq", "r"
    ) as fin:
        sparql_query = fin.read()

    sparql = SPARQLWrapper(PUBCHEM_ENDPOINT)
    sparql.setReturnFormat(JSON)
    sparql.setOnlyConneg(True)

    query_count = 0

    intermediate_df = pd.DataFrame()

    for protein_str in tqdm(query_protein_list, desc="Querying PubChem"):
        query_count += 1

        sparql_query_template = Template(sparql_query)
        substit_dict = dict(protein_list=protein_str)
        sparql_query_template_sub = sparql_query_template.substitute(substit_dict)

        sparql.setQuery(sparql_query_template_sub)
        res = sparql.queryAndConvert()

        df = pd.DataFrame(res["results"]["bindings"])
        for col in df:
            df[col] = df[col].map(lambda x: x["value"], na_action="ignore")

        intermediate_df = pd.concat([intermediate_df, df], ignore_index=True)

    # Record the end time
    end_time = datetime.datetime.now()

    """Metdata details"""
    # Get the current date and time
    current_date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    # Calculate the time elapsed
    time_elapsed = str(end_time - start_time)

    # Add the datasource, query, query time, and the date to metadata
    pubchem_metadata: Dict[str, Any] = {
        "datasource": PUBCHEM,
        "query": {
            "size": len(protein_list_str),
            "input_type": PUBCHEM_COMPOUND_INPUT_ID,
            "time": time_elapsed,
            "date": current_date,
            "url": PUBCHEM_ENDPOINT,
        },
    }

    if intermediate_df.empty:
        warnings.warn(
            f"There is no annotation for your input list in {PUBCHEM}.",
            stacklevel=2,
        )
        return pd.DataFrame(), pubchem_metadata

    # Organize the annotation results as an array of dictionaries
    assay_endpoint_types = {
        "http://www.bioassayontology.org/bao#BAO_0000034": "Kd",
        "http://www.bioassayontology.org/bao#BAO_0000186": "AC50",
        "http://www.bioassayontology.org/bao#BAO_0000187": "CC50",
        "http://www.bioassayontology.org/bao#BAO_0000188": "EC50",
        "http://www.bioassayontology.org/bao#BAO_0000190": "IC50",
        "http://www.bioassayontology.org/bao#BAO_0000192": "Ki",
        "http://www.bioassayontology.org/bao#BAO_0002146": "MIC",
    }

    # drop multitarget assays
    intermediate_df.rename(
        columns={
            "upProt": "target",
            "assay": "pubchem_assay_id",
            "SMILES": "smiles",
            "InChI": "inchi",
        },
        inplace=True,
    )

    # Want to keep compounds tested across multiple targets
    intermediate_df["target_count"] = intermediate_df["target_count"].astype(int)
    mask = intermediate_df["target_count"] > 1
    intermediate_df = intermediate_df[mask]

    intermediate_df.drop(columns=["target_count"], inplace=True)

    # Get identifiers from the URL
    intermediate_df["target"] = intermediate_df["target"].map(
        lambda x: x.split("http://purl.uniprot.org/uniprot/")[-1]
    )
    intermediate_df["pubchem_assay_id"] = intermediate_df["pubchem_assay_id"].map(
        lambda x: x.split("http://rdf.ncbi.nlm.nih.gov/pubchem/bioassay/")[-1]
    )

    intermediate_df["outcome"] = intermediate_df["outcome"].map(
        lambda x: x.split("http://rdf.ncbi.nlm.nih.gov/pubchem/vocabulary#")[1]
    )
    intermediate_df["compound_cid"] = intermediate_df["compound_cid"].map(
        lambda x: x.split("http://rdf.ncbi.nlm.nih.gov/pubchem/compound/")[-1]
    )
    intermediate_df["assay_type"] = intermediate_df["assay_type"].map(assay_endpoint_types)

    intermediate_df.drop_duplicates(
        subset=["target", "pubchem_assay_id", "compound_cid"], inplace=True
    )

    # Check if all keys in df match the keys in OUTPUT_DICT
    check_columns_against_constants(
        data_df=intermediate_df,
        output_dict=PUBCHEM_COMPOUND_OUTPUT_DICT,
        check_values_in=["outcome", "inchi"],
    )

    # Merge the two DataFrames on the target column
    merged_df = collapse_data_sources(
        data_df=data_df,
        source_namespace=PUBCHEM_COMPOUND_INPUT_ID,
        target_df=intermediate_df,
        common_cols=["target"],
        target_specific_cols=list(PUBCHEM_COMPOUND_OUTPUT_DICT.keys()),
        col_name=PUBCHEM_COMPOUND_ASSAYS_COL,
    )

    merged_df.reset_index(drop=True, inplace=True)

    """Update metadata"""
    # Calculate the number of new nodes
    num_new_nodes = intermediate_df["compound_cid"].nunique()
    # Calculate the number of new edges
    num_new_edges = intermediate_df.drop_duplicates(subset=["target", "compound_cid"]).shape[0]

    # Check the intermediate_df
    if num_new_edges != len(intermediate_df):
        warnings.warn(
            f"The intermediate_df in {PUBCHEM} annotatur should be checked, please create an issue on https://github.com/BioDataFuse/pyBiodatafuse/issues/.",
            stacklevel=2,
        )

    # Add the number of new nodes and edges to metadata
    pubchem_metadata["query"]["number_of_added_nodes"] = num_new_nodes
    pubchem_metadata["query"]["number_of_added_edges"] = num_new_edges

    return merged_df, pubchem_metadata


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
