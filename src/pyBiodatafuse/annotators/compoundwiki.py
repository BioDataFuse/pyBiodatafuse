# -*- coding: utf-8 -*-

"""Python file for querying the CompoundWiki database (https://compoundcloud.wikibase.cloud/)."""

import datetime
import logging
import os
import warnings
from string import Template

import pandas as pd
from SPARQLWrapper import JSON, SPARQLWrapper
from SPARQLWrapper.SPARQLExceptions import SPARQLWrapperException

from pyBiodatafuse.constants import COMPOUNDWIKI, COMPOUNDWIKI_COL, COMPOUNDWIKI_ENDPOINT
from pyBiodatafuse.utils import collapse_data_sources, get_identifier_of_interest


def check_endpoint_compoundwiki() -> bool:
    """Check the availability of the CompoundWiki SPARQL endpoint."""
    with open(os.path.dirname(__file__) + "/queries/compoundwiki-test.rq", "r") as fin:
        sparql_query = fin.read()

    sparql = SPARQLWrapper(COMPOUNDWIKI_ENDPOINT)
    sparql.setReturnFormat(JSON)
    sparql.setQuery(sparql_query)

    try:
        sparql.queryAndConvert()
        return True
    except SPARQLWrapperException:
        return False


def get_version_compoundwiki() -> dict:
    """Return metadata with a timestamp for versioning."""
    now = str(datetime.datetime.now())
    return {
        "metadata": {
            "data_version": {
                "dataVersion": {
                    "year": now[0:4],
                    "month": now[5:7],
                }
            },
        },
    }


def get_compound_annotations(bridgedb_df: pd.DataFrame):
    """Query CompoundWiki for annotations using PubChem CID.

    :param bridgedb_df: DataFrame with a 'PubChem CID' column
    :return: Annotated DataFrame and query metadata

    """
    # Check if the CompoundWiki API is available
    api_available = check_endpoint_compoundwiki()

    if not api_available:
        warnings.warn(
            f"{COMPOUNDWIKI} SPARQL endpoint is not available. Unable to retrieve data.",
            stacklevel=2,
        )
        return pd.DataFrame(), {}

    # Record the start time
    start_time = datetime.datetime.now()

    # Add version to metadata file
    compoundwiki_version = get_version_compoundwiki()

    data_df = get_identifier_of_interest(bridgedb_df, "PubChem Compound")

    cid_list = list(set(data_df["target"].tolist()))

    query_batches = []
    if len(cid_list) > 25:
        for i in range(0, len(cid_list), 25):
            tmp = cid_list[i : i + 25]
            query_batches.append(" ".join(f'"{cid}"' for cid in tmp))
    else:
        query_batches.append(" ".join(f'"{cid}"' for cid in cid_list))

    with open(os.path.dirname(__file__) + "/queries/compoundwiki-compounds.rq", "r") as fin:
        sparql_query_template = Template(fin.read())

    sparql = SPARQLWrapper(COMPOUNDWIKI_ENDPOINT)
    sparql.setReturnFormat(JSON)
    sparql.setMethod("POST")
    sparql.addCustomHttpHeader("User-Agent", "pyBiodatafuse/1.0")

    intermediate_df = pd.DataFrame()
    compound_to_cid = {}

    for _idx, batch in enumerate(query_batches):
        sparql_query = sparql_query_template.substitute(compound_list=batch)
        sparql.setQuery(sparql_query)

        res = sparql.queryAndConvert()
        bindings = res["results"]["bindings"]

        df = pd.DataFrame(bindings)

        for col in df:
            df[col] = df[col].map(lambda x: x["value"], na_action="ignore")

        for _, row in df.iterrows():
            if row.get("propLabel") == "PubChem CID":
                compound_to_cid[row["compound"]] = row["val"]

        intermediate_df = pd.concat([intermediate_df, df], ignore_index=True)

    if "compound" not in intermediate_df.columns:
        return pd.DataFrame(), {"datasource": COMPOUNDWIKI, "metadata": compoundwiki_version}

    intermediate_df["target"] = intermediate_df["compound"].map(compound_to_cid)
    intermediate_df["value"] = intermediate_df.apply(
        lambda row: row["valLabel"] if pd.notna(row.get("valLabel")) else row.get("val"), axis=1
    )

    intermediate_df = intermediate_df.rename(columns={"propLabel": "property"})
    intermediate_df = intermediate_df[intermediate_df["target"].notna()]

    annotations_df = (
        intermediate_df.groupby("target")
        .apply(lambda df: [{row["property"]: row["value"] for _, row in df.iterrows()}])
        .reset_index(name=COMPOUNDWIKI_COL)
    )

    merged_df = pd.merge(data_df, annotations_df, on="target", how="left")

    end_time = datetime.datetime.now()

    compoundwiki_metadata = {
        "datasource": COMPOUNDWIKI,
        "metadata": {"source_version": compoundwiki_version},
        "query": {
            "size": len(cid_list),
            "time": str(end_time - start_time),
            "date": end_time.strftime("%Y-%m-%d %H:%M:%S"),
            "url": COMPOUNDWIKI_ENDPOINT,
        },
    }

    return merged_df, compoundwiki_metadata
