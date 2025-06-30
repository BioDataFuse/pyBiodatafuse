# -*- coding: utf-8 -*-

"""Python file for querying the CompoundWiki database (https://compoundcloud.wikibase.cloud/)."""

import datetime
import logging
import os
import warnings
from string import Template
from typing import List

import pandas as pd
from SPARQLWrapper import JSON, SPARQLWrapper
from SPARQLWrapper.SPARQLExceptions import SPARQLWrapperException

import pyBiodatafuse.constants as Cons
from pyBiodatafuse import id_mapper
from pyBiodatafuse.utils import collapse_data_sources, get_identifier_of_interest


def check_endpoint_compoundwiki() -> bool:
    """Check the availability of the CompoundWiki SPARQL endpoint."""
    with open(os.path.dirname(__file__) + "/queries/compoundwiki-test.rq", "r") as fin:
        sparql_query = fin.read()

    sparql = SPARQLWrapper(Cons.COMPOUNDWIKI_ENDPOINT)
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


def query_compoundwiki(compound_ids: list[str]) -> dict:
    """Query CompoundWiki with a list of compound identifiers."""
    from string import Template

    query_batches = []
    if len(compound_ids) > 25:
        for i in range(0, len(compound_ids), 25):
            tmp = compound_ids[i : i + 25]
            query_batches.append(" ".join(f'"{cid}"' for cid in tmp))
    else:
        query_batches.append(" ".join(f'"{cid}"' for cid in compound_ids))

    with open(os.path.dirname(__file__) + "/queries/compoundwiki-compounds.rq", "r") as fin:
        sparql_query_template = Template(fin.read())

    sparql = SPARQLWrapper(Cons.COMPOUNDWIKI_ENDPOINT)
    sparql.setReturnFormat(JSON)
    sparql.setMethod("POST")
    sparql.addCustomHttpHeader("User-Agent", "pyBiodatafuse/1.0")

    result_df = pd.DataFrame()
    compound_to_cid = {}

    for batch in query_batches:
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

        result_df = pd.concat([result_df, df], ignore_index=True)

    if result_df.empty or "compound" not in result_df.columns:
        return {}

    result_df["target"] = result_df["compound"].map(compound_to_cid)
    result_df["value"] = result_df.apply(
        lambda row: row["valLabel"] if pd.notna(row.get("valLabel")) else row.get("val"), axis=1
    )

    result_df = result_df.rename(columns={"propLabel": "property"})
    result_df = result_df[result_df["target"].notna()]

    grouped = (
        result_df.groupby("target")
        .apply(
            lambda df: [
                {
                    **{row["property"]: row["value"] for _, row in df.iterrows()},
                    "compound label": (
                        df["compoundLabel"].iloc[0] if "compoundLabel" in df.columns else None
                    ),
                }
            ]
        )
        .to_dict()
    )

    return grouped



def inject_compoundwiki_annotations(df: pd.DataFrame, column_name: str, id_key: str, annotation_map: dict) -> pd.DataFrame:
    """Inject CompoundWiki annotations into a dataframe column with compound lists of dicts."""
    if column_name not in df.columns:
        return df

    def annotate_entry(entry):
        if isinstance(entry, list):
            for compound in entry:
                if isinstance(compound, dict) and id_key in compound:
                    cid = compound[id_key]
                    if cid in annotation_map:
                        compound[Cons.COMPOUNDWIKI_COL] = annotation_map[cid]
        return entry

    df[column_name] = df[column_name].apply(annotate_entry)
    return df


def get_compound_annotations(combined_df: pd.DataFrame, kegg_compound_df: pd.DataFrame = None):
    if not check_endpoint_compoundwiki():
        warnings.warn(f"{Cons.COMPOUNDWIKI} SPARQL endpoint is not available.", stacklevel=2)
        return pd.DataFrame(), {}

    start_time = datetime.datetime.now()
    compoundwiki_version = get_version_compoundwiki()
    annotation_map = {}

    # --- IntAct block ---
    if Cons.INTACT_INTERACT_COL in combined_df.columns:
        print("Processing IntAct column for CHEBI compounds...")

        chebi_ids = set()

        for idx, interactions in combined_df[Cons.INTACT_INTERACT_COL].dropna().items():
            if isinstance(interactions, list):
                for interaction in interactions:
                    if isinstance(interaction, dict):
                        for key in ["id_A", "id_B"]:
                            val = interaction.get(key)
                            if isinstance(val, str) and val.startswith("CHEBI:"):
                                chebi_ids.add(val)

        chebi_ids = list(chebi_ids)
        print(f"Found CHEBI IDs in IntAct: {chebi_ids}")

        if chebi_ids:
            pubchem_ids = query_bridgedb_for_pubchem(chebi_ids, "ChEBI")

            if pubchem_ids:
                intact_map = query_compoundwiki(pubchem_ids)
                print(f"Received CompoundWiki annotations for IntAct: {intact_map}")

                for cid in chebi_ids:
                    if cid not in intact_map:
                        intact_map[cid] = None

                annotation_map.update(intact_map)

                combined_df = inject_compoundwiki_annotations(
                    combined_df,
                    Cons.INTACT_INTERACT_COL,
                    "id",
                    intact_map,
                )
                print("Injected CompoundWiki annotations into combined_df IntAct column.")
            else:
                print("No PubChem IDs found to query CompoundWiki.")
        else:
            print("No CHEBI IDs found in IntAct column.")



    # --- KEGG block ---
    if kegg_compound_df is not None and Cons.KEGG_COMPOUND_OUTPUT_DICT_KEY in kegg_compound_df.columns:
        kegg_ids = [
            entry["KEGG_identifier"]
            for entry in kegg_compound_df[Cons.KEGG_COMPOUND_OUTPUT_DICT_KEY].dropna()
            if isinstance(entry, dict) and "KEGG_identifier" in entry
        ]
        kegg_ids = list(set(kegg_ids))
        if kegg_ids:
            kegg_map = query_compoundwiki(kegg_ids)
            annotation_map.update(kegg_map)
            combined_df = inject_compoundwiki_annotations(
                combined_df,
                Cons.KEGG_COMPOUND_OUTPUT_DICT_KEY,
                "KEGG_identifier",
                kegg_map,
            )

    end_time = datetime.datetime.now()

    compoundwiki_metadata = {
        "datasource": Cons.COMPOUNDWIKI,
        "metadata": {"source_version": compoundwiki_version},
        "query": {
            "size": sum(len(v) if v is not None else 0 for v in annotation_map.values()),
            "time": str(end_time - start_time),
            "date": end_time.strftime("%Y-%m-%d %H:%M:%S"),
            "url": Cons.COMPOUNDWIKI_ENDPOINT,
        },
    }

    return combined_df, compoundwiki_metadata


def query_bridgedb_for_pubchem(compound_ids: List[str], input_datatype: str) -> List[str]:

    print(f"Querying BridgeDb for PubChem IDs from {input_datatype}...")

    data_input = pd.DataFrame(compound_ids, columns=["identifier"])

    bridgedb_df, bridgedb_metadata = id_mapper.bridgedb_xref(
        identifiers=data_input,
        input_species="Human",
        input_datasource=input_datatype,
        output_datasource="All"
    )

    if bridgedb_df.empty:
        print("BridgeDb returned no results.")
        return []

    # Filter for rows with PubChem Compound in target.source
    pubchem_df = bridgedb_df[bridgedb_df["target.source"] == Cons.COMPOUNDWIKI_COMPOUND_INPUT_ID]

    if pubchem_df.empty:
        print("No PubChem Compound IDs found in BridgeDb result.")
        return []

    pubchem_ids = list(set(pubchem_df[Cons.TARGET_COL].tolist()))
    print(f"Found {len(pubchem_ids)} unique PubChem Compound IDs.")

    return pubchem_ids
