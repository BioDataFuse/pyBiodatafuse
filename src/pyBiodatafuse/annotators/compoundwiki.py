# -*- coding: utf-8 -*-

"""Python file for querying the CompoundWiki database (https://compoundcloud.wikibase.cloud/)."""

import datetime
import logging
import os
import warnings
from string import Template
from typing import List, Tuple

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


def query_compoundwiki(compound_ids) -> dict:
    """Query CompoundWiki with a list of compound identifiers.

    :param compound_ids: A list of compound identifiers (e.g., PubChem CIDs, KEGG IDs, etc.)
    :returns: A dictionary mapping each compound ID to its CompoundWiki annotations
    """
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
        print("Running SPARQL query:")
        print(sparql_query)
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
        print("A: ", result_df)

    if result_df.empty or "compound" not in result_df.columns:
        return {}

    result_df["target"] = result_df["compound"].map(compound_to_cid)
    result_df["value"] = result_df.apply(
        lambda row: row["valLabel"] if pd.notna(row.get("valLabel")) else row.get("val"), axis=1
    )

    result_df = result_df.rename(columns={"propLabel": "property"})
    result_df = result_df[result_df["target"].notna()]

    grouped = {}
    for target, df_group in result_df.groupby("target"):
        props = {row["property"]: row["value"] for _, row in df_group.iterrows()}
        label = df_group["compoundLabel"].iloc[0] if "compoundLabel" in df_group.columns else None
        input_id = df_group["target"].iloc[0]
        props.update({
            "compound label": label,
            "input_identifier": input_id,
        })
        grouped[target] = [props]

    return grouped


def inject_compoundwiki_annotations(
    df: pd.DataFrame,
    column_name: str,
    id: str,
    annotation_map: dict
) -> pd.DataFrame:
    """Inject CompoundWiki annotations into nested compound dictionaries in a DataFrame column.

    :param df: DataFrame containing a column with nested lists of compound dictionaries
    :param column_name: Name of the column to inject annotations into
    :param id: Key to extract the compound identifier from each nested dictionary (e.g., 'id', 'id_A')
    :param annotation_map: Dictionary mapping compound IDs to annotation dicts from CompoundWiki
    :returns: Updated DataFrame with CompoundWiki annotations injected
    """
    if column_name not in df.columns:
        return df

    updated_column = []

    for entry in df[column_name]:
        if isinstance(entry, list):
            for compound in entry:
                if isinstance(compound, dict):
                    for key in ["id", "id_A", "id_B"]:
                        cid = compound.get(key)
                        if cid in annotation_map:
                            print(f"Injecting CompoundWiki annotation for {cid}")
                            compound[Cons.COMPOUNDWIKI_COL] = annotation_map[cid]
                            break
        updated_column.append(entry)

    df[column_name] = updated_column
    return df


def get_compound_annotations(
    combined_df: pd.DataFrame,
    kegg_compound_df: pd.DataFrame = None
) -> Tuple[pd.DataFrame, dict]:
    """Annotate compounds in the input DataFrame using CompoundWiki data.

    :param combined_df: Main DataFrame containing compound identifiers and interaction columns
    :param kegg_compound_df: Optional DataFrame for KEGG annotations (if available)
    :returns: Tuple of (annotated DataFrame, metadata dictionary for provenance)
    """
    if not check_endpoint_compoundwiki():
        warnings.warn(f"{Cons.COMPOUNDWIKI} SPARQL endpoint is not available.", stacklevel=2)
        return pd.DataFrame(), {}

    start_time = datetime.datetime.now()
    compoundwiki_version = get_version_compoundwiki()
    annotation_map = {}
    
    # --- Input Identifiers ---
    if Cons.IDENTIFIER_COL in combined_df.columns and Cons.IDENTIFIER_SOURCE_COL in combined_df.columns:
        print("Processing input identifiers for compounds...")

        id_source_unique = combined_df[Cons.IDENTIFIER_SOURCE_COL].dropna().unique()

        if len(id_source_unique) == 1:
            id_source = id_source_unique[0]
            input_id_list = combined_df[Cons.IDENTIFIER_COL].dropna().astype(str).tolist()

            if id_source == "PubChem Compound":
                input_map = query_compoundwiki(input_id_list)

                combined_df[Cons.COMPOUNDWIKI_COL] = combined_df[Cons.IDENTIFIER_COL].map(lambda x: input_map.get(str(x), []))

            else:
                pubchem_ids, id_to_pubchem = query_bridgedb_for_pubchem(input_id_list, id_source)

                if pubchem_ids:
                    input_map = query_compoundwiki(pubchem_ids)
                    combined_df[Cons.COMPOUNDWIKI_COL] = combined_df[Cons.IDENTIFIER_COL].map(
                        lambda x: input_map.get(id_to_pubchem.get(str(x), ""), [])
                    )
                else:
                    combined_df[Cons.COMPOUNDWIKI_COL] = [[] for _ in range(len(combined_df))]

        else:
            combined_df[Cons.COMPOUNDWIKI_COL] = [[] for _ in range(len(combined_df))]

    # --- IntAct block ---
    for intact_col in [Cons.INTACT_INTERACT_COL, Cons.INTACT_COMPOUND_INTERACT_COL]:
        if intact_col in combined_df.columns:
            print(f"Processing IntAct column for compounds: {intact_col}")

            chebi_ids = set()

            # Extract CHEBI: IDs from nested interaction dicts
            for idx, interactions in combined_df[intact_col].dropna().items():
                if isinstance(interactions, list):
                    for interaction in interactions:
                        if isinstance(interaction, dict):
                            for key in ["id_A", "id_B"]:
                                val = interaction.get(key)
                                if isinstance(val, str) and val.startswith("CHEBI:"):
                                    chebi_ids.add(val)

            chebi_ids = list(chebi_ids)
            print(f"Found CHEBI IDs in {intact_col}: {chebi_ids}")

            if chebi_ids:
                # Convert CHEBI IDs to PubChem IDs
                pubchem_ids, id_to_pubchem = query_bridgedb_for_pubchem(chebi_ids, "ChEBI")

                if pubchem_ids:
                    intact_map_cid = query_compoundwiki(pubchem_ids)

                    print(f"Received CompoundWiki annotations for {intact_col}: {intact_map_cid}")

                    intact_map_original = {}
                    for original_id, pubchem_id in id_to_pubchem.items():
                        if pubchem_id in intact_map_cid:
                            intact_map_original[original_id] = intact_map_cid[pubchem_id]

                    annotation_map.update(intact_map_original)

                    combined_df = inject_compoundwiki_annotations(
                        combined_df,
                        intact_col,
                        "id",
                        intact_map_original,
                    )
                    print(f"Injected CompoundWiki annotations into {intact_col}.")
                else:
                    print(f"No PubChem IDs found for {intact_col}.")
            else:
                print(f"No CHEBI IDs found in {intact_col}.")


    # --- KEGG block ---
    if kegg_compound_df is not None and Cons.KEGG_PATHWAY_COL in combined_df.columns:
        print("Processing KEGG column for compounds...")
        kegg_ids = []

        for entry in kegg_compound_df[Cons.KEGG_COMPOUND_COL].dropna():
            if isinstance(entry, list):
                for item in entry:
                    if isinstance(item, dict) and "KEGG_identifier" in item:
                        kegg_ids.append(item["KEGG_identifier"])

        kegg_ids = list(set(kegg_ids))

        if kegg_ids:
            kegg_map = query_compoundwiki(kegg_ids)
            annotation_map.update(kegg_map)

            combined_df = inject_compoundwiki_annotations(
                combined_df,
                Cons.KEGG_PATHWAY_COL,
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


def query_bridgedb_for_pubchem(
    compound_ids: List[str],
    input_datatype: str
) -> Tuple[List[str], dict]:
    """Query BridgeDb to convert compound identifiers to PubChem Compound IDs.

    :param compound_ids: List of original identifiers (e.g., ChEBI, KEGG)
    :param input_datatype: Data source of input IDs (e.g., 'ChEBI', 'KEGG Compound')
    :returns: A tuple of (list of PubChem IDs, mapping from original ID to PubChem ID)
    """

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
    id_to_pubchem = dict(pubchem_df[['identifier', Cons.TARGET_COL]].values)
    print("id_to_pubchem ", id_to_pubchem)

    return pubchem_ids, id_to_pubchem
