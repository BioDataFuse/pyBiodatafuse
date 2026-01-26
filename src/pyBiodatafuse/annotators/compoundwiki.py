# -*- coding: utf-8 -*-

"""Python file for querying the CompoundWiki database (https://compoundcloud.wikibase.cloud/)."""

import datetime
import os
import warnings
from string import Template
from typing import List, Tuple

import numpy as np
import pandas as pd
import requests
from SPARQLWrapper import JSON, SPARQLWrapper
from SPARQLWrapper.SPARQLExceptions import SPARQLWrapperException

import pyBiodatafuse.constants as Cons
from pyBiodatafuse import id_mapper


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
    """Get version of CompoundWiki from the MediaWiki API.

    Fetches the most recent change from CompoundWiki to get version information
    including revision ID (revid) and timestamp.

    :returns: a dictionary containing the version information with revid as source_version
    """
    # CompoundWiki MediaWiki API endpoint for recent changes
    api_url = "https://compoundcloud.wikibase.cloud/w/api.php"
    params = {
        "action": "query",
        "list": "recentchanges",
        "rcprop": "title|ids|sizes|flags|user|timestamp",
        "rclimit": "1",
        "format": "json",
    }

    try:
        response = requests.get(api_url, params=params, timeout=10)
        response.raise_for_status()
        data = response.json()

        # Extract version information from the API response
        if "query" in data and "recentchanges" in data["query"]:
            recent_changes = data["query"]["recentchanges"]
            if recent_changes:
                recent_change = recent_changes[0]
                revid = recent_change.get("revid", "")
                timestamp = recent_change.get("timestamp", "")

                # Use revid as the source_version, with timestamp as additional info
                # Format: "revid_30029 (2026-01-25)"
                if timestamp:
                    dt = datetime.datetime.fromisoformat(timestamp.replace("Z", "+00:00"))
                    date_str = dt.strftime("%Y-%m-%d")
                    source_version = f"revid_{revid} ({date_str})"
                else:
                    source_version = f"revid_{revid}"

                return {"source_version": source_version}

    except (requests.RequestException, KeyError, ValueError) as e:
        warnings.warn(
            f"Could not fetch CompoundWiki version from API: {e}. Using current timestamp.",
            stacklevel=2,
        )

    # Fallback to current timestamp if API call fails
    now = datetime.datetime.now()
    return {"source_version": f"unknown ({now.strftime('%Y-%m-%d')})"}


def query_compoundwiki(compound_ids) -> dict:
    """Query CompoundWiki with a list of compound identifiers.

    :param compound_ids: A list of compound identifiers (e.g., PubChem CIDs, KEGG IDs, etc.)
    :returns: A dictionary mapping each compound ID to its CompoundWiki annotations
    """
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

    grouped = {}
    for target, df_group in result_df.groupby("target"):
        props = {row["property"]: row["value"] for _, row in df_group.iterrows()}
        label = df_group["compoundLabel"].iloc[0] if "compoundLabel" in df_group.columns else None
        input_id = df_group["target"].iloc[0]
        props.update(
            {
                "compound label": label,
                "input_identifier": input_id,
            }
        )
        grouped[target] = [props]

    return grouped


def inject_compoundwiki_annotations(
    df: pd.DataFrame, column_name: str, id_key: str, annotation_map: dict
) -> pd.DataFrame:
    """Inject CompoundWiki annotations into nested compound dictionaries in a DataFrame column.

    :param df: DataFrame containing a column with nested lists of compound dictionaries
    :param column_name: Name of the column to inject annotations into
    :param id_key: Key to extract the compound identifier from each nested dictionary (e.g., 'id', 'id_A')
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
                    compound_id = compound.get(id_key)
                    if compound_id in annotation_map:
                        compound[Cons.COMPOUNDWIKI_COL] = annotation_map[compound_id]
        updated_column.append(entry)

    df[column_name] = updated_column
    return df


def get_compound_annotations(
    combined_df: pd.DataFrame, kegg_compound_df: pd.DataFrame = None
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

    empty_annotation = [{key: np.nan for key in Cons.COMPOUNDWIKI_OUTPUT_DICT.keys()}]

    # --- Input Identifiers ---
    if (
        Cons.IDENTIFIER_COL in combined_df.columns
        and Cons.IDENTIFIER_SOURCE_COL in combined_df.columns
    ):

        id_source_unique = combined_df[Cons.IDENTIFIER_SOURCE_COL].dropna().unique()

        if len(id_source_unique) == 1:
            id_source = id_source_unique[0]
            input_id_list = combined_df[Cons.IDENTIFIER_COL].dropna().astype(str).tolist()

            annotations = []

            if id_source == "PubChem-compound":
                input_map = query_compoundwiki(input_id_list)
                annotations = [
                    input_map.get(str(x), empty_annotation)
                    for x in combined_df[Cons.IDENTIFIER_COL]
                ]
            else:
                pubchem_ids, id_to_pubchem = query_bridgedb_for_pubchem(
                    list(input_id_list), id_source
                )
                if pubchem_ids:
                    input_map = query_compoundwiki(pubchem_ids)
                    annotations = [
                        input_map.get(id_to_pubchem.get(str(x), ""), empty_annotation)
                        for x in combined_df[Cons.IDENTIFIER_COL]
                    ]
                else:
                    annotations = [empty_annotation for _ in range(len(combined_df))]

            if any(annotation != empty_annotation for annotation in annotations):
                combined_df[Cons.COMPOUNDWIKI_COL] = annotations

    # --- OpenTargets ---
    if Cons.OPENTARGETS_GENE_COMPOUND_COL in combined_df.columns:
        print(
            f"Processing {Cons.OPENTARGETS} column for compounds: {Cons.OPENTARGETS_GENE_COMPOUND_COL}"
        )

        chembl_ids = []

        # Collect identifiers
        for entry in combined_df[Cons.OPENTARGETS_GENE_COMPOUND_COL].dropna():
            if isinstance(entry, list):
                for compound in entry:
                    if isinstance(compound, dict) and Cons.CHEMBL_ID in compound:
                        chembl_id = compound[Cons.CHEMBL_ID]
                        if isinstance(chembl_id, str) and chembl_id.startswith("CHEMBL:"):
                            chembl_ids.append(chembl_id.replace("CHEMBL:", ""))

        # Deduplicate
        chembl_ids = list(set(chembl_ids))

        if chembl_ids:
            # Convert identifiers to PubChem identifiers
            pubchem_ids, chembl_to_pubchem = query_bridgedb_for_pubchem(
                list(chembl_ids), "ChEMBL compound"
            )

            if pubchem_ids:
                # Query CompoundWiki with the PubChem identifiers
                opentargets_map = query_compoundwiki(pubchem_ids)
                annotation_map.update(opentargets_map)

                # Prepare mapping for injection using the original identifiers
                injection_map = {
                    f"CHEMBL:{chembl_id}": opentargets_map[pcid]
                    for chembl_id, pcid in chembl_to_pubchem.items()
                    if pcid in opentargets_map
                }

                # Inject annotations into the nested dictionaries
                combined_df = inject_compoundwiki_annotations(
                    combined_df,
                    Cons.OPENTARGETS_GENE_COMPOUND_COL,
                    Cons.CHEMBL_ID,
                    injection_map,
                )

    # --- IntAct ---
    for intact_col in [Cons.INTACT_INTERACT_COL, Cons.INTACT_COMPOUND_INTERACT_COL]:
        if intact_col in combined_df.columns:
            print(f"Processing {Cons.INTACT} column for compounds: {intact_col}")
            chebi_ids = set()

            # Collect all identifiers
            for _idx, interactions in combined_df[intact_col].dropna().items():
                if isinstance(interactions, list):
                    for interaction in interactions:
                        if isinstance(interaction, dict):
                            for key in [Cons.INTACT_ID_A, Cons.INTACT_ID_B]:
                                val = interaction.get(key)
                                if isinstance(val, str) and val.startswith("CHEBI:"):
                                    chebi_ids.add(val)

            chebi_ids_list = list(chebi_ids)

            if chebi_ids_list:
                pubchem_ids, id_to_pubchem = query_bridgedb_for_pubchem(chebi_ids_list, "ChEBI")
                if pubchem_ids:
                    intact_map_cid = query_compoundwiki(pubchem_ids)

                    intact_map_original = {}
                    for chebi_id, pubchem_id in id_to_pubchem.items():
                        if pubchem_id in intact_map_cid:
                            annotated = []
                            for ann in intact_map_cid[pubchem_id]:
                                ann_copy = ann.copy()
                                ann_copy["input_identifier"] = chebi_id
                                annotated.append(ann_copy)

                            intact_map_original[chebi_id] = annotated

                    annotation_map.update(intact_map_original)

                    for key in [Cons.INTACT_ID_A, Cons.INTACT_ID_B]:
                        combined_df = inject_compoundwiki_annotations(
                            combined_df,
                            intact_col,
                            key,
                            intact_map_original,
                        )

    # --- KEGG ---
    if kegg_compound_df is not None and Cons.KEGG_PATHWAY_COL in combined_df.columns:
        print(f"Processing {Cons.KEGG} column for compounds: {Cons.KEGG_PATHWAY_COL}")
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

    # --- PubChem ---
    if Cons.PUBCHEM_COMPOUND_ASSAYS_COL in combined_df.columns:
        print(f"Processing {Cons.PUBCHEM} column for compounds: {Cons.PUBCHEM_COMPOUND_ASSAYS_COL}")

        pubchem_ids = []

        for entry in combined_df[Cons.PUBCHEM_COMPOUND_ASSAYS_COL].dropna():
            if isinstance(entry, list):
                for compound in entry:
                    if (
                        isinstance(compound, dict)
                        and "compound_cid" in compound
                        and isinstance(compound["compound_cid"], str)
                    ):
                        cid = compound["compound_cid"]
                        if cid.startswith("CID:"):
                            cid = cid.replace("CID:", "")
                        pubchem_ids.append(cid)

        pubchem_ids = list(set(pubchem_ids))

        if pubchem_ids:
            pubchem_map = query_compoundwiki(pubchem_ids)
            annotation_map.update(pubchem_map)

            combined_df = inject_compoundwiki_annotations(
                combined_df,
                Cons.PUBCHEM_COMPOUND_ASSAYS_COL,
                "compound_cid",
                {
                    compound["compound_cid"]: pubchem_map[cid]
                    for entry in combined_df[Cons.PUBCHEM_COMPOUND_ASSAYS_COL].dropna()
                    for compound in (entry if isinstance(entry, list) else [])
                    if (
                        isinstance(compound, dict)
                        and "compound_cid" in compound
                        and isinstance(compound["compound_cid"], str)
                        and (cid := compound["compound_cid"].replace("CID:", "")) in pubchem_map
                    )
                },
            )

    # --- MolMeDB ---
    for molmedb_col in [Cons.MOLMEDB_PROTEIN_COMPOUND_COL, Cons.MOLMEDB_COMPOUND_PROTEIN_COL]:
        if molmedb_col in combined_df.columns:
            print(f"Processing {Cons.MOLMEDB} column for compounds: {molmedb_col}")

            inchikeys = []

            for entry in combined_df[molmedb_col].dropna():
                if isinstance(entry, list):
                    for compound in entry:
                        if isinstance(compound, dict) and Cons.MOLMEDB_INCHIKEY in compound:
                            inchikey = compound[Cons.MOLMEDB_INCHIKEY]
                            if isinstance(inchikey, str) and inchikey:
                                inchikeys.append(inchikey)

            inchikeys = list(set(inchikeys))

            if inchikeys:
                pubchem_ids, inchikey_to_pubchem = query_bridgedb_for_pubchem(
                    list(inchikeys), "InChIKey"
                )

                if pubchem_ids:
                    molmedb_map = query_compoundwiki(pubchem_ids)
                    annotation_map.update(molmedb_map)

                    combined_df = inject_compoundwiki_annotations(
                        combined_df,
                        molmedb_col,
                        Cons.MOLMEDB_INCHIKEY,
                        {
                            inchikey: molmedb_map[pcid]
                            for inchikey, pcid in inchikey_to_pubchem.items()
                            if pcid in molmedb_map
                        },
                    )

    end_time = datetime.datetime.now()

    compoundwiki_metadata = {
        "datasource": Cons.COMPOUNDWIKI,
        "metadata": compoundwiki_version,
        "query": {
            "size": sum(len(v) if v is not None else 0 for v in annotation_map.values()),
            "time": str(end_time - start_time),
            "date": end_time.strftime("%Y-%m-%d %H:%M:%S"),
            "url": Cons.COMPOUNDWIKI_ENDPOINT,
        },
    }

    return combined_df, compoundwiki_metadata


def query_bridgedb_for_pubchem(
    compound_ids: List[str], input_datatype: str  # type: ignore
) -> Tuple[List[str], dict]:
    """Query BridgeDb to convert compound identifiers to PubChem Compound IDs.

    :param compound_ids: List of original identifiers (e.g., ChEBI, KEGG)
    :param input_datatype: Data source of input IDs (e.g., 'ChEBI', 'KEGG Compound')
    :returns: A tuple of (list of PubChem IDs, mapping from original ID to PubChem ID)
    """
    data_input = pd.DataFrame(compound_ids, columns=["identifier"])

    bridgedb_df, bridgedb_metadata = id_mapper.bridgedb_xref(
        identifiers=data_input,
        input_species="Human",
        input_datasource=input_datatype,  # type: ignore
        output_datasource=["All"],
    )

    if bridgedb_df.empty:
        return [], {}

    # Filter for rows with PubChem Compound in target.source
    pubchem_df = bridgedb_df[bridgedb_df["target.source"] == Cons.COMPOUNDWIKI_COMPOUND_INPUT_ID]

    if pubchem_df.empty:
        return [], {}

    pubchem_ids = list(set(pubchem_df[Cons.TARGET_COL].tolist()))
    id_to_pubchem = dict(pubchem_df[["identifier", Cons.TARGET_COL]].values)

    return pubchem_ids, id_to_pubchem
