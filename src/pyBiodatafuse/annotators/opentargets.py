# coding: utf-8

"""Python file for querying the OpenTargets database (https://www.opentargets.org/)."""

import datetime
import warnings
from typing import List, Tuple

import numpy as np
import pandas as pd
import requests

from pyBiodatafuse import id_mapper
from pyBiodatafuse.constants import (
    OPENTARGETS,
    OPENTARGETS_COMPOUND_COL,
    OPENTARGETS_COMPOUND_INPUT_ID,
    OPENTARGETS_COMPOUND_OUTPUT_DICT,
    OPENTARGETS_COMPOUND_QUERY_INPUT_ID,
    OPENTARGETS_DISEASE_COL,
    OPENTARGETS_DISEASE_OUTPUT_DICT,
    OPENTARGETS_ENDPOINT,
    OPENTARGETS_GENE_INPUT_ID,
    OPENTARGETS_GO_COL,
    OPENTARGETS_GO_OUTPUT_DICT,
    OPENTARGETS_IGNORE_DISEASE_IDS,
    OPENTARGETS_REACTOME_COL,
    OPENTARGETS_REACTOME_OUTPUT_DICT,
)
from pyBiodatafuse.utils import (
    check_columns_against_constants,
    collapse_data_sources,
    get_identifier_of_interest,
)


def check_endpoint_opentargets() -> bool:
    """Check the availability of the OpenTargets API endpoint.

    :returns: a dictionary containing the version information
    """
    query = """
        query MetaInfo {
            meta{
                name
                apiVersion{
                        x
                        y
                        z
                }
                dataVersion{
                        year
                        month
                }
            }
        }"""
    r = requests.post(OPENTARGETS_ENDPOINT, json={"query": query}).json()

    if not r["data"]:
        return False

    return True


def get_version_opentargets() -> dict:
    """Get version of OpenTargets API.

    :returns: a dictionary containing the version information
    """
    query = """
        query MetaInfo {
            meta{
                name
                apiVersion{
                        x
                        y
                        z
                }
                dataVersion{
                        year
                        month
                }
            }
        }"""
    r = requests.post(OPENTARGETS_ENDPOINT, json={"query": query}).json()

    metadata = {
        "datasource": r["data"]["meta"]["name"],
        "metadata": {
            "source_version": {
                "apiVersion": {
                    "x": r["data"]["meta"]["apiVersion"]["x"],
                    "y": r["data"]["meta"]["apiVersion"]["y"],
                    "z": r["data"]["meta"]["apiVersion"]["z"],
                }
            },
            "data_version": {
                "dataVersion": {
                    "year": r["data"]["meta"]["dataVersion"]["year"],
                    "month": r["data"]["meta"]["dataVersion"]["month"],
                }
            },
        },
    }

    return metadata


def get_gene_go_process(
    bridgedb_df: pd.DataFrame,
) -> Tuple[pd.DataFrame, dict]:
    """Get information about GO pathways associated with a genes of interest.

    :param bridgedb_df: BridgeDb output for creating the list of gene ids to query
    :returns: a DataFrame containing the OpenTargets output and dictionary of the query metadata.
    """
    # Check if the API is available
    api_available = check_endpoint_opentargets()
    if not api_available:
        warnings.warn(
            f"{OPENTARGETS} GraphQL endpoint is not available. Unable to retrieve data.",
            stacklevel=2,
        )
        return pd.DataFrame(), {}

    data_df = get_identifier_of_interest(bridgedb_df, OPENTARGETS_GENE_INPUT_ID)
    gene_ids = data_df["target"].tolist()

    # Record the start time
    opentargets_version = get_version_opentargets()
    start_time = datetime.datetime.now()

    query_string = """
      query targetPathways {
        targets (ensemblIds: $ids){
          id
          geneOntology {
            term {
              id
              name
            }
            aspect
          }
        }
      }
    """
    query_string = query_string.replace("$ids", str(gene_ids).replace("'", '"'))

    r = requests.post(OPENTARGETS_ENDPOINT, json={"query": query_string}).json()

    # Record the end time
    end_time = datetime.datetime.now()

    """Metdata details"""
    # Get the current date and time
    current_date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    # Calculate the time elapsed
    time_elapsed = str(end_time - start_time)

    # Add version, datasource, query, query time, and the date to metadata
    opentargets_version["query"] = {
        "size": len(gene_ids),
        "input_type": OPENTARGETS_GENE_INPUT_ID,
        "time": time_elapsed,
        "date": current_date,
        "url": OPENTARGETS_ENDPOINT,
    }

    # Generate the OpenTargets DataFrame
    intermediate_df = pd.DataFrame()

    for gene in r["data"]["targets"]:
        terms = [i["term"] for i in gene["geneOntology"]]
        types = [i["aspect"] for i in gene["geneOntology"]]
        path_df = pd.DataFrame(terms)
        path_df["go_type"] = types
        path_df = path_df.drop_duplicates()
        path_df["target"] = gene["id"]
        intermediate_df = pd.concat([intermediate_df, path_df], ignore_index=True)

    if intermediate_df.empty:
        warnings.warn(
            f"There is no annotation for your input list in {OPENTARGETS_GO_COL}.",
            stacklevel=2,
        )
        return pd.DataFrame(), opentargets_version

    intermediate_df.rename(
        columns={
            "id": "go_id",
            "name": "go_name",
        },
        inplace=True,
    )

    # Check if all keys in df match the keys in OUTPUT_DICT
    check_columns_against_constants(
        data_df=intermediate_df,
        output_dict=OPENTARGETS_GO_OUTPUT_DICT,
        check_values_in=["go_id"],
    )

    # Merge the two DataFrames on the target column
    merged_df = collapse_data_sources(
        data_df=data_df,
        source_namespace=OPENTARGETS_GENE_INPUT_ID,
        target_df=intermediate_df,
        common_cols=["target"],
        target_specific_cols=list(OPENTARGETS_GO_OUTPUT_DICT.keys()),
        col_name=OPENTARGETS_GO_COL,
    )

    """Update metadata"""
    # Calculate the number of new nodes
    num_new_nodes = intermediate_df["go_id"].nunique()
    # Calculate the number of new edges
    num_new_edges = intermediate_df.drop_duplicates(subset=["target", "go_id"]).shape[0]

    # Check the intermediate_df
    if num_new_edges != len(intermediate_df):
        warnings.warn(
            f"The intermediate_df in {OPENTARGETS_GO_COL} annotatur should be checked, please create an issue on https://github.com/BioDataFuse/pyBiodatafuse/issues/.",
            stacklevel=2,
        )

    # Add the number of new nodes and edges to metadata
    opentargets_version["query"]["number_of_added_nodes"] = num_new_nodes
    opentargets_version["query"]["number_of_added_edges"] = num_new_edges

    return merged_df, opentargets_version


def get_gene_reactome_pathways(
    bridgedb_df: pd.DataFrame,
) -> Tuple[pd.DataFrame, dict]:
    """Get information about Reactome pathways associated with a gene.

    :param bridgedb_df: BridgeDb output for creating the list of gene ids to query
    :returns: a DataFrame containing the OpenTargets output and dictionary of the query metadata.
    """
    # Check if the API is available
    api_available = check_endpoint_opentargets()
    if not api_available:
        warnings.warn(
            f"{OPENTARGETS} GraphQL endpoint is not available. Unable to retrieve data.",
            stacklevel=2,
        )
        return pd.DataFrame(), {}

    data_df = get_identifier_of_interest(bridgedb_df, OPENTARGETS_GENE_INPUT_ID)
    gene_ids = data_df["target"].tolist()

    # Record the start time
    opentargets_version = get_version_opentargets()
    start_time = datetime.datetime.now()

    query_string = """
      query targetPathways {
        targets (ensemblIds: $ids){
          id
          pathways {
            pathway
            pathwayId
          }
        }
      }
    """
    query_string = query_string.replace("$ids", str(gene_ids).replace("'", '"'))

    r = requests.post(OPENTARGETS_ENDPOINT, json={"query": query_string}).json()

    # Record the end time
    end_time = datetime.datetime.now()
    
    """Metdata details"""
    # Get the current date and time
    current_date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    # Calculate the time elapsed
    time_elapsed = str(end_time - start_time)

    # Add version, datasource, query, query time, and the date to metadata
    opentargets_version["query"] = {
        "size": len(gene_ids),
        "input_type": OPENTARGETS_GENE_INPUT_ID,
        "time": time_elapsed,
        "date": current_date,
        "url": OPENTARGETS_ENDPOINT,
    }

    # Generate the OpenTargets DataFrame
    intermediate_df = pd.DataFrame()

    for gene in r["data"]["targets"]:
        path_df = pd.DataFrame(gene["pathways"])
        path_df = path_df.drop_duplicates()
        path_df["target"] = gene["id"]
        intermediate_df = pd.concat([intermediate_df, path_df], ignore_index=True)

    if intermediate_df.empty:
        warnings.warn(
            f"There is no annotation for your input list in {OPENTARGETS_REACTOME_COL}.",
            stacklevel=2,
        )
        return pd.DataFrame(), opentargets_version

    intermediate_df.rename(
        columns={"pathway": "pathway_label", "pathwayId": "pathway_id"}, inplace=True
    )

    # Check if all keys in df match the keys in OUTPUT_DICT
    check_columns_against_constants(
        data_df=intermediate_df,
        output_dict=OPENTARGETS_REACTOME_OUTPUT_DICT,
        check_values_in=["pathway_id"],
    )

    # Merge the two DataFrames on the target column
    merged_df = collapse_data_sources(
        data_df=data_df,
        source_namespace=OPENTARGETS_GENE_INPUT_ID,
        target_df=intermediate_df,
        common_cols=["target"],
        target_specific_cols=["pathway_label", "pathway_id"],
        col_name=OPENTARGETS_REACTOME_COL,
    )

    """Update metadata"""
    # Calculate the number of new nodes
    num_new_nodes = intermediate_df["pathway_id"].nunique()
    # Calculate the number of new edges
    num_new_edges = intermediate_df.drop_duplicates(subset=["target", "pathway_id"]).shape[0]

    # Check the intermediate_df
    if num_new_edges != len(intermediate_df):
        warnings.warn(
            f"The intermediate_df in {OPENTARGETS_REACTOME_COL} annotatur should be checked, please create an issue on https://github.com/BioDataFuse/pyBiodatafuse/issues/.",
            stacklevel=2,
        )
        
    # Add the number of new nodes and edges to metadata
    opentargets_version["query"]["number_of_added_nodes"] = num_new_nodes
    opentargets_version["query"]["number_of_added_edges"] = num_new_edges

    return merged_df, opentargets_version


# TODO: Look into the utility of this function while applying filters
def get_gene_tractability(
    bridgedb_df: pd.DataFrame,
) -> Tuple[pd.DataFrame, dict]:
    """Get tractability information about a gene.

    :param bridgedb_df: BridgeDb output for creating the list of gene ids to query
    :returns: a DataFrame containing the OpenTargets output and dictionary of the query metadata.

    More metadata here- https://platform-docs.opentargets.org/target/tractability
    """
    # Check if the API is available
    api_available = check_endpoint_opentargets()
    if not api_available:
        warnings.warn(
            f"{OPENTARGETS} GraphQL endpoint is not available. Unable to retrieve data.",
            stacklevel=2,
        )
        return pd.DataFrame(), {}

    query_string = """
      query targetTractability {
        targets (ensemblIds: $ids){
          id
          tractability {
            label
            modality
            value
          }
        }
      }
    """

    data_df = get_identifier_of_interest(bridgedb_df, OPENTARGETS_GENE_INPUT_ID)
    gene_ids = data_df["target"].tolist()

    # Record the start time
    opentargets_version = get_version_opentargets()
    start_time = datetime.datetime.now()

    query_string = query_string.replace("$ids", str(gene_ids).replace("'", '"'))

    r = requests.post(OPENTARGETS_ENDPOINT, json={"query": query_string}).json()

    # Record the end time
    end_time = datetime.datetime.now()

    intermediate_df = pd.DataFrame()

    for gene in r["data"]["targets"]:
        tract_df = pd.DataFrame(gene["tractability"])
        tract_df = tract_df.drop_duplicates()
        tract_df["ensembl_id"] = gene["id"]
        intermediate_df = pd.concat([intermediate_df, tract_df], ignore_index=True)

    if intermediate_df.empty:
        return pd.DataFrame(), opentargets_version

    intermediate_df = intermediate_df[intermediate_df["value"] == True]
    intermediate_df.drop(columns=["value"], inplace=True)

    # TODO: Check if all keys in df match the keys in OUTPUT_DICT
    # TODO: Merge data wih main data_df

    """Metdata details"""
    # Get the current date and time
    current_date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    # Calculate the time elapsed
    time_elapsed = str(end_time - start_time)
    # Add version, datasource, query, query time, and the date to metadata
    opentargets_version["query"] = {
        "size": len(gene_ids),
        "input_type": OPENTARGETS_GENE_INPUT_ID,
        "time": time_elapsed,
        "date": current_date,
        "url": OPENTARGETS_ENDPOINT,
    }

    return intermediate_df, opentargets_version


def get_gene_compound_interactions(
    bridgedb_df: pd.DataFrame,
) -> Tuple[pd.DataFrame, dict]:
    """Get information about drugs associated with a genes of interest.

    :param bridgedb_df: BridgeDb output for creating the list of gene ids to query
    :returns: a DataFrame containing the OpenTargets output and dictionary of the query metadata.
    """
    # Check if the API is available
    api_available = check_endpoint_opentargets()
    if not api_available:
        warnings.warn(
            f"{OPENTARGETS} GraphQL endpoint is not available. Unable to retrieve data.",
            stacklevel=2,
        )
        return pd.DataFrame(), {}

    data_df = get_identifier_of_interest(bridgedb_df, OPENTARGETS_GENE_INPUT_ID)
    gene_ids = data_df["target"].tolist()

    # Record the start time
    opentargets_version = get_version_opentargets()
    start_time = datetime.datetime.now()

    query_string = """
      query targetDrugs {
        targets (ensemblIds: $ids){
          id
          knownDrugs {
            rows {
              mechanismOfAction
              drug {
                id
                name
                isApproved
                crossReferences {
                    source
                    reference
                }
                adverseEvents{
                    count
                    rows{
                      name
                    }
                  }
              }
            }
          }
        }
      }
    """
    query_string = query_string.replace("$ids", str(gene_ids).replace("'", '"'))

    r = requests.post(OPENTARGETS_ENDPOINT, json={"query": query_string}).json()

    # Record the end time
    end_time = datetime.datetime.now()

    """Metdata details"""
    # Get the current date and time
    current_date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    # Calculate the time elapsed
    time_elapsed = str(end_time - start_time)

    # Add version, datasource, query, query time, and the date to metadata
    opentargets_version["query"] = {
        "size": len(gene_ids),
        "input_type": OPENTARGETS_GENE_INPUT_ID,
        "time": time_elapsed,
        "date": current_date,
        "url": OPENTARGETS_ENDPOINT,
    }

    # Generate the OpenTargets DataFrame
    intermediate_df = pd.DataFrame()

    for gene in r["data"]["targets"]:
        if not gene["knownDrugs"]:
            continue

        drug_info = gene["knownDrugs"]["rows"]
        drug_df = pd.DataFrame(drug_info)

        if drug_df.empty:
            continue

        drug_df[
            ["chembl_id", "compound_name", "is_approved", "cross_references", "adverse_events"]
        ] = drug_df["drug"].apply(pd.Series)
        drug_df.drop(columns=["drug"], inplace=True)

        drug_df["target"] = gene["id"]

        drug_df["mechanismOfAction"] = drug_df["mechanismOfAction"].apply(
            lambda x: "inhibits" if "antagonist" in x else "activates"
        )
        drug_df.rename(columns={"mechanismOfAction": "relation"}, inplace=True)

        drug_df["drugbank_id"] = drug_df["cross_references"].apply(
            lambda x: (
                next((ref["reference"][0] for ref in x if ref["source"] == "drugbank"), None)
                if x
                else None
            )
        )

        drug_df[["adverse_effect_count", "adverse_effect"]] = drug_df.apply(
            lambda row: (
                pd.Series([row["adverse_events"]["count"], row["adverse_events"]["rows"]])
                if row["adverse_events"]
                else pd.Series([None, None])
            ),
            axis=1,
        )

        intermediate_df = pd.concat([intermediate_df, drug_df], ignore_index=True)
        intermediate_df.drop(["cross_references", "adverse_events"], axis=1, inplace=True)
        intermediate_df = intermediate_df.drop_duplicates(
            subset=[col for col in intermediate_df.columns if col != "adverse_effect"]
        )

    if intermediate_df.empty:
        warnings.warn(
            f"There is no annotation for your input list in {OPENTARGETS_COMPOUND_COL}.",
            stacklevel=2,
        )
        return pd.DataFrame(), opentargets_version

    # Fixing chembl_id to pubchem_id
    mapped_df, _ = id_mapper.pubchem_xref(
        identifiers=intermediate_df["chembl_id"], identifier_type="name"
    )
    intermediate_df["compound_cid"] = mapped_df["target"]

    # Check if all keys in df match the keys in OUTPUT_DICT
    check_columns_against_constants(
        data_df=intermediate_df,
        output_dict=OPENTARGETS_COMPOUND_OUTPUT_DICT,
        check_values_in=["chembl_id", "drugbank_id", "relation"],
    )

    # Merge the two DataFrames on the target column
    merged_df = collapse_data_sources(
        data_df=data_df,
        source_namespace=OPENTARGETS_GENE_INPUT_ID,
        target_df=intermediate_df,
        common_cols=["target"],
        target_specific_cols=list(OPENTARGETS_COMPOUND_OUTPUT_DICT.keys()),
        col_name=OPENTARGETS_COMPOUND_COL,
    )

    """Update metadata"""
    # Calculate the number of new nodes
    num_new_nodes = intermediate_df["chembl_id"].nunique()
    # Calculate the number of new edges
    num_new_edges = intermediate_df.drop_duplicates(subset=["target", "chembl_id"]).shape[0]

    # Check the intermediate_df
    if num_new_edges != len(intermediate_df):
        warnings.warn(
            f"The intermediate_df in {OPENTARGETS_COMPOUND_COL} annotatur should be checked, please create an issue on https://github.com/BioDataFuse/pyBiodatafuse/issues/.",
            stacklevel=2,
        )

    # Add the number of new nodes and edges to metadata
    opentargets_version["query"]["number_of_added_nodes"] = num_new_nodes
    opentargets_version["query"]["number_of_added_edges"] = num_new_edges

    return merged_df, opentargets_version


def _process_disease_xref(row) -> List[str]:
    tmp = []
    for val in row:
        namespace, idx = val.split(":")
        if namespace.lower() in OPENTARGETS_IGNORE_DISEASE_IDS:  # skipping non essential ones
            continue
        if namespace.lower() in ["mesh", "msh"]:
            tmp.append(f"MESH_{idx}")
        elif namespace.lower() in ["ncit"]:
            tmp.append(f"NCI_{idx}")
        elif namespace.lower() in ["omim"]:
            tmp.append(f"OMIM_{idx}")
        elif namespace.lower() in ["mondo"]:
            tmp.append(f"MONDO_{idx}")
        elif namespace.lower() in ["efo"]:
            tmp.append(f"EFO_{idx}")
        elif namespace.lower() in ["doid"]:
            tmp.append(f"DO_{idx}")
        elif namespace.lower() in ["umls"]:
            tmp.append(f"UMLS_{idx}")
        elif namespace.lower() in ["hp"]:
            tmp.append(f"HPO_HP:{idx}")
        elif namespace.lower() in ["orphanet", "ordo"]:
            tmp.append(f"ORDO_{idx}")
        else:
            raise ValueError(f"Unknown namespace: {namespace}")
    return tmp


def get_compound_disease_interactions(
    bridgedb_df: pd.DataFrame,
) -> Tuple[pd.DataFrame, dict]:
    """Get information about drugs associated with diseases of interest.

    :param bridgedb_df: BridgeDb output for creating the list of gene ids to query
    :returns: a DataFrame containing the OpenTargets output and dictionary of the query metadata.
    """
    # Check if the API is available
    api_available = check_endpoint_opentargets()
    if not api_available:
        warnings.warn(
            f"{OPENTARGETS} GraphQL endpoint is not available. Unable to retrieve data.",
            stacklevel=2,
        )
        return pd.DataFrame(), {}

    if (
        OPENTARGETS_COMPOUND_QUERY_INPUT_ID in bridgedb_df["target.source"].values
    ):  # for chembl_id in col
        data_df = bridgedb_df[bridgedb_df["target.source"] == OPENTARGETS_COMPOUND_QUERY_INPUT_ID]
        chembl_gene_map = data_df.set_index("target")[
            "identifier"
        ].to_dict()  # Dict of chembl_id:gene_id
        chembl_cid_map = None
        chembl_ids = set(data_df["target"].tolist())
    else:
        data_df = get_identifier_of_interest(bridgedb_df, OPENTARGETS_COMPOUND_INPUT_ID)
        pubchem_ids = data_df["target"].tolist()
        chembl_gene_map = None
        chembl_cid_map = id_mapper.cid2chembl(pubchem_ids)  # Dict of chembl_id:pubchem_id
        chembl_ids = set(list(chembl_cid_map.keys()))

    # remove nan entries
    chembl_ids = [x for x in chembl_ids if str(x) != "nan"]  # type: ignore

    # Record the start time
    opentargets_version = get_version_opentargets()
    start_time = datetime.datetime.now()

    query_string = """
    query KnownDrugsQuery{
        drugs (chemblIds: $chemblIds) {
            id
            name
            maximumClinicalTrialPhase
            hasBeenWithdrawn
            linkedDiseases {
                rows {
                    id
                    name
                    dbXRefs
                    therapeuticAreas {
                        id
                        name
                    }
                }
            }
        }
    }"""
    query_string = query_string.replace("$chemblIds", str(chembl_ids).replace("'", '"'))

    r = requests.post(OPENTARGETS_ENDPOINT, json={"query": query_string}).json()
    # Record the end time
    end_time = datetime.datetime.now()
    
    """Metdata details"""
    # Get the current date and time
    current_date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    # Calculate the time elapsed
    time_elapsed = str(end_time - start_time)

    # Add version, datasource, query, query time, and the date to metadata
    opentargets_version["query"] = {
        "size": len(chembl_ids),
        "time": time_elapsed,
        "date": current_date,
        "url": OPENTARGETS_ENDPOINT,
    }
    
    # Generate the OpenTargets DataFrame
    intermediate_df = pd.DataFrame()

    if r["data"]["drugs"] is None:
        warnings.warn(
            f"There is no annotation for your input list in {OPENTARGETS_DISEASE_COL}.",
            stacklevel=2,
        )
        return pd.DataFrame(), opentargets_version

    for drug in r["data"]["drugs"]:
        if not drug["linkedDiseases"]:
            continue

        # Based on clinical trial data
        disease_info = drug["linkedDiseases"]["rows"]
        disease_df = pd.DataFrame(disease_info)

        if disease_df.empty:
            continue

        disease_df["therapeutic_areas"] = disease_df["therapeuticAreas"].apply(
            lambda x: ", ".join([f"{i['id']}:{i['name']}" for i in x])
        )  # Multiple values separated by comma

        disease_df.rename(
            columns={"id": "opentarget_disease_id", "name": "disease_name"}, inplace=True
        )

        disease_df["disease_id"] = disease_df.apply(
            lambda row: ", ".join(
                ["umls:" + i.split(":")[1] for i in row["dbXRefs"] if i.startswith("UMLS:")]
                or [row["opentarget_disease_id"]]
            ),
            axis=1,
        )

        # Fixing the xrefs
        xrefs = []  # type: ignore

        for row in disease_df["dbXRefs"]:
            if len(row) == 0:
                xrefs.append([])
                continue
            tmp = _process_disease_xref(row)
            xrefs.append(tmp)
        disease_df["disease_xrefs"] = xrefs

        if chembl_cid_map:
            disease_df["target"] = chembl_cid_map.get(drug["id"], None)
        else:
            disease_df["identifier"] = chembl_gene_map[drug["id"]]

        disease_df["drug_name"] = drug["name"]
        disease_df["max_clinical_trial_phase"] = drug["maximumClinicalTrialPhase"]
        disease_df["is_withdrawn"] = drug["hasBeenWithdrawn"]

        disease_df.drop(columns=["therapeuticAreas", "dbXRefs"], inplace=True)

        intermediate_df = pd.concat([intermediate_df, disease_df], ignore_index=True)

    if intermediate_df.empty:
        warnings.warn(
            f"There is no annotation for your input list in {OPENTARGETS_DISEASE_COL}.",
            stacklevel=2,
        )
        return pd.DataFrame(), opentargets_version

    # Check if all keys in df match the keys in OUTPUT_DICT
    check_columns_against_constants(
        data_df=intermediate_df,
        output_dict=OPENTARGETS_DISEASE_OUTPUT_DICT,
        check_values_in=["disease_id", "drug_id"],
    )

    # Merge the two DataFrames on the target column
    if (
        OPENTARGETS_COMPOUND_QUERY_INPUT_ID in bridgedb_df["target.source"].values
    ):  # for chembl_id in col (with gene-compoud function)
        merged_df = collapse_data_sources(
            data_df=bridgedb_df,
            source_namespace=OPENTARGETS_COMPOUND_QUERY_INPUT_ID,
            target_df=intermediate_df,
            common_cols=["identifier"],
            target_specific_cols=list(OPENTARGETS_DISEASE_OUTPUT_DICT.keys()),
            col_name=OPENTARGETS_DISEASE_COL,
        )
        opentargets_version["query"]["input_type"] = OPENTARGETS_COMPOUND_QUERY_INPUT_ID

    else:
        merged_df = collapse_data_sources(
            data_df=bridgedb_df,
            source_namespace=OPENTARGETS_COMPOUND_INPUT_ID,
            target_df=intermediate_df,
            common_cols=["target"],
            target_specific_cols=list(OPENTARGETS_DISEASE_OUTPUT_DICT.keys()),
            col_name=OPENTARGETS_DISEASE_COL,
        )
        opentargets_version["query"]["input_type"] = OPENTARGETS_COMPOUND_INPUT_ID

    # Calculate the number of new nodes
    num_new_nodes = intermediate_df["opentarget_disease_id"].nunique()
    # Calculate the number of new edges
    num_new_edges = intermediate_df.drop_duplicates(subset=["target", "opentarget_disease_id"]).shape[0]

    # Check the intermediate_df
    if num_new_edges != len(intermediate_df):
        warnings.warn(
            f"The intermediate_df in {OPENTARGETS_DISEASE_COL} annotatur should be checked, please create an issue on https://github.com/BioDataFuse/pyBiodatafuse/issues/.",
            stacklevel=2,
        )

    # Add the number of new nodes and edges to metadata
    opentargets_version["query"]["number_of_added_nodes"] = num_new_nodes
    opentargets_version["query"]["number_of_added_edges"] = num_new_edges

    return merged_df, opentargets_version


def get_gene_disease_interactions(
    bridgedb_df: pd.DataFrame,
) -> Tuple[pd.DataFrame, dict]:
    """Get information about disease connected to drugs associated with a genes of interest.

    :param bridgedb_df: BridgeDb output for creating the list of gene ids to query
    :returns: a DataFrame containing the OpenTargets output and dictionary of the query metadata.
    """
    # Check if the API is available
    api_available = check_endpoint_opentargets()
    if not api_available:
        warnings.warn(
            f"{OPENTARGETS} GraphQL endpoint is not available. Unable to retrieve data.",
            stacklevel=2,
        )
        return pd.DataFrame(), {}

    # Get gene-compound interactions
    cmp_merged_df, opentargets_version = get_gene_compound_interactions(bridgedb_df)

    if cmp_merged_df.empty:
        return pd.DataFrame(), opentargets_version

    # Making a bridgeDb dataframe for the compound ids
    gene_cmpd_data = []
    for row in cmp_merged_df.values:
        (
            gene,
            gene_namespace,
            tar,
            target_namespace,
            cmpds,
        ) = row
        if len(cmpds) < 1:
            continue

        for cmpd in cmpds:
            gene_cmpd_data.append(
                {
                    "identifier": gene,
                    "identifier.source": gene_namespace,
                    "target": cmpd["chembl_id"],
                    "target.source": OPENTARGETS_COMPOUND_QUERY_INPUT_ID,
                }
            )

    gene_cmpd_df = pd.DataFrame(gene_cmpd_data)

    # Get compound-disease interactions
    dis_merged_df, opentargets_version = get_compound_disease_interactions(gene_cmpd_df)

    # Fixing merge to look like gene-disease
    merged_df = bridgedb_df.copy(deep=True)

    disease_col = []
    for gene in merged_df["identifier"]:
        tmp = dis_merged_df[dis_merged_df["identifier"] == gene]
        if tmp.empty:
            disease_col.append([{i: np.nan for i in OPENTARGETS_DISEASE_OUTPUT_DICT}])
        else:
            vals = tmp.iloc[0][OPENTARGETS_DISEASE_COL]
            disease_col.append(vals)

    merged_df[OPENTARGETS_DISEASE_COL] = disease_col

    return merged_df, opentargets_version
