# coding: utf-8

"""Python file for querying the OpenTargets database (https://www.opentargets.org/)."""

import datetime
import warnings
from typing import Tuple

import pandas as pd
import requests

from pyBiodatafuse import id_mapper
from pyBiodatafuse.constants import (
    OPENTARGETS,
    OPENTARGETS_COMPOUND_COL,
    OPENTARGETS_COMPOUND_INPUT_ID,
    OPENTARGETS_COMPOUND_OUTPUT_DICT,
    OPENTARGETS_DISEASE_COL,
    OPENTARGETS_DISEASE_OUTPUT_DICT,
    OPENTARGETS_ENDPOINT,
    OPENTARGETS_GENE_INPUT_ID,
    OPENTARGETS_GO_COL,
    OPENTARGETS_GO_OUTPUT_DICT,
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
        return pd.DataFrame(), opentargets_version

    intermediate_df.rename(
        columns={
            "id": "go_id",
            "name": "go_name",
        },
        inplace=True,
    )

    print(intermediate_df)
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

    """Metdata details"""
    # Get the current date and time
    current_date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    # Calculate the time elapsed
    time_elapsed = str(end_time - start_time)
    # Calculate the number of new nodes
    num_new_nodes = intermediate_df["go_id"].nunique()
    # Calculate the number of new edges
    num_edges = len(intermediate_df)

    # Add version, datasource, query, query time, and the date to metadata
    opentargets_version["query"] = {
        "size": len(gene_ids),
        "input_type": OPENTARGETS_GENE_INPUT_ID,
        "number_of_added_nodes": num_new_nodes,
        "number_of_added_edges": num_edges,
        "time": time_elapsed,
        "date": current_date,
        "url": OPENTARGETS_ENDPOINT,
    }

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

    # Generate the OpenTargets DataFrame
    intermediate_df = pd.DataFrame()

    for gene in r["data"]["targets"]:
        path_df = pd.DataFrame(gene["pathways"])
        path_df = path_df.drop_duplicates()
        path_df["target"] = gene["id"]
        intermediate_df = pd.concat([intermediate_df, path_df], ignore_index=True)

    if intermediate_df.empty:
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

    return merged_df, opentargets_version


def get_compound_disease_interactions(
    data_df: pd.DataFrame,
) -> Tuple[pd.DataFrame, dict]:
    """Get information about drugs associated with diseases of interest.

    :param data_df: Dataframe with list of ChEMBL ids to query.
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

    chembl_ids = data_df["target"].tolist()

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

    intermediate_df = pd.DataFrame()

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
        disease_df.drop(
            columns=["therapeuticAreas", "dbXRefs"], inplace=True
        )  # Xrefs has all other ids for diseases

        disease_df["target"] = drug["id"]
        disease_df["drug_name"] = drug["name"]
        disease_df["max_clinical_trial_phase"] = drug["maximumClinicalTrialPhase"]
        disease_df["is_withdrawn"] = drug["hasBeenWithdrawn"]

        intermediate_df = pd.concat([intermediate_df, disease_df], ignore_index=True)

    if intermediate_df.empty:
        return pd.DataFrame(), opentargets_version

    # Check if all keys in df match the keys in OUTPUT_DICT
    check_columns_against_constants(
        data_df=intermediate_df,
        output_dict=OPENTARGETS_DISEASE_OUTPUT_DICT,
        check_values_in=["disease_id", "drug_id"],
    )

    # Merge the two DataFrames on the target column
    merged_df = collapse_data_sources(
        data_df=data_df,
        source_namespace=OPENTARGETS_COMPOUND_INPUT_ID,
        target_df=intermediate_df,
        common_cols=["target"],
        target_specific_cols=list(OPENTARGETS_DISEASE_OUTPUT_DICT.keys()),
        col_name=OPENTARGETS_DISEASE_COL,
    )

    """Metdata details"""
    # Get the current date and time
    current_date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    # Calculate the time elapsed
    time_elapsed = str(end_time - start_time)
    # Add version, datasource, query, query time, and the date to metadata
    opentargets_version["query"] = {
        "size": len(chembl_ids),
        "input_type": OPENTARGETS_COMPOUND_INPUT_ID,
        "time": time_elapsed,
        "date": current_date,
        "url": OPENTARGETS_ENDPOINT,
    }

    return merged_df, opentargets_version
