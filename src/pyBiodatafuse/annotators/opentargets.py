# coding: utf-8

"""Python file for querying the OpenTargets database (https://www.opentargets.org/)."""

import datetime
import math
import warnings
from typing import Tuple

import pandas as pd
import requests

from pyBiodatafuse.constants import (
    OPENTARGETS,
    OPENTARGETS_COMPOUND_COL,
    OPENTARGETS_COMPOUND_OUTPUT_DICT,
    OPENTARGETS_DISEASE_COL,
    OPENTARGETS_DISEASE_OUTPUT_DICT,
    OPENTARGETS_ENDPOINT,
    OPENTARGETS_GO_COL,
    OPENTARGETS_GO_OUTPUT_DICT,
    OPENTARGETS_INPUT_ID,
    OPENTARGETS_LOCATION_COL,
    OPENTARGETS_LOCATION_OUTPUT_DICT,
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


def get_gene_location(
    bridgedb_df: pd.DataFrame,
) -> Tuple[pd.DataFrame, dict]:
    """Get location of gene in human body.

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

    data_df = get_identifier_of_interest(bridgedb_df, OPENTARGETS_INPUT_ID)
    gene_ids = data_df["target"].tolist()

    # Record the start time
    opentargets_version = get_version_opentargets()
    start_time = datetime.datetime.now()

    query_string = """
      query targetLocation {
        targets (ensemblIds: $ids){
          id
          subcellularLocations {
            location
            termSL
            labelSL
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
        gene_id = gene["id"]
        loc_df = pd.DataFrame(gene["subcellularLocations"])
        loc_df = loc_df.drop_duplicates()
        loc_df["target"] = gene_id
        intermediate_df = pd.concat([intermediate_df, loc_df], ignore_index=True)

    if intermediate_df.empty:  # If no data is returned, return empty DataFrame
        return pd.DataFrame(), opentargets_version

    intermediate_df.dropna(
        subset=["termSL"], inplace=True
    )  # Drop rows where termSL is not available
    intermediate_df.rename(
        columns={
            "termSL": "location_id",
            "labelSL": "subcellular_location",
        },
        inplace=True,
    )

    # Check if all keys in df match the keys in OUTPUT_DICT
    check_columns_against_constants(
        data_df=intermediate_df,
        output_dict=OPENTARGETS_LOCATION_OUTPUT_DICT,
        check_values_in=["location_id"],
    )

    # Merge the two DataFrames on the target column
    merged_df = collapse_data_sources(
        data_df=data_df,
        source_namespace=OPENTARGETS_INPUT_ID,
        target_df=intermediate_df,
        common_cols=["target"],
        target_specific_cols=list(OPENTARGETS_LOCATION_OUTPUT_DICT.keys()),
        col_name=OPENTARGETS_LOCATION_COL,
    )

    """Metdata details"""
    # Get the current date and time
    current_date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    # Calculate the time elapsed
    time_elapsed = str(end_time - start_time)
    # Add version, datasource, query, query time, and the date to metadata
    opentargets_version["query"] = {
        "size": len(gene_ids),
        "input_type": OPENTARGETS_INPUT_ID,
        "time": time_elapsed,
        "date": current_date,
        "url": OPENTARGETS_ENDPOINT,
    }

    return merged_df, opentargets_version


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

    data_df = get_identifier_of_interest(bridgedb_df, OPENTARGETS_INPUT_ID)
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

    # Check if all keys in df match the keys in OUTPUT_DICT
    check_columns_against_constants(
        data_df=intermediate_df,
        output_dict=OPENTARGETS_GO_OUTPUT_DICT,
        check_values_in=["go_id"],
    )

    # Merge the two DataFrames on the target column
    merged_df = collapse_data_sources(
        data_df=data_df,
        source_namespace=OPENTARGETS_INPUT_ID,
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
    # Add version, datasource, query, query time, and the date to metadata
    opentargets_version["query"] = {
        "size": len(gene_ids),
        "input_type": OPENTARGETS_INPUT_ID,
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

    data_df = get_identifier_of_interest(bridgedb_df, OPENTARGETS_INPUT_ID)
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
        source_namespace=OPENTARGETS_INPUT_ID,
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
        "input_type": OPENTARGETS_INPUT_ID,
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

    data_df = get_identifier_of_interest(bridgedb_df, OPENTARGETS_INPUT_ID)
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
        "input_type": OPENTARGETS_INPUT_ID,
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

    data_df = get_identifier_of_interest(bridgedb_df, OPENTARGETS_INPUT_ID)
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
                crossReferences{source,reference}
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
        drug_df["compound_cid"] = drug_df["cross_references"].apply(
            lambda x: (
                next((ref["reference"][0] for ref in x if ref["source"] == "PubChem"), None)
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

    # Check if all keys in df match the keys in OUTPUT_DICT
    check_columns_against_constants(
        data_df=intermediate_df,
        output_dict=OPENTARGETS_COMPOUND_OUTPUT_DICT,
        check_values_in=["chembl_id", "relation"],
    )

    # TODO: Covert the ChEMBL ids to Pubchem using BridgeDb
    # Merge the two DataFrames on the target column
    merged_df = collapse_data_sources(
        data_df=data_df,
        source_namespace=OPENTARGETS_INPUT_ID,
        target_df=intermediate_df,
        common_cols=["target"],
        target_specific_cols=list(OPENTARGETS_COMPOUND_OUTPUT_DICT.keys()),
        col_name=OPENTARGETS_COMPOUND_COL,  # TODO: Cross-check if correct name
    )

    """Metdata details"""
    # Get the current date and time
    current_date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    # Calculate the time elapsed
    time_elapsed = str(end_time - start_time)
    # Add version, datasource, query, query time, and the date to metadata
    opentargets_version["query"] = {
        "size": len(gene_ids),
        "input_type": OPENTARGETS_INPUT_ID,
        "time": time_elapsed,
        "date": current_date,
        "url": OPENTARGETS_ENDPOINT,
    }

    return merged_df, opentargets_version


def get_gene_disease_associations(
    bridgedb_df: pd.DataFrame,
) -> Tuple[pd.DataFrame, dict]:
    """Get information about diseases associated with genes based on OpenTargets.

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

    data_df = get_identifier_of_interest(bridgedb_df, OPENTARGETS_INPUT_ID)
    gene_ids = data_df["target"].tolist()

    # Record the start time
    opentargets_version = get_version_opentargets()
    start_time = datetime.datetime.now()

    query_string = """
      query targetDiseases {
        targets (ensemblIds: $ids){
          id
          knownDrugs {
            rows {
              disease {
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
        disease_info = gene["knownDrugs"]["rows"]
        disease_df = pd.DataFrame(disease_info)
        disease_df = disease_df.drop_duplicates()

        if disease_df.empty:
            continue

        disease_df[["mondo_id", "disease_name", "dbXRefs", "therapeutic_area"]] = disease_df[
            "disease"
        ].apply(pd.Series)

        disease_df["therapeutic_areas"] = disease_df["therapeutic_area"].apply(
            lambda x: ", ".join([f"{i['id']}:{i['name']}" for i in x])
        )

        disease_df["disease_id"] = disease_df.apply(
            lambda row: ", ".join(
                ["umls:" + i.split(":")[1] for i in row["dbXRefs"] if i.startswith("UMLS:")]
                or [row["mondo_id"]]
            ),
            axis=1,
        )

        disease_df.drop(columns=["disease", "therapeutic_area", "dbXRefs"], inplace=True)

        disease_df["target"] = gene["id"]

        intermediate_df = pd.concat([intermediate_df, disease_df], ignore_index=True)

    if intermediate_df.empty:
        return pd.DataFrame(), opentargets_version

    # Check if all keys in df match the keys in OUTPUT_DICT
    check_columns_against_constants(
        data_df=intermediate_df,
        output_dict=OPENTARGETS_DISEASE_OUTPUT_DICT,
        check_values_in=["disease_id"],
    )

    # Merge the two DataFrames on the target column
    merged_df = collapse_data_sources(
        data_df=data_df,
        source_namespace=OPENTARGETS_INPUT_ID,
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
        "size": len(gene_ids),
        "input_type": OPENTARGETS_INPUT_ID,
        "time": time_elapsed,
        "date": current_date,
        "url": OPENTARGETS_ENDPOINT,
    }

    return merged_df, opentargets_version


# TODO: Fix the function; Modify to same format
# TODO: Check if all keys in df match the keys in OUTPUT_DICT
def get_drug_disease_interactions(
    drug_df: pd.DataFrame, disgenet_result: pd.DataFrame
) -> Tuple[pd.DataFrame, dict]:
    """Get information about drugs associated with diseases of interest.

    :param drug_df: get_gene_disease_associations output for creating the list of diseases ids to query
    :param disgenet_result: DisGeNET diseases dadtaframe
    :returns: a DataFrame containing the drug-disease relationships between drug and diseases included in your graph
    """
    # Check if the API is available
    api_available = check_endpoint_opentargets()
    if not api_available:
        warnings.warn(
            "OpenTargets GraphQL endpoint is not available. Unable to retrieve data.", stacklevel=2
        )
        return pd.DataFrame(), {}

    base_url = "https://api.platform.opentargets.org/api/v4/graphql"

    # Iterate through the dictionary and remove entries with NaN values for drugs and diseases
    # Drugs
    drugs_all = dict(drug_df["ChEMBL_Drugs"])
    drugs_list = []
    for _key, value in drugs_all.items():
        for lists in value:
            drugs_list.append(lists["chembl_id"])

    drugs_ids = [
        value
        for value in set(drugs_list)
        if not (math.isnan(value) if isinstance(value, float) else False)
    ]

    # Diseases
    diseases_all = dict(disgenet_result["DisGeNET"])
    diseases_list = []
    for _key, value in diseases_all.items():
        for lists in value:
            diseases_list.append(lists["diseaseid"])

    diseases_ids = [
        value
        for value in set(diseases_list)
        if not (math.isnan(value) if isinstance(value, float) else False)
    ]

    diseases_drugs_list = []
    for x in drugs_ids:
        query_string = """
                        query KnownDrugsQuery(
              $cursor: String
              $freeTextQuery: String
              $size: Int = 10
            ) {
              drug(chemblId: $chemblId) {
                id
                name
                knownDrugs(cursor: $cursor, freeTextQuery: $freeTextQuery, size: $size) {
                  count
                  rows {
                    disease {
                      id
                      name
                      dbXRefs
                    }
                    target {
                      id
                      approvedName
                      approvedSymbol
                    }
                  }
                }
              }
            }

            """
        query_string = query_string.replace("$chemblId", '"' + x + '"')
        r = requests.post(base_url, json={"query": query_string}).json()
        diseases_drugs_list.append(r["data"])

    # Extracting and simplifying the structure
    col_names = ["identifier", "compound_name", "drug_diseases"]
    drug_disease_df = pd.DataFrame(columns=col_names)

    for entry in diseases_drugs_list:
        entries = entry["drug"]
        drugs = entries["knownDrugs"]["rows"]
        drugs_list = []
        dict_new = None  # Initialize dict_new outside the loop

        for entry in drugs:
            drugs_fixed = entry["disease"]
            other_ids = drugs_fixed["dbXRefs"]

            umls_id = None  # Initialize umls_id before the loop

            for value in other_ids:
                if "UMLS" in value:
                    result_list = value.split(":")
                    umls_id = result_list[1]
                    break  # Once we find the UMLS ID, we can exit the loop

            if umls_id is None:
                umls_id = math.nan

            del drugs_fixed["dbXRefs"]
            drugs_fixed["umls"] = umls_id

            if drugs_fixed["umls"] in diseases_ids:  # only include diseases in my graph
                drugs_list.append(drugs_fixed)

        if drugs_list:
            dict_new = {
                "identifier": entries["id"],
                "compound_name": entries["name"],
                "drug_diseases": drugs_list,
            }

        if dict_new is not None:
            dict_new_df = pd.DataFrame([dict_new])
            drug_disease_df = pd.concat([drug_disease_df, dict_new_df], ignore_index=True)

    return drug_disease_df
