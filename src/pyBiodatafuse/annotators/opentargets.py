# coding: utf-8

"""Python file for querying the OpenTargets database (https://www.opentargets.org/)."""

import datetime
import math
from typing import Tuple

import pandas as pd
import requests

from pyBiodatafuse.utils import collapse_data_sources, get_identifier_of_interest

# URL of OpenTarget's GraphQL API endpoint
base_url = "https://api.platform.opentargets.org/api/v4/graphql"


def get_version() -> dict:
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
    r = requests.post(base_url, json={"query": query}).json()

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


def get_gene_location(bridgedb_df: pd.DataFrame) -> Tuple[pd.DataFrame, dict]:
    """Get location of gene in human body.

    :param bridgedb_df: BridgeDb output for creating the list of gene ids to query
    :returns: a DataFrame containing the OpenTargets output and dictionary of the query metadata.
    """
    data_df = get_identifier_of_interest(bridgedb_df, "Ensembl")
    gene_ids = data_df["target"].tolist()

    # Record the start time
    version_metadata = get_version()
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

    r = requests.post(base_url, json={"query": query_string}).json()

    # Record the end time
    end_time = datetime.datetime.now()

    # Update the metadata file
    version_metadata["query"] = {
        "size": len(gene_ids),
        "input_type": "Ensembl",
        "time": str(end_time - start_time),
        "date": str(datetime.datetime.now()),
        "url": base_url,
    }

    # Generate the OpenTargets DataFrame
    data = []

    for gene in r["data"]["targets"]:
        gene_id = gene["id"]
        loc_df = pd.DataFrame(gene["subcellularLocations"])
        loc_df["target"] = gene_id
        data.append(loc_df)

    if len(data) == 0:  # If no data is returned, return empty DataFrame
        return pd.DataFrame(), version_metadata

    opentargets_df = pd.concat(data)
    opentargets_df.dropna(
        subset=["termSL"], inplace=True
    )  # Drop rows where termSL is not available
    opentargets_df.rename(
        columns={
            "termSL": "loc_identifier",
            "labelSL": "subcellular_loc",
        },
        inplace=True,
    )

    # Merge the two DataFrames on the target column
    merged_df = collapse_data_sources(
        data_df=data_df,
        source_namespace="Ensembl",
        target_df=opentargets_df,
        common_cols=["target"],
        target_specific_cols=["loc_identifier", "subcellular_loc", "location"],
        col_name="OpenTargets_Location",
    )

    return merged_df, version_metadata


def get_gene_go_process(bridgedb_df: pd.DataFrame) -> Tuple[pd.DataFrame, dict]:
    """Get information about GO pathways associated with a genes of interest.

    :param bridgedb_df: BridgeDb output for creating the list of gene ids to query
    :returns: a DataFrame containing the OpenTargets output and dictionary of the query metadata.
    """
    data_df = get_identifier_of_interest(bridgedb_df, "Ensembl")
    gene_ids = data_df["target"].tolist()

    # Record the start time
    version_metadata = get_version()
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
          }
        }
      }
    """
    query_string = query_string.replace("$ids", str(gene_ids).replace("'", '"'))

    r = requests.post(base_url, json={"query": query_string}).json()

    # Record the end time
    end_time = datetime.datetime.now()

    # Update the metadata file
    version_metadata["query"] = {
        "size": len(gene_ids),
        "input_type": "Ensembl",
        "time": str(end_time - start_time),
        "date": str(datetime.datetime.now()),
        "url": base_url,
    }

    # Generate the OpenTargets DataFrame
    data = []

    for gene in r["data"]["targets"]:
        terms = [i["term"] for i in gene["geneOntology"]]
        path_df = pd.DataFrame(terms)
        path_df["target"] = gene["id"]
        data.append(path_df)

    if len(data) == 0:
        return pd.DataFrame(), version_metadata

    opentargets_df = pd.concat(data)

    opentargets_df.rename(
        columns={
            "id": "go_id",
            "name": "go_name",
        },
        inplace=True,
    )

    # Merge the two DataFrames on the target column
    merged_df = collapse_data_sources(
        data_df=data_df,
        source_namespace="Ensembl",
        target_df=opentargets_df,
        common_cols=["target"],
        target_specific_cols=["go_id", "go_name"],
        col_name="GO_Process",  # TODO: Cross-check if correct name
    )

    return merged_df, version_metadata


def get_gene_reactome_pathways(bridgedb_df: pd.DataFrame) -> Tuple[pd.DataFrame, dict]:
    """Get information about Reactome pathways associated with a gene.

    :param bridgedb_df: BridgeDb output for creating the list of gene ids to query
    :returns: a DataFrame containing the OpenTargets output and dictionary of the query metadata.
    """
    data_df = get_identifier_of_interest(bridgedb_df, "Ensembl")
    gene_ids = data_df["target"].tolist()

    # Record the start time
    version_metadata = get_version()
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

    r = requests.post(base_url, json={"query": query_string}).json()

    # Record the end time
    end_time = datetime.datetime.now()

    # Update the metadata file
    version_metadata["query"] = {
        "size": len(gene_ids),
        "input_type": "Ensembl",
        "time": str(end_time - start_time),
        "date": str(datetime.datetime.now()),
        "url": base_url,
    }

    # Generate the OpenTargets DataFrame
    data = []

    for gene in r["data"]["targets"]:
        path_df = pd.DataFrame(gene["pathways"])
        path_df["target"] = gene["id"]
        data.append(path_df)

    if len(data) == 0:
        return pd.DataFrame(), version_metadata

    opentargets_df = pd.concat(data)

    opentargets_df.rename(
        columns={"pathway": "pathway_name", "pathwayId": "pathway_id"}, inplace=True
    )

    # Merge the two DataFrames on the target column
    merged_df = collapse_data_sources(
        data_df=data_df,
        source_namespace="Ensembl",
        target_df=opentargets_df,
        common_cols=["target"],
        target_specific_cols=["pathway_name", "pathway_id"],
        col_name="Reactome_Pathways",  # TODO: Cross-check if correct name
    )

    return merged_df, version_metadata


# TODO: Look into the utility of this function while applying filters
def get_gene_tractability(bridgedb_df: pd.DataFrame) -> Tuple[pd.DataFrame, dict]:
    """Get tractability information about a gene.

    :param bridgedb_df: BridgeDb output for creating the list of gene ids to query
    :returns: a DataFrame containing the OpenTargets output and dictionary of the query metadata.

    More metadata here- https://platform-docs.opentargets.org/target/tractability
    """
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

    data_df = get_identifier_of_interest(bridgedb_df, "Ensembl")
    gene_ids = data_df["target"].tolist()

    # Record the start time
    version_metadata = get_version()
    start_time = datetime.datetime.now()

    query_string = query_string.replace("$ids", str(gene_ids).replace("'", '"'))

    r = requests.post(base_url, json={"query": query_string}).json()

    # Record the end time
    end_time = datetime.datetime.now()

    # Update the metadata file
    version_metadata["query"] = {
        "size": len(gene_ids),
        "input_type": "Ensembl",
        "time": str(end_time - start_time),
        "date": str(datetime.datetime.now()),
        "url": base_url,
    }

    data = []

    for gene in r["data"]["targets"]:
        tract_df = pd.DataFrame(gene["tractability"])
        tract_df["ensembl_id"] = gene["id"]
        data.append(tract_df)

    if len(data) == 0:
        return pd.DataFrame(), version_metadata

    opentargets_df = pd.concat(data)

    opentargets_df = opentargets_df[opentargets_df["value"] == True]
    opentargets_df.drop(columns=["value"], inplace=True)

    return opentargets_df, version_metadata


def get_gene_drug_interactions(bridgedb_df: pd.DataFrame) -> Tuple[pd.DataFrame, dict]:
    """Get information about drugs associated with a genes of interest.

    :param bridgedb_df: BridgeDb output for creating the list of gene ids to query
    :returns: a DataFrame containing the OpenTargets output and dictionary of the query metadata.
    """
    data_df = get_identifier_of_interest(bridgedb_df, "Ensembl")
    gene_ids = data_df["target"].tolist()

    # Record the start time
    version_metadata = get_version()
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
              }
            }
          }
        }
      }
    """
    query_string = query_string.replace("$ids", str(gene_ids).replace("'", '"'))

    r = requests.post(base_url, json={"query": query_string}).json()

    # Record the end time
    end_time = datetime.datetime.now()

    # Update the metadata file
    version_metadata["query"] = {
        "size": len(gene_ids),
        "input_type": "Ensembl",
        "time": str(end_time - start_time),
        "date": str(datetime.datetime.now()),
        "url": base_url,
    }

    data = []

    for gene in r["data"]["targets"]:
        if not gene["knownDrugs"]:
            continue

        drug_info = gene["knownDrugs"]["rows"]
        drug_df = pd.DataFrame(drug_info)

        if drug_df.empty:
            continue

        drug_df[["chembl_id", "drug_name"]] = drug_df["drug"].apply(pd.Series)
        drug_df.drop(columns=["drug"], inplace=True)

        drug_df["target"] = gene["id"]

        drug_df["mechanismOfAction"] = drug_df["mechanismOfAction"].apply(
            lambda x: "inhibits" if "antagonist" in x else "activates"
        )
        drug_df.rename(columns={"mechanismOfAction": "relation"}, inplace=True)
        data.append(drug_df)

    if len(data) == 0:
        return pd.DataFrame(), version_metadata

    opentargets_df = pd.concat(data)

    # Merge the two DataFrames on the target column
    merged_df = collapse_data_sources(
        data_df=data_df,
        source_namespace="Ensembl",
        target_df=opentargets_df,
        common_cols=["target"],
        target_specific_cols=["chembl_id", "drug_name", "relation"],
        col_name="ChEMBL_Drugs",  # TODO: Cross-check if correct name
    )

    return merged_df, version_metadata


def get_targetgene_disease_associations(bridgedb_df: pd.DataFrame) -> Tuple[pd.DataFrame, dict]:
    """Get information about diseases associated with genes based on OpenTargets.

    :param bridgedb_df: BridgeDb output for creating the list of gene ids to query
    :returns: a DataFrame containing the OpenTargets output and dictionary of the query metadata.
    """
    data_df = get_identifier_of_interest(bridgedb_df, "Ensembl")
    gene_ids = data_df["target"].tolist()

    # Record the start time
    version_metadata = get_version()
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

    r = requests.post(base_url, json={"query": query_string}).json()

    # Record the end time
    end_time = datetime.datetime.now()

    # Update the metadata file
    version_metadata["query"] = {
        "size": len(gene_ids),
        "input_type": "Ensembl",
        "time": str(end_time - start_time),
        "date": str(datetime.datetime.now()),
        "url": base_url,
    }

    data = []

    for gene in r["data"]["targets"]:
        if not gene["knownDrugs"]:
            continue
        disease_info = gene["knownDrugs"]["rows"]
        disease_df = pd.DataFrame(disease_info)

        if disease_df.empty:
            continue

        disease_df[["disease_id", "disease_name", "therapeutic_area"]] = disease_df[
            "disease"
        ].apply(pd.Series)
        disease_df["therapeutic_areas"] = disease_df["therapeutic_area"].apply(
            lambda x: ", ".join([f"{i['id']}:{i['name']}" for i in x])
        )
        disease_df.drop(columns=["disease", "therapeutic_area"], inplace=True)

        disease_df["target"] = gene["id"]

        # disease_df['therapeuticAreas'] = disease_df['therapeuticAreas'].apply(
        #   lambda x: ', '.join([i['therapeuticArea'] for i in x])
        # )
        data.append(disease_df)

    if len(data) == 0:
        return pd.DataFrame(), version_metadata

    opentargets_df = pd.concat(data)

    # Merge the two DataFrames on the target column
    merged_df = collapse_data_sources(
        data_df=data_df,
        source_namespace="Ensembl",
        target_df=opentargets_df,
        common_cols=["target"],
        target_specific_cols=["disease_id", "disease_name", "therapeutic_areas"],
        col_name="OpenTargets_Diseases",
    )

    return merged_df, version_metadata


def get_drug_disease_interactions(
    drug_df: pd.DataFrame, disgenet_result: pd.DataFrame
) -> Tuple[pd.DataFrame, dict]:
    """Get information about drugs associated with diseases of interest.

    :param drug_df: get_gene_disease_associations output for creating the list of diseases ids to query
    :param disgenet_result: DisGeNET diseases dadtaframe
    :returns: a DataFrame containing the drug-disease relationships between drug and diseases included in your graph
    """
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
    col_names = ["identifier", "drug_name", "drug_diseases"]
    drug_disease_df = pd.DataFrame(columns=col_names)

    # Extracting and simplifying the structure
    col_names = ["identifier", "drug_name", "drug_diseases"]
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
                "drug_name": entries["name"],
                "drug_diseases": drugs_list,
            }

        if dict_new is not None:
            dict_new_df = pd.DataFrame([dict_new])
            drug_disease_df = pd.concat([drug_disease_df, dict_new_df], ignore_index=True)

    return drug_disease_df
