# coding: utf-8

"""Python file for querying the OpenTargets database (https://www.opentargets.org/)."""

import datetime
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


def get_gene_disease_associations(bridgedb_df: pd.DataFrame) -> Tuple[pd.DataFrame, dict]:
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
