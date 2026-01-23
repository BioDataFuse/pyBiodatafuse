# coding: utf-8

"""Python file for querying the OpenTargets database (https://www.opentargets.org/)."""

import datetime
import warnings
from typing import Dict, Literal, Tuple

import numpy as np
import pandas as pd
import requests
from tqdm import tqdm

import pyBiodatafuse.constants as Cons
from pyBiodatafuse import id_mapper
from pyBiodatafuse.utils import (
    check_columns_against_constants,
    collapse_data_sources,
    get_identifier_of_interest,
    give_annotator_warning,
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
    r = requests.post(Cons.OPENTARGETS_ENDPOINT, json={"query": query}).json()

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
    r = requests.post(Cons.OPENTARGETS_ENDPOINT, json={"query": query}).json()

    year = r["data"]["meta"]["dataVersion"]["year"]
    month = r["data"]["meta"]["dataVersion"]["month"]
    api_version_x = r["data"]["meta"]["apiVersion"]["x"]
    api_version_y = r["data"]["meta"]["apiVersion"]["y"]
    api_version_z = r["data"]["meta"]["apiVersion"]["z"]

    metadata = {
        Cons.DATASOURCE: r["data"]["meta"]["name"],
        Cons.METADATA: {
            "source_version": {
                "apiVersion": f"{api_version_x}.{api_version_y}.{api_version_z}",
            },
            "data_version": f"{year}-{month}",
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
            f"{Cons.OPENTARGETS} GraphQL endpoint is not available. Unable to retrieve data.",
            stacklevel=2,
        )
        return pd.DataFrame(), {}

    data_df = get_identifier_of_interest(bridgedb_df, Cons.OPENTARGETS_GENE_INPUT_ID)
    gene_ids = data_df[Cons.TARGET_COL].tolist()

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

    r = requests.post(Cons.OPENTARGETS_ENDPOINT, json={"query": query_string}).json()

    # Record the end time
    end_time = datetime.datetime.now()

    """Metadata details"""
    # Get the current date and time
    current_date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    # Calculate the time elapsed
    time_elapsed = str(end_time - start_time)

    # Add version, datasource, query, query time, and the date to metadata
    opentargets_version["query"] = {
        "size": len(gene_ids),
        "input_type": Cons.OPENTARGETS_GENE_INPUT_ID,
        "time": time_elapsed,
        "date": current_date,
        "url": Cons.OPENTARGETS_ENDPOINT,
    }

    # Generate the OpenTargets DataFrame
    intermediate_df = pd.DataFrame()

    for gene in tqdm(r["data"]["targets"], desc="Processing gene annotation"):
        terms = [i["term"] for i in gene["geneOntology"]]
        types = [i["aspect"] for i in gene["geneOntology"]]
        path_df = pd.DataFrame(terms)
        path_df[Cons.OPENTARGETS_GO_TYPE] = types
        path_df = path_df.drop_duplicates()
        path_df[Cons.TARGET_COL] = gene["id"]
        intermediate_df = pd.concat([intermediate_df, path_df], ignore_index=True)

    if intermediate_df.empty:
        warnings.warn(
            f"There is no annotation for your input list in {Cons.OPENTARGETS_GO_COL}.",
            stacklevel=2,
        )
        return pd.DataFrame(), opentargets_version

    intermediate_df.rename(
        columns={
            "id": Cons.OPENTARGETS_GO_ID,
            "name": Cons.OPENTARGETS_GO_NAME,
        },
        inplace=True,
    )

    # Check if all keys in df match the keys in OUTPUT_DICT
    check_columns_against_constants(
        data_df=intermediate_df,
        output_dict=Cons.OPENTARGETS_GO_OUTPUT_DICT,
        check_values_in=[Cons.GO],
    )

    # Merge the two DataFrames on the target column
    merged_df = collapse_data_sources(
        data_df=data_df,
        source_namespace=Cons.OPENTARGETS_GENE_INPUT_ID,
        target_df=intermediate_df,
        common_cols=[Cons.TARGET_COL],
        target_specific_cols=list(Cons.OPENTARGETS_GO_OUTPUT_DICT.keys()),
        col_name=Cons.OPENTARGETS_GO_COL,
    )

    """Update metadata"""
    # Calculate the number of new nodes
    num_new_nodes = intermediate_df[Cons.OPENTARGETS_GO_ID].nunique()
    # Calculate the number of new edges
    num_new_edges = intermediate_df.drop_duplicates(
        subset=[Cons.TARGET_COL, Cons.OPENTARGETS_GO_ID]
    ).shape[0]

    # Check the intermediate_df
    if num_new_edges != len(intermediate_df):
        give_annotator_warning(Cons.OPENTARGETS_GO_COL)

    # Add the number of new nodes and edges to metadata
    opentargets_version[Cons.QUERY][Cons.NUM_NODES] = num_new_nodes
    opentargets_version[Cons.QUERY][Cons.NUM_EDGES] = num_new_edges

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
            f"{Cons.OPENTARGETS} GraphQL endpoint is not available. Unable to retrieve data.",
            stacklevel=2,
        )
        return pd.DataFrame(), {}

    data_df = get_identifier_of_interest(bridgedb_df, Cons.OPENTARGETS_GENE_INPUT_ID)
    gene_ids = data_df[Cons.TARGET_COL].tolist()

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

    r = requests.post(Cons.OPENTARGETS_ENDPOINT, json={"query": query_string}).json()

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
        "input_type": Cons.OPENTARGETS_GENE_INPUT_ID,
        "time": time_elapsed,
        "date": current_date,
        "url": Cons.OPENTARGETS_ENDPOINT,
    }

    # Generate the OpenTargets DataFrame
    intermediate_df = pd.DataFrame()

    for gene in tqdm(r["data"]["targets"], desc="Processing gene-pathway interactions"):
        path_df = pd.DataFrame(gene[Cons.PATHWAYS])
        path_df = path_df.drop_duplicates()
        path_df[Cons.TARGET_COL] = gene["id"]
        intermediate_df = pd.concat([intermediate_df, path_df], ignore_index=True)

    if intermediate_df.empty:
        warnings.warn(
            f"There is no annotation for your input list in {Cons.OPENTARGETS_REACTOME_COL}.",
            stacklevel=2,
        )
        return pd.DataFrame(), opentargets_version

    intermediate_df.rename(
        columns={
            "pathway": Cons.PATHWAY_LABEL,
            "pathwayId": Cons.PATHWAY_ID,
        },
        inplace=True,
    )

    # Fixing the pathway_id
    new_ids = []
    for idx in intermediate_df[Cons.PATHWAY_ID]:
        if idx.startswith("R-"):
            new_ids.append(f"{Cons.REACTOME}:{idx}")
        elif idx.startswith("WP"):
            new_ids.append(f"{Cons.WP}:{idx}")
        else:
            print(idx)  # TODO: if this occures, then we need to check the data
            new_ids.append(idx)
    intermediate_df[Cons.PATHWAY_ID] = new_ids

    # Check if all keys in df match the keys in OUTPUT_DICT
    check_columns_against_constants(
        data_df=intermediate_df,
        output_dict=Cons.OPENTARGETS_REACTOME_OUTPUT_DICT,
        check_values_in=[Cons.OPENTARGETS_POSSIBLE_PATHWAY_IDS],
    )

    # Merge the two DataFrames on the target column
    merged_df = collapse_data_sources(
        data_df=data_df,
        source_namespace=Cons.OPENTARGETS_GENE_INPUT_ID,
        target_df=intermediate_df,
        common_cols=[Cons.TARGET_COL],
        target_specific_cols=[Cons.PATHWAY_LABEL, Cons.PATHWAY_ID],
        col_name=Cons.OPENTARGETS_REACTOME_COL,
    )

    """Update metadata"""
    # Calculate the number of new nodes
    num_new_nodes = intermediate_df[Cons.PATHWAY_ID].nunique()
    # Calculate the number of new edges
    num_new_edges = intermediate_df.drop_duplicates(
        subset=[Cons.TARGET_COL, Cons.PATHWAY_ID]
    ).shape[0]

    # Check the intermediate_df
    if num_new_edges != len(intermediate_df):
        give_annotator_warning(Cons.OPENTARGETS_REACTOME_COL)

    # Add the number of new nodes and edges to metadata
    opentargets_version[Cons.QUERY][Cons.NUM_NODES] = num_new_nodes
    opentargets_version[Cons.QUERY][Cons.NUM_EDGES] = num_new_edges

    return merged_df, opentargets_version


def _process_compounds(
    drug_df: pd.DataFrame, target_id: str, target_type: Literal["gene", "disease"]
) -> pd.DataFrame:
    """Process the compounds data from OpenTargets."""
    if target_type == "gene":
        cols = Cons.OPENTARET_COMPOUND_COLS + ["mechanisms_of_action"]
    else:
        cols = Cons.OPENTARET_COMPOUND_COLS
    drug_df[cols] = drug_df["drug"].apply(pd.Series)
    drug_df[Cons.TARGET_COL] = target_id
    drug_df[Cons.CHEMBL_ID] = drug_df[Cons.CHEMBL_ID].astype(str)

    if target_type == "gene":
        drug_df[Cons.OPENTARGETS_COMPOUND_RELATION] = drug_df["mechanismOfAction"].apply(
            lambda x: "inhibits" if "antagonist" in x else "activates"
        )

    drug_df[Cons.DRUGBANK_ID] = drug_df["cross_references"].apply(
        lambda x: (
            next((ref["ids"][0] for ref in x if ref["source"] == "drugbank"), None) if x else None
        )
    )
    drug_df[Cons.DRUGBANK_ID] = drug_df[Cons.DRUGBANK_ID].apply(
        lambda x: f"{Cons.DRUGBANK}:{x}" if x else None
    )

    drug_df[[Cons.OPENTARGETS_ADVERSE_EFFECT_COUNT, Cons.OPENTARGETS_ADVERSE_EFFECT]] = (
        drug_df.apply(
            lambda row: (
                pd.Series([row["adverse_events"]["count"], row["adverse_events"]["rows"]])
                if row["adverse_events"]
                else pd.Series([0, None])
            ),
            axis=1,
        )
    )
    drug_df[Cons.OPENTARGETS_ADVERSE_EFFECT_COUNT] = drug_df[
        Cons.OPENTARGETS_ADVERSE_EFFECT_COUNT
    ].astype(int)

    drug_df.drop(columns=["drug", "cross_references", "adverse_events"], inplace=True)
    return drug_df


def get_gene_compound_interactions(
    bridgedb_df: pd.DataFrame,
    cache_pubchem_cid: bool = True,
) -> Tuple[pd.DataFrame, dict]:
    """Get information about drugs associated with a genes of interest.

    :param bridgedb_df: BridgeDb output for creating the list of gene ids to query
    :param cache_pubchem_cid: whether to cache the PubChem CID for the ChEMBL ID
    :returns: a DataFrame containing the OpenTargets output and dictionary of the query metadata.
    """
    # Check if the API is available
    api_available = check_endpoint_opentargets()
    if not api_available:
        warnings.warn(
            f"{Cons.OPENTARGETS} GraphQL endpoint is not available. Unable to retrieve data.",
            stacklevel=2,
        )
        return pd.DataFrame(), {}

    data_df = get_identifier_of_interest(bridgedb_df, Cons.OPENTARGETS_GENE_INPUT_ID)
    gene_ids = data_df[Cons.TARGET_COL].tolist()

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
                  maximumClinicalTrialPhase
                  crossReferences {
                    source
                    ids
                  }
                  adverseEvents {
                    count
                    rows {
                      name
                    }
                  }
                  mechanismsOfAction {
                    uniqueActionTypes
                    uniqueTargetTypes
                  }
                }
              }
            }
          }
        }

    """
    query_string = query_string.replace("$ids", str(gene_ids).replace("'", '"'))
    r = requests.post(Cons.OPENTARGETS_ENDPOINT, json={"query": query_string}).json()
    # Record the end time
    end_time = datetime.datetime.now()

    """Metadata details"""
    # Get the current date and time
    current_date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    # Calculate the time elapsed
    time_elapsed = str(end_time - start_time)

    # Add version, datasource, query, query time, and the date to metadata
    opentargets_version["query"] = {
        "size": len(gene_ids),
        "input_type": Cons.OPENTARGETS_GENE_INPUT_ID,
        "time": time_elapsed,
        "date": current_date,
        "url": Cons.OPENTARGETS_ENDPOINT,
    }

    # Generate the OpenTargets DataFrame
    intermediate_df = pd.DataFrame()
    if r.get("data") and r.get("data").get("targets"):
        for gene in tqdm(r["data"]["targets"], desc="Processing gene-drug interactions"):
            if not gene["knownDrugs"]:
                continue

            drug_info = gene["knownDrugs"]["rows"]
            drug_df = pd.DataFrame(drug_info)

            if drug_df.empty:
                continue

            drug_df = _process_compounds(drug_df, gene["id"], "gene")

            intermediate_df = pd.concat([intermediate_df, drug_df], ignore_index=True)
            intermediate_df = intermediate_df.drop_duplicates(
                subset=[
                    col
                    for col in intermediate_df.columns
                    if col not in [Cons.OPENTARGETS_ADVERSE_EFFECT, "mechanisms_of_action"]
                ]
            )

    if intermediate_df.empty:
        warnings.warn(
            f"There is no annotation for your input list in {Cons.OPENTARGETS_GENE_COMPOUND_COL}.",
            stacklevel=2,
        )
        return pd.DataFrame(), opentargets_version

    # Fixing chembl_id to pubchem_id
    chembl_ids = intermediate_df[Cons.CHEMBL_ID].values.tolist()
    mapped_df, _ = id_mapper.pubchem_xref(
        identifiers=chembl_ids,
        identifier_type="name",
        cache_res=cache_pubchem_cid,
    )
    mapped_df = mapped_df[[Cons.IDENTIFIER_COL, Cons.TARGET_COL]]
    mapped_dict = mapped_df.set_index(Cons.IDENTIFIER_COL).to_dict()[Cons.TARGET_COL]
    intermediate_df["compound_cid"] = intermediate_df[Cons.CHEMBL_ID].map(mapped_dict)
    intermediate_df[Cons.CHEMBL_ID] = f"{Cons.CHEMBL}:" + intermediate_df[Cons.CHEMBL_ID]

    # Check if all keys in df match the keys in OUTPUT_DICT
    check_columns_against_constants(
        data_df=intermediate_df,
        output_dict=Cons.OPENTARGETS_COMPOUND_OUTPUT_DICT,
        check_values_in=Cons.OPENTARGETS_COMPOUND_VALUE_CHECK_LIST,
    )

    # Merge the two DataFrames on the target column
    merged_df = collapse_data_sources(
        data_df=data_df,
        source_namespace=Cons.OPENTARGETS_GENE_INPUT_ID,
        target_df=intermediate_df,
        common_cols=[Cons.TARGET_COL],
        target_specific_cols=list(Cons.OPENTARGETS_COMPOUND_OUTPUT_DICT.keys()),
        col_name=Cons.OPENTARGETS_GENE_COMPOUND_COL,
    )

    # Ensure the column is correctly assigned
    if Cons.OPENTARGETS_GENE_COMPOUND_COL not in merged_df.columns:
        merged_df[Cons.OPENTARGETS_GENE_COMPOUND_COL] = pd.Series(
            [[] for _ in range(len(merged_df))]
        )

    """Update metadata"""
    # Calculate the number of new nodes
    num_new_nodes = intermediate_df[Cons.CHEMBL_ID].nunique()
    # Calculate the number of new edges
    num_new_edges = intermediate_df.drop_duplicates(subset=[Cons.TARGET_COL, Cons.CHEMBL_ID]).shape[
        0
    ]

    # Check the intermediate_df
    if num_new_edges != len(intermediate_df):
        give_annotator_warning(Cons.OPENTARGETS_GENE_COMPOUND_COL)

    # Add the number of new nodes and edges to metadata
    opentargets_version[Cons.QUERY][Cons.NUM_NODES] = num_new_nodes
    opentargets_version[Cons.QUERY][Cons.NUM_EDGES] = num_new_edges
    return merged_df, opentargets_version


def get_disease_compound_interactions(
    bridgedb_df: pd.DataFrame,
    cache_pubchem_cid: bool = False,
) -> Tuple[pd.DataFrame, dict]:
    """Get information about drugs associated with diseases of interest.

    :param bridgedb_df: BridgeDb output for creating the list of gene ids to query.
    :param cache_pubchem_cid: If True, the PubChem CID will be cached for future use.
    :returns: a DataFrame containing the OpenTargets output and dictionary of the query metadata.
    """
    # Check if the API is available
    api_available = check_endpoint_opentargets()
    if not api_available:
        warnings.warn(
            f"{Cons.OPENTARGETS} GraphQL endpoint is not available. Unable to retrieve data.",
            stacklevel=2,
        )
        return pd.DataFrame(), {}

    if bridgedb_df.empty:
        warnings.warn(
            "There is no input.",
            stacklevel=2,
        )
        return pd.DataFrame(), {}

    data_df = bridgedb_df[bridgedb_df[Cons.TARGET_SOURCE_COL] == Cons.OPENTARGETS_DISEASE_INPUT_ID]
    efo_ids = data_df[Cons.TARGET_COL].tolist()

    # Record the start time
    opentargets_version = get_version_opentargets()
    start_time = datetime.datetime.now()

    query_string = """
    query DiseaseDrugs{
        diseases (efoIds: $efoIds) {
            id
            name
            knownDrugs {
                rows {
                    drug {
                        id
                        name
                        isApproved
                        maximumClinicalTrialPhase
                        crossReferences {
                            ids
                            source
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
    }"""
    final_data = []

    # query in batches of 25
    for i in range(0, len(efo_ids), 25):
        batch_ids = efo_ids[i : i + 25]
        query_string = query_string.replace("$efoIds", str(batch_ids).replace("'", '"'))
        r = requests.post(Cons.OPENTARGETS_ENDPOINT, json={"query": query_string}).json()
        if not r.get("data", {}).get("diseases") or r is None:
            continue

        final_data.append(r)

    # Record the end time
    end_time = datetime.datetime.now()

    """Metdata details"""
    # Get the current date and time
    current_date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    # Calculate the time elapsed
    time_elapsed = str(end_time - start_time)

    # Add version, datasource, query, query time, and the date to metadata
    opentargets_version["query"] = {
        "size": len(efo_ids),
        "input_type": Cons.OPENTARGETS_DISEASE_INPUT_ID,
        "time": time_elapsed,
        "date": current_date,
        "url": Cons.OPENTARGETS_ENDPOINT,
    }

    # Generate the OpenTargets DataFrame
    intermediate_df = pd.DataFrame()

    if len(final_data) == 0:
        warnings.warn(
            f"There is no annotation for your input list in {Cons.OPENTARGETS_DISEASE_COMPOUND_COL}.",
            stacklevel=2,
        )
        return pd.DataFrame(), opentargets_version

    for data in tqdm(final_data, desc="Processing diseases-drug interactions"):
        for disease in data["data"]["diseases"]:
            if not disease["knownDrugs"]:
                continue

            # Based on clinical trial data
            drug_info = disease["knownDrugs"]["rows"]
            drug_df = pd.DataFrame(drug_info)

            if drug_df.empty:
                continue

            drug_df = _process_compounds(drug_df, disease["id"], "disease")

            intermediate_df = pd.concat([intermediate_df, drug_df], ignore_index=True)
            intermediate_df = intermediate_df.drop_duplicates(
                subset=[
                    col
                    for col in intermediate_df.columns
                    if col not in [Cons.OPENTARGETS_ADVERSE_EFFECT, "mechanisms_of_action"]
                ]
            )

    if intermediate_df.empty:
        warnings.warn(
            f"There is no annotation for your input list in {Cons.OPENTARGETS_DISEASE_COMPOUND_COL}.",
            stacklevel=2,
        )
        return pd.DataFrame(), opentargets_version

    # Fixing chembl_id to pubchem_id
    chembl_ids = intermediate_df[Cons.CHEMBL_ID].values.tolist()
    mapped_df, _ = id_mapper.pubchem_xref(
        identifiers=chembl_ids,
        identifier_type="name",
        cache_res=cache_pubchem_cid,
    )

    mapped_df = mapped_df[[Cons.IDENTIFIER_COL, Cons.TARGET_COL]]
    mapped_dict = mapped_df.set_index(Cons.IDENTIFIER_COL).to_dict()[Cons.TARGET_COL]
    intermediate_df["compound_cid"] = intermediate_df[Cons.CHEMBL_ID].map(mapped_dict)
    intermediate_df[Cons.CHEMBL_ID] = f"{Cons.CHEMBL}:" + intermediate_df[Cons.CHEMBL_ID]
    intermediate_df[Cons.OPENTARGETS_COMPOUND_RELATION] = Cons.OPENTARGETS_COMPOUND_DISEASE_RELATION

    # Check if all keys in df match the keys in OUTPUT_DICT
    check_columns_against_constants(
        data_df=intermediate_df,
        output_dict=Cons.OPENTARGETS_COMPOUND_OUTPUT_DICT,
        check_values_in=Cons.OPENTARGETS_COMPOUND_VALUE_CHECK_LIST,
    )

    # Merge the two DataFrames on the target column
    merged_df = collapse_data_sources(
        data_df=data_df,
        source_namespace=Cons.OPENTARGETS_DISEASE_INPUT_ID,
        target_df=intermediate_df,
        common_cols=[Cons.TARGET_COL],
        target_specific_cols=list(Cons.OPENTARGETS_COMPOUND_OUTPUT_DICT.keys()),
        col_name=Cons.OPENTARGETS_DISEASE_COMPOUND_COL,
    )

    # Fill in missing values for columns in bridgedb_df
    if len(merged_df) != len(bridgedb_df):
        subset_df = bridgedb_df[~bridgedb_df[Cons.TARGET_COL].isin(merged_df[Cons.TARGET_COL])]
        fill_data = [[{i: np.nan for i in Cons.OPENTARGETS_COMPOUND_OUTPUT_DICT.keys()}]] * len(
            subset_df
        )
        subset_df[Cons.OPENTARGETS_DISEASE_COMPOUND_COL] = fill_data
        merged_df = pd.concat([merged_df, subset_df], ignore_index=True)

    """Update metadata"""
    # Calculate the number of new nodes
    num_new_nodes = intermediate_df[Cons.CHEMBL_ID].nunique()
    # Calculate the number of new edges
    num_new_edges = intermediate_df.drop_duplicates(subset=[Cons.TARGET_COL, Cons.CHEMBL_ID]).shape[
        0
    ]

    # Check the intermediate_df
    if num_new_edges != len(intermediate_df):
        give_annotator_warning(Cons.OPENTARGETS_DISEASE_COMPOUND_COL)

    # Add the number of new nodes and edges to metadata
    opentargets_version[Cons.QUERY][Cons.NUM_NODES] = num_new_nodes
    opentargets_version[Cons.QUERY][Cons.NUM_EDGES] = num_new_edges
    return merged_df, opentargets_version
