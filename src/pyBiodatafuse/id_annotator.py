# coding: utf-8

"""Python file to annotate an input list with selected data sources."""

from collections import defaultdict
from typing import Callable, DefaultDict, Dict, Optional, Tuple

import pandas as pd

from pyBiodatafuse.annotators import (
    bgee,
    disgenet,
    minerva,
    molmedb,
    opentargets,
    pubchem,
    stringdb,
    wikipathways,
)
from pyBiodatafuse.utils import combine_sources


def run_gene_selected_sources(
    bridgedb_df: pd.DataFrame, selected_sources_list: list, api_key: Optional[str] = None
) -> Tuple[pd.DataFrame, DefaultDict[str, dict]]:
    """Query the selected databases and convert the output to a dataframe.

    :param bridgedb_df: BridgeDb output for creating the list of gene ids to query.
    :param selected_sources_list: list of selected databases
    :param api_key: DisGeNET API key (more details can be found at https://disgenet.com/plans)
    :returns: a DataFrame containing the combined output and dictionary of the metadata.
    :raises ValueError: If 'disgenet' is in the selected_sources_list and api_key is not provided.
    """
    # Check if 'disgenet' is in the selected sources and if API key is provided
    if "disgenet" in selected_sources_list and not api_key:
        raise ValueError("API key is required for the 'disgenet' data source.")

    # Initialize variables
    combined_data = pd.DataFrame()
    combined_metadata: DefaultDict[str, dict] = defaultdict(dict)

    # Dictionary to map the datasource names to their corresponding functions
    data_source_functions: Dict[str, Callable[[pd.DataFrame], Tuple[pd.DataFrame, dict]]] = {
        "bgee": bgee.get_gene_expression,
        "disgenet": _get_gene_disease_api_function(api_key),
        "minerva": minerva.get_gene_minerva_pathways,
        "molmedb": molmedb.get_gene_compound_inhibitor,
        "opentarget.gene_ontology": opentargets.get_gene_go_process,
        "opentarget.reactome": opentargets.get_gene_reactome_pathways,
        "opentarget.drug_interactions": opentargets.get_gene_compound_interactions,
        "opentarget.disease_associations": opentargets.get_gene_disease_associations,
        "pubchem": pubchem.get_protein_molecule_screened,
        "string": stringdb.get_ppi,
        "wikipathways": wikipathways.get_gene_wikipathways
        # TODO: "wikidata"
    }
    warnings = []  # Initialize empty list for warnings

    for source in selected_sources_list:
        if source in data_source_functions:
            tmp_data, tmp_metadata = data_source_functions[source](bridgedb_df)
            combined_metadata[source] = tmp_metadata
            if tmp_data.empty:
                warnings.append(f"No annotation available for {source}")
            if not tmp_data.empty:
                combined_data = combine_sources([combined_data, tmp_data])

    return combined_data, combined_metadata


def _get_gene_disease_api_function(
    api_key: Optional[str] = None,
) -> Callable[[pd.DataFrame], Tuple[pd.DataFrame, dict]]:
    """Extract gene-disease data from DisGeNET using the provided API key.

    :param api_key: DisGeNET API key (more details can be found at https://disgenet.com/plans)
    :returns: A function that takes a DataFrame and returns a tuple containing the annotated DataFrame and metadata dictionary.
    """

    def wrapper(bridgedb_df: pd.DataFrame) -> Tuple[pd.DataFrame, dict]:
        """Extract gene-disease data from DisGeNET using the provided API key.

        :param bridgedb_df: BridgeDb output for creating the list of gene ids to query.
        :returns: a DataFrame containing the DisGeNET output and dictionary of the DisGeNET metadata.
        """
        return disgenet.get_gene_disease(bridgedb_df, api_key)

    return wrapper
