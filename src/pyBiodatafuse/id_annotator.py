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
    bridgedb_df: pd.DataFrame,
    selected_sources_list: list,
    api_key: Optional[str] = None,
    map_name: Optional[str] = None,
) -> Tuple[pd.DataFrame, DefaultDict[str, dict]]:
    """Query the selected databases and convert the output to a dataframe.

    :param bridgedb_df: BridgeDb output for creating the list of gene ids to query.
    :param selected_sources_list: list of selected databases.
    :param api_key: DisGeNET API key (more details can be found at https://disgenet.com/plans).
    :param map_name: name of the map you want to retrieve the information from. The extensive list \
        can be found at https://minerva-net.lcsb.uni.lu/table.html.
    :returns: a DataFrame containing the combined output and dictionary of the metadata.
    :raises ValueError: If 'disgenet' is in the selected_sources_list and api_key is not provided. \
        Or if 'minerva' is in the selected sources and if map name is not provided.
    """
    # Check if 'disgenet' is in the selected sources and if API key is provided
    if "disgenet" in selected_sources_list and not api_key:
        raise ValueError("API key is required for the 'disgenet' data source.")

    # Check if 'minerva' is in the selected sources and if map name is provided
    if "minerva" in selected_sources_list and not map_name:
        raise ValueError(
            "Map name is required for the 'minerva' data source. See here: https://minerva-net.lcsb.uni.lu/table.html"
        )

    # Initialize variables
    combined_data = pd.DataFrame()
    combined_metadata: DefaultDict[str, dict] = defaultdict(dict)

    # Dictionary to map the datasource names to their corresponding functions
    data_source_functions: Dict[str, Callable[[pd.DataFrame], Tuple[pd.DataFrame, dict]]] = {
        "bgee.gene_expression": bgee.get_gene_expression,
        "disgenet.gene_disease": _get_gene_disease_disgenet_wrapper(api_key),
        "minerva.gene_minerva_pathways": _get_gene_minerva_pathway_wrapper(
            map_name=map_name or "COVID19 Disease Map",
            get_elements=True,
            get_reactions=True,
        ),
        "molmedb.gene_compound": molmedb.get_gene_compound_inhibitor,
        "opentarget.gene_go": opentargets.get_gene_go_process,
        "opentarget.gene_reactome": opentargets.get_gene_reactome_pathways,
        "opentarget.gene_compound": opentargets.get_gene_compound_interactions,
        "opentarget.disease_compound": opentargets.get_disease_compound_interactions,
        "pubchem.protein_compound": pubchem.get_protein_compound_screened,
        "string.protein_protein": stringdb.get_ppi,
        "wikipathways.gene_wikipathways": wikipathways.get_gene_wikipathways,
        # TODO: "wikidata"
    }
    warnings = []

    for source in selected_sources_list:
        if source in data_source_functions:
            tmp_data, tmp_metadata = data_source_functions[source](bridgedb_df)
            combined_metadata[source] = tmp_metadata
            if tmp_data.empty:
                warnings.append(f"No annotation available for {source}")
            if not tmp_data.empty:
                combined_data = combine_sources(bridgedb_df, [combined_data, tmp_data])

    return combined_data, combined_metadata


def _get_gene_disease_disgenet_wrapper(
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


def _get_gene_minerva_pathway_wrapper(
    map_name: str,
    get_elements: Optional[bool] = True,
    get_reactions: Optional[bool] = True,
) -> Callable[[pd.DataFrame], Tuple[pd.DataFrame, dict]]:
    """Create a function to extract gene-minerva pathways data with default parameters.

    :param map_name: name of the map you want to retrieve the information from. The extensive list
        can be found at https://minerva-net.lcsb.uni.lu/table.html.
    :param get_elements: boolean to get elements of the chosen diagram.
    :param get_reactions: if get_reactions = boolean to get reactions of the chosen diagram.
    :returns: A function that takes a DataFrame and returns a tuple containing the annotated DataFrame and metadata dictionary.
    """

    def wrapper(bridgedb_df: pd.DataFrame) -> Tuple[pd.DataFrame, dict]:
        """Extract gene-minerva pathways data.

        :param bridgedb_df: BridgeDb output for creating the list of gene ids to query
        :returns: a tuple containing MINERVA outputs and dictionary of the MINERVA metadata.
        """
        return minerva.get_gene_pathways(
            bridgedb_df,
            map_name=map_name,
            get_elements=get_elements,
            get_reactions=get_reactions,
        )

    return wrapper
