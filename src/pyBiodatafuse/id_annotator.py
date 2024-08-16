from collections import defaultdict

import pandas as pd

from pyBiodatafuse.annotators import disgenet, opentargets, stringdb, wikipathways
from pyBiodatafuse.utils import combine_sources


def process_selected_sources(
    bridgedb_df: pd.DataFrame, selected_sources_list: list
) -> pd.DataFrame:
    """Query the selected databases and convert the output to a dataframe.

    :param bridgedb_df: BridgeDb output for creating the list of gene ids to query.
    :param selected_sources_list: list of selected databases
    :returns: a DataFrame containing the combined output and dictionary of the metadata.
    """
    # Initialize variables
    combined_data = pd.DataFrame()
    combined_metadata = defaultdict(lambda: defaultdict(str))
    # Dictionary to map the datasource names to their corresponding functions
    data_source_functions = {
        # TODO: "bgee": bgee.get_gene_expression,
        "disgenet": disgenet.get_gene_disease,
        # TODO: "minerva": minera.get,
        # TODO: "molmedb": molmedb.get_mol_gene_inhibitor,
        # TODO: "pubchem"
        "opentarget.gene_ontology": opentargets.get_gene_go_process,
        "opentarget.reactome": opentargets.get_gene_reactome_pathways,
        "opentarget.drug_interactions": opentargets.get_gene_compound_interactions,
        "opentarget.disease_associations": opentargets.get_gene_disease_interactions,
        "string": stringdb.get_ppi,
        "wikipathways": wikipathways.get_gene_wikipathway
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
