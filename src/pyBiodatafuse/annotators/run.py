# coding: utf-8

"""Python file for running the annotators."""

import pandas as pd
from tqdm import tqdm

from pyBiodatafuse import id_mapper
from pyBiodatafuse.annotators import *


def run_gene_annotators(bridgdb_df, bridgdb_metadata, selected_annotators):
    combined_data = []

    for annotator_name in tqdm(selected_annotators):
        if annotator_name == "bgee":
            data_df, metadata = annotator_name.get_gene_expression(
                bridgedb_df=bridgdb_df, anatomical_entities=pd.DataFrame()  # TODO
            )

        elif annotator_name == "disgenet":
            data_df, metadata = annotator_name.get_gene_disease(
                bridgedb_df=bridgdb_df,
            )

        elif annotator_name == "opentargets":
            data_df, metadata = annotator_name.get_gene_disease(
                bridgedb_df=bridgdb_df,
            )

        combined_data.append(data_df)

    return pd.DataFrame(combined_data)


def run_disease_annotators():
    pass


def run_compound_annotators():
    pass


def run_annotator(
    data_input: pd.DataFrame,
    input_specie: str = "Human",
    input_type: str = "HGNC",
    selected_annotators: list = ["disgenet", "opentargets", "bgee"],
):
    assert isinstance(
        data_input, pd.DataFrame
    ), "data_input parameter must be a pandas DataFrame. Please check your input."

    assert input_specie in [
        "Human"
    ], "specie_name parameter must be one of the following: 'Human'. Please check your input."

    # TODO: Add here the check for the input_type parameter (something from BridgeDb)

    bridgdb_df, bridgdb_metadata = id_mapper.bridgedb_xref(
        identifiers=data_input,
        input_species=input_specie,
        input_datasource=input_type,
        output_datasource="All",
    )

    if input_type in ["HGNC"]:
        run_gene_annotators(bridgdb_df, bridgdb_metadata, selected_annotators)

    elif input_type in ["Pubchem"]:
        run_compound_annotators(
            data_input=pd.DataFrame(),
            input_specie="Human",
            input_type="Ensembl",
            annotators=["disgenet", "opentargets", "bgee"],
            annotator_args={},
        )
