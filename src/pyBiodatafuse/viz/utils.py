# coding: utf-8

"""Python utils file for plotting functions."""

import matplotlib.pyplot as plt
import pandas as pd
import plotly as px
import seaborn as sns


def plot_pie_chart():
    pass


if __name__ == "__main__":
    from pyBiodatafuse import id_mapper

    chemicals_of_interest = """124211730
    99997592
    99995871
    99990871"""
    chem_list = chemicals_of_interest.split("\n")

    data_input = pd.DataFrame(chem_list, columns=["identifier"])

    bridgdb_df, bridgdb_metadata = id_mapper.bridgedb_xref(
        identifiers=data_input["identifier"],
        input_species="Human",
        input_datasource="HGNC",
        output_datasource="All",
    )
