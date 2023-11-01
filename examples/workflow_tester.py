# coding: utf-8

"""Python file for testing the workflow package."""

import pandas as pd

from pyBiodatafuse import id_mapper
from pyBiodatafuse.viz.patent_data import _process_data_for_plot, get_patent_data
from pyBiodatafuse.viz.utils import plot_barplot_chart, plot_pie_chart


def main():
    """Run the workflow."""

    chemicals_of_interest = """glucose
    Aspirin"""
    chem_list = chemicals_of_interest.split("\n")

    pubchem_df, pubchem_metadata = id_mapper.pubchem_xref(
        identifiers=chem_list, indentifier_type="name"
    )

    m = get_patent_data(pubchem_df)
    df = _process_data_for_plot(m)

    entries = df.index.to_list()
    test = df.loc[entries[0]]
    # print(test)
    plot_barplot_chart(test)


if __name__ == "__main__":
    main()
