# coding: utf-8

"""Python file for extracting patent data.

The module contains special functions that are server expensive and can only be performed for smaller datasets.
"""

import time
from typing import Literal, Union

import pandas as pd
import requests
from tqdm import tqdm
import plotly.express as px
import matplotlib.pyplot as plt
from pyBiodatafuse.analyzer.utils import (
    plot_hbarplot_chart,
    plot_pie_chart,
    plotly_barplot_chart,
    plotly_pie_chart,
)
from pyBiodatafuse.constants import PATENT_INPUT_ID
from pyBiodatafuse.utils import get_identifier_of_interest


def process_patent_data(patent_dict: dict) -> pd.DataFrame:
    """Process patent data into a DataFrame.

    :param patent_dict: A dictionary with the patent data
    :returns: A DataFrame with the processed patent data
    """
    data_df = []

    for cmpd_id, patent_set in patent_dict.items():
        for patent_country, patents in patent_set.items():
            data_df.append(
                {
                    "compound_id": cmpd_id,
                    "label": patent_country,
                    "value": len(patents),
                }
            )

    data_df = pd.DataFrame(data_df)
    data_df.set_index("compound_id", inplace=True)
    return data_df


def get_patent_from_pubchem(bridgedb_df: pd.DataFrame) -> pd.DataFrame:
    """Get patent data summary from PubChem compounds.

    The output is the following: {CID: ["US: X", "EP: X", "WO: X", "Others: X"]}
    :param bridgedb_df: A dataframe with the BridgeDb or Pubchem harmonized output
    :returns: A dictionary with the PubChem Compound ID as key and the patent counts as value
    """
    # Get column of interest
    data_df = get_identifier_of_interest(bridgedb_df, PATENT_INPUT_ID)

    cid_pat_dict = {}

    # TODO: Neeed some pre-processing:
    # 1. Remove duplicates (WO-03078408-A1 against WO03078408A1, Not classic Patent offices such as AR, AU againts WO)
    # 2. Adding Granted, non-granted patent counts

    c = 0

    for cid in tqdm(data_df["target"], desc="Getting patent data from PubChem"):
        patent_detail_dict = {"US": set(), "EP": set(), "WO": set(), "Others": set()}  # type: ignore[var-annotated]

        try:
            url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid.split(':')[1]}/xrefs/PatentID/TXT"
            patent_list = requests.get(url)
        except requests.exceptions.ConnectionError:
            cid_pat_dict[cid] = patent_detail_dict
            continue

        if patent_list.status_code != 200 or patent_list.text.strip() == "":
            cid_pat_dict[cid] = patent_detail_dict
            continue

        patent_df = pd.DataFrame(patent_list.text.splitlines(), columns=["PatentID"])
        for patent in patent_df["PatentID"].tolist():
            p = patent.replace("-", "")
            if p.startswith("US"):
                patent_detail_dict["US"].add(p)
            elif p.startswith("EP"):
                patent_detail_dict["EP"].add(p)
            elif p.startswith("WO"):
                patent_detail_dict["WO"].add(p)
            else:
                patent_detail_dict["Others"].add(p)

        c += 1
        if c == 100:
            time.sleep(3)
            c = 0
        cid_pat_dict[cid] = patent_detail_dict

    return process_patent_data(cid_pat_dict)


def plot_patent_summary(
    data_df: pd.DataFrame,
    compound_id: str = "",
    fig_size: tuple = (10, 5),
    interactive: bool = True,
    plot_style: Literal["bar", "pie"] = "bar",
) -> Union[px.Figure, plt.Figure, str]:
    """Plot patent summary data.

    :param data_df: A dataframe with two columns: "label" and "value"
    :param compound_id: The compound identifier for the title
    :param fig_size: A tuple with the size of the figure
    :param interactive: Whether to create an interactive plotly plot or a static matplotlib plot
    :param plot_style: The style of the plot, either "bar" or "pie"
    :returns: A plotly or matplotlib plot depending on the interactive parameter
    """
    if compound_id == "":
        return "Please provide a compound identifier for the title."

    if compound_id not in data_df.index:
        return f"{compound_id} not found in the patent data."

    patent_df = data_df.loc[compound_id]

    # Drop entries with zero values
    patent_df = patent_df[patent_df["value"] > 0]

    if patent_df.empty:
        return "No patent data available to plot."

    if plot_style not in ["bar", "pie"]:
        return f"Plot style {plot_style} not recognized. Please use 'bar' or 'pie'."

    if plot_style == "pie" and interactive:
        fig = plotly_pie_chart(
            template_df=patent_df,
            fig_size=fig_size,
            fig_title=f"Patent summary for {compound_id}",
        )
        fig.update_layout(title_text=f"Patent summary for {compound_id}", showlegend=True)

        return fig

    elif plot_style == "bar" and interactive:
        fig = plotly_barplot_chart(
            template_df=patent_df,
            x_label="Number of Patents",
            y_label="Patent Offices",
            fig_size=fig_size,
            fig_title=f"Patent summary for {compound_id}",
        )
        fig.update_layout(title_text=f"Patent summary for {compound_id}", showlegend=False)

        return fig

    elif plot_style == "pie" and not interactive:
        return plot_pie_chart(
            template_df=patent_df,
            fig_size=fig_size,
            fig_title=f"Patent summary for {compound_id}",
        )

    elif plot_style == "bar" and not interactive:
        return plot_hbarplot_chart(
            template_df=patent_df,
            x_label="Number of Patents",
            y_label="Patent Offices",
            fig_size=fig_size,
            fig_title=f"Patent summary for {compound_id}",
        )

    return "No patent data available to plot."
