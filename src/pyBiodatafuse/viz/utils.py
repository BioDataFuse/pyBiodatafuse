# coding: utf-8

"""Python utils file for plotting functions."""

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import plotly.express as px

"""Matplotlib plotting functions."""

def plot_pie_chart(
    template_df: pd.DataFrame,
    fig_size: tuple = (10, 10),
) -> plt:
    """Plot a pie chart.

    :param template_df: A dataframe with two columns: "label" and "value"
    :param fig_size: A tuple with the size of the figure
    :returns: A pie chart
    """
    plt.figure(figsize=fig_size)
    colors = sns.color_palette("pastel")[:len(template_df)]
    plt.pie(template_df["value"], labels=template_df["label"], colors=colors)
    return plt.show()


def plot_hbarplot_chart(
    template_df: pd.DataFrame,
    x_label: str = "Label",
    y_label: str = "Value",
    fig_size: tuple = (10, 10),
) -> plt:
    """Plot a bar plot.

    :param template_df: A dataframe with two columns: "label" and "value"
    :param x_label: The x-axis label
    :param y_label: The y-axis label
    :param fig_size: A tuple with the size of the figure
    :returns: A bar plot
    """
    plt.figure(figsize=fig_size)

    ax = sns.barplot(
        data=template_df,
        x="value",
        y="label",
        color="#c7c263",
        orient="h",
    )

    ax.bar_label(ax.containers[0])

    plt.xlabel("Number of patent", size=12)
    plt.ylabel("Patent countries", size=12)
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)

    plt.tight_layout()

    return plt.show()



"""PLotly based functions"""

def plotly_pie_chart(
    template_df: pd.DataFrame,
    x_label: str = "label",
    y_label: str = "value",
    fig_size: tuple = (10, 10),
) -> px:
    """Plot a pie chart using Plotly.

    :param template_df: A dataframe with two columns: "label" and "value"
    :param fig_size: A tuple with the size of the figure
    :returns: A plotly pie chart
    """
    fig = px.pie(template_df, values=y_label, names=x_label)
    return fig.show()


def plotly_barplot_chart(
    template_df: pd.DataFrame,
    x_label: str = "Label",
    y_label: str = "Value",
    fig_size: tuple = (10, 10),
) -> px:
    """Plot a bar plot using Plotly.

    :param template_df: A dataframe with two columns: "label" and "value"
    :param x_label: The x-axis label
    :param y_label: The y-axis label
    :param fig_size: A tuple with the size of the figure
    :returns: A bar plot
    """
    fig = px.bar(template_df, x="value", y="label").update_layout(
        xaxis_title=x_label, yaxis_title=y_label
    )

    return fig.show()