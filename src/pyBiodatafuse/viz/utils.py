# coding: utf-8

"""Python utils file for plotting functions."""

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def plot_pie_chart(template_df: pd.DataFrame) -> plt:
    """Plot a pie chart.

    :param template_df: A dataframe with two columns: "label" and "value"
    :returns: A pie chart
    """
    plt.figure(figsize=(10, 10))
    colors = sns.color_palette("pastel")[0:5]
    plt.pie(template_df["value"], labels=template_df["label"], colors=colors)
    return plt.show()


def plot_barplot_chart(template_df: pd.DataFrame) -> plt:
    """Plot a pie chart.

    :param template_df: A dataframe with two columns: "label" and "value"
    :returns: A pie chart
    """
    sns.barplot(data=template_df, x="label", y="value")
    return plt.show()
