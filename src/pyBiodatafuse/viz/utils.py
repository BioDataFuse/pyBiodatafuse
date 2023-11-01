# coding: utf-8

"""Python utils file for plotting functions."""

import matplotlib.pyplot as plt


def plot_pie_chart():
    """Plot a pie chart."""
    labels = "Frogs", "Hogs", "Dogs", "Logs"
    sizes = [15, 30, 45, 10]

    fig, ax = plt.subplots()
    ax.pie(sizes, labels=labels)
    plt.show()
