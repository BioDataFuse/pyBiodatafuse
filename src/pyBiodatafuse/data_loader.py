# coding: utf-8

"""Python script to conver a list of identifiers to a dataframe."""

import re

import pandas as pd


def create_df_from_file(file_path: str) -> pd.DataFrame:
    """Create a DataFrame from a file containing a list of identifiers.

    :param file_path: path to the file containing the list of identifiers
    :returns: a DataFrame containing the list of identifiers
    """
    # Initialize an empty list to store the data
    data = []

    # Open the file and read its contents
    with open(file_path, "r") as file:
        content = file.read()
        # Split the content using regular expressions to handle multiple delimiters (',' and '\n')
        identifiers = [val.strip() for val in re.split(r"[,\n]+", content) if val.strip()]
        data.extend(identifiers)

    # Create a DataFrame using pandas
    df = pd.DataFrame(data, columns=["identifier"])

    return df


def create_df_from_text(text_input: str) -> pd.DataFrame:
    """Create a DataFrame from a text containing a list of identifiers.

    :param text_input: text containing the list of identifiers with each identifier on a new line.
    :returns: a DataFrame containing the list of identifiers
    """
    # Initialize an empty list to store the data
    data = []

    # Split the text using newline characters to create a list of identifiers
    identifiers = [val.strip() for val in text_input.split("\n") if val.strip()]
    data.extend(identifiers)

    # Create a DataFrame using pandas
    df = pd.DataFrame(data, columns=["identifier"])

    return df
