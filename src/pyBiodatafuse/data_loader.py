# coding: utf-8

"""Python script to conver a list of identifiers to a dataframe."""

import re
from typing import Optional

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


def create_df_from_dea(file_path: str) -> pd.DataFrame:
    """Read a dataframe containing the result of the differential expression analysis (DEA).

    :param file_path: path to the file containing the result of DEA
    :returns: the DEA dataframe with proper column name
    :raises ValueError: if the file is not value
    """
    # Get the file extension
    file_extension = file_path.split(".")[-1].lower()
    if file_extension == "xlsx":
        # Read Excel file (xlsx)
        try:
            df = pd.read_excel(file_path)
            df = df.rename(columns={df.columns[0]: "identifier"})
            return df
        except Exception as e:
            raise ValueError(f"Error reading Excel file: {str(e)}")
    if file_extension == "xls":
        # Read Excel file (xls)
        try:
            df = pd.read_excel(file_path, engine="xlrd")
            df = df.rename(columns={df.columns[0]: "identifier"})
            return df
        except Exception as e:
            raise ValueError(f"Error reading Excel file: {str(e)}")
    elif file_extension == "csv" or file_extension == "txt":
        # Read CSV or text file
        try:
            delimiter = "," if file_extension == "csv" else "\t"
            df = pd.read_csv(file_path, sep=delimiter)
            df = df.rename(columns={df.columns[0]: "identifier"})
            return df
        except Exception as e:
            raise ValueError(f"Error reading CSV/text file: {str(e)}")
    else:
        raise ValueError("Unsupported file format. Please provide an Excel, CSV, or TXT file.")


def filter_dea(
    data: pd.DataFrame,
    column_name: str,
    min_value: Optional[float] = None,
    max_value: Optional[float] = None,
    abs_value: Optional[float] = None,
) -> pd.DataFrame:
    """Filter the differential expression analysis (DEA) table.

    :param data: DEA dataframe
    :param column_name: the column to filter
    :param min_value: the minimum value
    :param max_value: the maximum value
    :param abs_value: the absolute value (when filtering for LogFoldChange)
    :returns: the filtered DEA dataframe
    :raises ValueError: if the paramaters are invalid
    """
    if (min_value is not None or max_value is not None) and abs_value is not None:
        raise ValueError(
            "When providing abs_value, min_value and max_value should not be specified"
        )
    elif min_value is None and max_value is None:
        filtered_abs_df = data[abs(data[column_name]) >= abs_value]
        return filtered_abs_df
    if abs_value is None:
        filtered_max_df = pd.DataFrame()
        filtered_abs_df = pd.DataFrame()
        if min_value is not None:
            filtered_min_df = data[data[column_name] >= min_value]
        if max_value is not None:
            filtered_max_df = data[data[column_name] <= max_value]

        filtered_df = pd.concat([filtered_min_df, filtered_max_df])

        return filtered_df
