# coding: utf-8

"""
Python file for MitoCarta ETL process.

You can download the MitoCarta dataset from **MitoCarta**. Please visit the following page for the download:

[MitoCarta Download Page](https://personal.broadinstitute.org/scalvo/MitoCarta3.0/)

The datasets you need can be downloaded from the following links:

**For Humans (Homo sapiens):**
- [Human.MitoCarta3.0.xls](https://personal.broadinstitute.org/scalvo/MitoCarta3.0/Human.MitoCarta3.0.xls)

**For Mice (Mus musculus):**
- [Mouse.MitoCarta3.0.xls](https://personal.broadinstitute.org/scalvo/MitoCarta3.0/Mouse.MitoCarta3.0.xls)

These files contain the MitoCarta data in a simple format for each species.
"""

import os
from datetime import datetime
from typing import Tuple

import pandas as pd
import requests

import pyBiodatafuse.constants as Cons
from pyBiodatafuse.utils import collapse_data_sources, get_identifier_of_interest


def download_mitocarta_dataset(
    mitocarta_file: str, filename: str, sheet_name: str = "A Human MitoCarta3.0"
) -> Tuple[pd.DataFrame, dict]:
    """Download, save, and read a MitoCarta dataset.

    :param mitocarta_file: The MitoCarta dataset to download. Human "Human.MitoCarta3.0.xls".
    :param sheet_name: The name of the sheet in the Excel file to read. Default is "A Human MitoCarta3.0".
    :param filename: The local file path to save the downloaded dataset.
    :returns: A MitoCarta DataFrame and dictionary of the MitoCarta metadata.
    :raises ValueError: If the file cannot be downloaded.
    """
    # Dowonload the TF-Target dataset
    url = f"{Cons.MITOCARTA_DOWNLOAD_URL}/{mitocarta_file}"
    if not os.path.exists(filename):
        response = requests.get(url)
        try:
            response.raise_for_status()
        except requests.HTTPError as e:
            raise ValueError(f"Failed to download file. HTTP Error: {e}")
        else:
            with open(filename, "wb") as file:
                file.write(response.content)

    mitocarta_df = pd.read_excel(filename, sheet_name=sheet_name)

    if mitocarta_df is not None:
        # Add version
        mitocarta_metadata = {
            "datasource": Cons.MITOCARTA,
            "metadata": {
                "download date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                "download link": f"{Cons.MITOCARTA_DOWNLOAD_URL}/{mitocarta_file}",
            },
        }

        return mitocarta_df, mitocarta_metadata

    # Return empty DataFrame and metadata if mitocarta_df is None
    empty_metadata = {
        "datasource": Cons.MITOCARTA,
        "metadata": {
            "download date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "download link": f"{Cons.MITOCARTA_DOWNLOAD_URL}/{mitocarta_file}",
            "note": "No data found in the downloaded file",
        },
    }
    return pd.DataFrame(), empty_metadata


def process_mitocarta(mitocarta_df: pd.DataFrame, species: str = "hsapiens") -> pd.DataFrame:
    """Add targets and TFs to each row (gene).

    :param mitocarta_df: The mitocarta dataset.
    :param species: The species to process the data for; defaults to "hsapiens".
    :returns: mitocarta_df with targets and TFs in each row.
    :raises ValueError: If species is not supported.
    """
    # Select relevant columns for inclusion in the graph
    if species == "hsapiens":
        selected_columns = Cons.MITO_SELECTED_COLUMNS["human"]
    elif species == "mmusculus":
        selected_columns = Cons.MITO_SELECTED_COLUMNS["mouse"]

    else:
        raise ValueError(f"Species {species} not supported.")

    # rename columns
    mitocarta_subset = mitocarta_df[selected_columns]
    mitocarta_subset.rename(columns=Cons.MITOCART_COL_MAPPER, inplace=True, errors="ignore")

    mitocarta_subset[Cons.MITO_PATHWAYS] = (
        mitocarta_subset[Cons.MITO_PATHWAYS]
        .str.split(">")
        .str[-1]
        .str.split("|")
        .str[0]
        .str.strip()
    )

    return mitocarta_subset


def get_gene_mito_pathways(
    bridgedb_df: pd.DataFrame,
    mitocarta_file: str,
    filename: str,
    species: str = "hsapiens",
    sheet_name: str = "A Human MitoCarta3.0",
) -> Tuple[pd.DataFrame, dict]:
    """Get gene and mitochondia pathways from MitoCarta.

    :param bridgedb_df: BridgeDb output for creating the list of gene ids to query.
    :param mitocarta_file: Name of the remote MitoCarta file to download.
    :param filename: The local file path to save the downloaded dataset.
    :param species: Species for which to process the data; defaults to "hsapiens".
    :param sheet_name: Excel sheet name to read from the file; defaults to "A Human MitoCarta3.0".
    :returns: A tuple containing the processed DataFrame and a metadata dictionary.
    """
    # Download dataset and get metadata
    mitocarta_df, mitocarta_metadata = download_mitocarta_dataset(
        mitocarta_file=mitocarta_file, filename=filename, sheet_name=sheet_name
    )

    # Subset the dataset according to species
    subset_df = process_mitocarta(mitocarta_df=mitocarta_df, species=species)

    # Merge the processed DataFrame with the original bridgedb_df
    data_df = get_identifier_of_interest(bridgedb_df, Cons.MITOCARTA_GENE_INPUT_ID)

    merged_df = collapse_data_sources(
        data_df=data_df,
        source_namespace=Cons.MITOCARTA_GENE_INPUT_ID,
        target_df=subset_df,
        common_cols=[Cons.TARGET_COL],
        target_specific_cols=Cons.MITOCART_OUTPUT,
        col_name=Cons.MITOCART_PATHWAY_COL,
    )

    return merged_df, mitocarta_metadata
