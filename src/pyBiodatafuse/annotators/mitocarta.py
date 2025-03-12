# coding: utf-8

"""Python file for MitoCarta ETL process:
You can download the MitoCarta dataset from **MitoCarta**. Please visit the following page for the download:

[MitoCarta Download Page](https://personal.broadinstitute.org/scalvo/MitoCarta3.0/)

The datasets you need can be downloaded from the following links:

**For Humans (Homo sapiens):**
- [Human.MitoCarta3.0.xls](https://personal.broadinstitute.org/scalvo/MitoCarta3.0/Human.MitoCarta3.0.xls)

**For Mice (Mus musculus):**
- [Mouse.MitoCarta3.0.xls](https://personal.broadinstitute.org/scalvo/MitoCarta3.0/Mouse.MitoCarta3.0.xls)

These files contain the MitoCarta data in a simple format for each species.
"""

import gzip
import os
from datetime import datetime
from typing import Tuple

import pandas as pd
import requests

from pyBiodatafuse.constants import MITOCARTA, MITOCARTA_DOWNLOAD_URL, MITOCARTA_GENE_INPUT_ID
from pyBiodatafuse.utils import collapse_data_sources, get_identifier_of_interest


def download_mitocarta_dataset(mitocarta_file: str, filename: str, sheet_name:str = "A Human MitoCarta3.0") -> Tuple[pd.DataFrame, dict]:
    """Downloads, saves and reads a MitoCarta dataset.

    :param mitocarta_file: The MitoCarta dataset to download. Human "Human.MitoCarta3.0.xls".
    :param sheet_name: The name of the sheet in the Excel file to read. Default is "A Human MitoCarta3.0".
    :param filename: The local file path to save the downloaded dataset.
    :returns: A MitoCarta DataFrame and dictionary of the MitoCarta metadata.
    """
    # Dowonload the TF-Target dataset
    url = f"{MITOCARTA_DOWNLOAD_URL}/{mitocarta_file}"
    if not os.path.exists(filename):
        response = requests.get(url)
        try:
            response.raise_for_status()
        except requests.HTTPError as e:
            print(f"Failed to download file. HTTP Error: {e}")
        else:
            with open(filename, "wb") as file:
                file.write(response.content)
    
    with gzip.open(filename, "rt") as f:
        mitocarta_df = pd.read_excel(url, sheet_name=sheet_name)

    if mitocarta_df is not None:
        # Add version
        mitocarta_metadata = {
            "datasource": MITOCARTA,
            "metadata": {
                "download date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                "download link": f"{MITOCARTA_DOWNLOAD_URL}/{mitocarta_file}",
            },
        }

        return mitocarta_df, mitocarta_metadata


def process_mitocarta(mitocarta_df: pd.DataFrame, species: str = "hsapiens") -> pd.DataFrame:
    """Add targets and TFs to each row (gene).

    :param ncbi_df: BridgeDb output with ncbi id as target source.
    :param mitocarta_df: The mitocarta dataset.
    :returns: ncbi_df with targets and TFs in each row.
    """
    # Select relevant columns for inclusion in the graph
    if species == "hsapiens":
        selected_columns = [
            "EnsemblGeneID_mapping_version_20200130",
            "Description",
            "MitoCarta3.0_Evidence",
            "MitoCarta3.0_SubMitoLocalization",
            "MitoCarta3.0_MitoPathways",
            "HPA_Main_Location_2020 (Reliability)",
            "Tissues",
        ]
        mitocarta_subset = mitocarta_df[selected_columns]
        # Rename columns for clarity
        mitocarta_subset.rename(
            columns={
                "EnsemblGeneID_mapping_version_20200130": "ensembl_id",
                "Description": "gene_description",
                "MitoCarta3.0_Evidence": "evidence",
                "MitoCarta3.0_SubMitoLocalization": "sub_mito_localization",
                "MitoCarta3.0_MitoPathways": "mito_pathways",
                "HPA_Main_Location_2020 (Reliability)": "hpa_location",
                "Tissues": "tissue_expression",
            },
            inplace=True,
        )
    elif species == "mmusculus":
        selected_columns = [
            "EnsemblGeneID",
            "Description",
            "MitoCarta3.0_Evidence",
            "MitoCarta3.0_SubMitoLocalization",
            "MitoCarta3.0_MitoPathways",
            "HPA_Main_Location_2020 (Reliability)",
            "Tissues",
        ]
        mitocarta_subset = mitocarta_df[selected_columns]
        mitocarta_subset.rename(
            columns={
                "EnsemblGeneID": "ensembl_id",
                "Description": "gene_description",
                "MitoCarta3.0_Evidence": "evidence",
                "MitoCarta3.0_SubMitoLocalization": "sub_mito_localization",
                "MitoCarta3.0_MitoPathways": "mito_pathways",
                "HPA_Main_Location_2020 (Reliability)": "hpa_location",
                "Tissues": "tissue_expression",
            },
            inplace=True,
        )
    # 
    mitocarta_subset["mito_pathways"] = (
        mitocarta_subset["mito_pathways"]
        .str.split(">")
        .str[-1]
        .str.split("|")
        .str[0]
        .str.strip()
    )
    mitocarta_melted = mitocarta_subset.apply(
        lambda row: [row["ensembl_id"], [row.drop("ensembl_id").to_dict()]], axis=1
        )
    intermediate_df = pd.DataFrame(mitocarta_melted.tolist(), columns=["target", MITOCARTA])

    return intermediate_df


def get_gene_mito_pathways(
    bridgedb_df: pd.DataFrame,
    mitocarta_file: str,
    filename: str,
    species: str = "hsapiens",
    sheet_name: str = "A Human MitoCarta3.0"
) -> Tuple[pd.DataFrame, dict]:
    
    """get gene and mitochondia pathways from MitoCarta.

    :param bridgedb_df: BridgeDb output for creating the list of gene ids to query.
    :param mitocarta_file: Name of the remote MitoCarta file to download.
    :param filename: The local file path to save the downloaded dataset.
    :param species: Species for which to process the data; defaults to "hsapiens".
    :param sheet_name: Excel sheet name to read from the file; defaults to "A Human MitoCarta3.0".
    :returns: A tuple containing the processed DataFrame and a metadata dictionary.
    """
    # Download dataset and get metadata
    mitocarta_df, mitocarta_metadata = download_mitocarta_dataset(
        mitocarta_file=mitocarta_file,
        filename=filename,
        sheet_name=sheet_name
    )
    # Process the dataset according to species
    intermediate_df = process_mitocarta(
        mitocarta_df=mitocarta_df,
        species=species
    )
    
    # Merge the processed DataFrame with the original bridgedb_df
    data_df = bridgedb_df[bridgedb_df["target.source"] == MITOCARTA_GENE_INPUT_ID]
    # merged_df = collapse_data_sources(
    #     data_df=bridgedb_df,
    #     source_namespace=MITOCARTA_GENE_INPUT_ID,
    #     target_df=intermediate_df,
    #     common_cols=["target"],
    #     target_specific_cols=list([]),
    #     col_name=MITOCARTA,
    # )
    merged_df = pd.merge(  # TODO: check why the output is not correct  when using collapse_data_sources function
        data_df,
        intermediate_df,
        on="target",
        how="left",
    )

    return merged_df, mitocarta_metadata
