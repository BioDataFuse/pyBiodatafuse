# coding: utf-8

"""
Python file for TFLink ETL process.

You can download the Gene-TF interactions dataset from **TFLink**. Please visit the following page for the download: [TFLink Download Page](https://tflink.net/download/)

The datasets you need can be downloaded from the following links:

**For Humans (Homo sapiens):**
- [TFLink_Homo_sapiens_interactions_All_simpleFormat_v1.0.tsv.gz](https://cdn.netbiol.org/tflink/download_files/TFLink_Homo_sapiens_interactions_All_simpleFormat_v1.0.tsv.gz)

**For Mice (Mus musculus):**
- [TFLink_Mus_musculus_interactions_All_simpleFormat_v1.0.tsv.gz](https://cdn.netbiol.org/tflink/download_files/TFLink_Mus_musculus_interactions_All_simpleFormat_v1.0.tsv.gz)

These files contain the TF-target interaction data in a simple format for each species.
"""

import gzip
import os
from datetime import datetime
from typing import Optional, Tuple

import pandas as pd
import requests

import pyBiodatafuse.constants as Cons
from pyBiodatafuse.utils import get_identifier_of_interest


def load_tflink_dataset(tf_file: str) -> Tuple[pd.DataFrame, dict]:
    """Load a TFLink dataset.

    :param tf_file: Path of the TF-Target dataset. See the links in the docstring for downloading.
    :returns: A TFLink DataFrame and dictionary of the TFLink metadata.
    """
    with gzip.open(tf_file, "rt") as f:
        tflink_df = pd.read_csv(f, sep="\t")

    # Add version
    tflink_metadata = {
        "datasource": Cons.TFLINK,
        "metadata": {
            "download date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "download link": f"{Cons.TFLINK_DOWNLOAD_URL}/{tf_file}",
        },
    }

    return tflink_df, tflink_metadata


def add_target_and_tf_interaction(
    ncbi_df: pd.DataFrame,
    tflink_df: pd.DataFrame,
    padj_filter: Optional[float],
    padj_colname: Optional[str],
) -> pd.DataFrame:
    """Add targets and TFs to each row (gene).

    :param ncbi_df: BridgeDb output with ncbi id as target source.
    :param tflink_df: The TF-Target dataset.
    :param padj_filter: The adjusted p-value threshold for filtering DEGs.
    :param padj_colname: The name of the column containing adjusted p-values.
    :returns: ncbi_df with targets and TFs in each row.
    """
    ncbi_df[Cons.ITS_TARGET_COL] = None
    ncbi_df[Cons.ITS_TF_COL] = None

    for index, row in ncbi_df.iterrows():
        if not pd.isna(row[Cons.IS_TF_COL]):
            targets = tflink_df[tflink_df[Cons.TFLINK_GENE_ID_TF] == row[Cons.TARGET_COL]]
            if targets.empty:
                ncbi_df.at[index, Cons.ITS_TARGET_COL] = []
            else:
                target_info_list = targets[Cons.TF_COLS_TO_KEEP].to_dict(orient="records")
                ncbi_df.at[index, Cons.ITS_TARGET_COL] = target_info_list

        if row[Cons.IS_TARGET_COL]:
            if padj_colname is None:
                continue

            if row[f"{padj_colname}_dea"] <= padj_filter:
                tf = tflink_df[tflink_df[Cons.TFLINK_GENE_ID_TARGET] == row[Cons.TARGET_COL]]
                tf_info_list = tf[Cons.TARGET_COLS_TO_KEEP].to_dict(orient="records")
                ncbi_df.at[index, Cons.ITS_TF_COL] = tf_info_list
            else:
                ncbi_df.at[index, Cons.ITS_TF_COL] = []

    return ncbi_df


def get_tf_target(
    tf_file: str,
    bridgedb_df: pd.DataFrame,
    filter_deg: bool,
    padj_filter: Optional[float] = 0.01,
    padj_colname: Optional[str] = None,
) -> Tuple[pd.DataFrame, dict]:
    """Add tfs and targets from tflink.

    :param tf_file: Path of the TF-Target dataset. See the links in the docstring for downloading.
    :param bridgedb_df: BridgeDb output for creating the list of gene ids to query.
    :param filter_deg: Filter the data based on DEA output, if true, makes sure the column to filter and threshold should be checked.
    :param padj_filter: The adjusted p-value threshold for filtering DEGs (default is 0.01).
    :param padj_colname: The name of the column containing adjusted p-values (default is None).
    :returns: A TFLink DataFrame and dictionary of the TFLink metadata.
    :raises ValueError: If the specified column for filtering DEGs is not found in the DataFrame.
    """
    # Dowanload TFLink dataset and metadata
    tflink_df, tflink_metadata = load_tflink_dataset(tf_file)

    # Extract the "target" values in bridgedb_df
    data_df = get_identifier_of_interest(bridgedb_df, Cons.TFLINK_GENE_INPUT_ID)

    # Filter rows where both TF and target exist in bridgedb_df
    tflink_df = tflink_df[tflink_df[Cons.TFLINK_GENE_ID_TF].isin(data_df[Cons.TARGET_COL])]
    tflink_df = tflink_df[tflink_df[Cons.TFLINK_GENE_ID_TARGET].isin(data_df[Cons.TARGET_COL])]

    # Extract the TF and the TF targets
    tf_list = list(tflink_df[Cons.TFLINK_GENE_ID_TF])
    target_list = list(tflink_df[Cons.TFLINK_GENE_ID_TARGET])

    # Add 'is_tf' and 'is_target' columns to bridgedb_df
    data_df[Cons.IS_TF_COL] = data_df[Cons.TARGET_COL].isin(tf_list)
    data_df[Cons.IS_TARGET_COL] = data_df[Cons.TARGET_COL].isin(target_list)

    # kepp only rows where the target is a DEG
    if filter_deg and padj_colname is not None:
        if f"{padj_colname}_dea" not in data_df.columns:
            raise ValueError(f"Column '{padj_colname}' not found in the DataFrame.")
        deg_df = data_df[data_df[f"{padj_colname}_dea"] <= padj_filter][Cons.TARGET_COL]
        tflink_df = tflink_df[tflink_df[Cons.TFLINK_GENE_ID_TARGET].isin(deg_df)]

    # add ensembl gene id for the TF and target
    ncbi_to_ensembl = dict(zip(data_df[Cons.TARGET_COL], data_df[Cons.IDENTIFIER_COL]))
    tflink_df[Cons.ENSEMBL_GENE_ID_TF] = tflink_df[Cons.TFLINK_GENE_ID_TF].map(ncbi_to_ensembl)
    tflink_df[Cons.ENSEMBL_GENE_ID_TARGET] = tflink_df[Cons.TFLINK_GENE_ID_TARGET].map(
        ncbi_to_ensembl
    )

    merged_df = add_target_and_tf_interaction(data_df, tflink_df, padj_filter, padj_colname)

    return merged_df, tflink_metadata
