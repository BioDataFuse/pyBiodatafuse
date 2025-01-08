# coding: utf-8

"""Python file for TFLink ETL process (https://tflink.net/download/)."""

import gzip
import os
from datetime import datetime
from typing import Tuple
from xmlrpc.client import boolean

import pandas as pd
import requests

from pyBiodatafuse.constants import TFLINK, TFLINK_DOWNLOAD_URL, TFLINK_GENE_INPUT_ID
from pyBiodatafuse.utils import get_identifier_of_interest


def download_tflink_dataset(tf_file: str, filename: str) -> Tuple[pd.DataFrame, dict]:
    """Downloads, saves and reads a tfLink dataset.

    :param tf_file: The TF-Target dataset to download. Human "TFLink_Homo_sapiens_interactions_All_simpleFormat_v1.0.tsv.gz"
    :param filename: The local file path to save the downloaded dataset.
    :returns: A TFLink DataFrame and dictionary of the TFLink metadata.
    """
    # Dowonload the TF-Target dataset
    if not os.path.exists(filename):
        url = f"{TFLINK_DOWNLOAD_URL}/{tf_file}"
        response = requests.get(url)
        try:
            response.raise_for_status()
        except requests.HTTPError as e:
            print(f"Failed to download file. HTTP Error: {e}")
        else:
            with open(filename, "wb") as file:
                file.write(response.content)
    else:
        with gzip.open(filename, "rt") as f:
            tflink_df = pd.read_csv(f, sep="\t")

    # Add version
    tflink_metadata = {
        "datasource": TFLINK,
        "metadata": {
            "download date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "download link": f"{TFLINK_DOWNLOAD_URL}/{tf_file}",
        },
    }

    return tflink_df, tflink_metadata


def add_target_and_tf_interaction(ncbi_df: pd.DataFrame, tflink_df: pd.DataFrame) -> pd.DataFrame:
    """Add targets and TFs to each row (gene).

    :param ncbi_df: BridgeDb output with ncbi id as target source.
    :param tflink_df: The TF-Target dataset.
    :returns: ncbi_df with targets and TFs in each row.
    """
    ncbi_df["its_target"] = None
    for index, row in ncbi_df.iterrows():
        if row["is_tf"]:
            targets = tflink_df[tflink_df["NCBI.GeneID.TF"] == row["target"]]
            if not targets.empty:
                target_info_list = targets[
                    [
                        "NCBI.GeneID.Target",
                        "Ensembl.GeneID.Target",
                        "Name.Target",
                        "UniprotID.Target",
                        "Target.TFLink.ortho",
                        "Target.nonTFLink.ortho",
                        "Detection.method",
                        "PubmedID",
                        "Source.database",
                        "Small-scale.evidence",
                    ]
                ].to_dict(orient="records")
                ncbi_df.at[index, "its_target"] = target_info_list
            else:
                ncbi_df.at[index, "its_target"] = []

    ncbi_df["its_tf"] = None
    for index, row in ncbi_df.iterrows():
        if row["is_target"] and row["padj_dea"] <= 0.01:
            tf = tflink_df[tflink_df["NCBI.GeneID.Target"] == row["target"]]
            if not tf.empty:
                tf_info_list = tf[
                    [
                        "NCBI.GeneID.TF",
                        "Ensembl.GeneID.TF",
                        "Name.TF",
                        "UniprotID.TF",
                        "TF.TFLink.ortho",
                        "TF.nonTFLink.ortho",
                        "Detection.method",
                        "PubmedID",
                        "Source.database",
                        "Small-scale.evidence",
                    ]
                ].to_dict(orient="records")
                ncbi_df.at[index, "its_tf"] = tf_info_list
            else:
                # if row["padj_dea"] <= 0.01:
                #     ncbi_df.at[index, "its_tf"] = None
                # else:
                ncbi_df.at[index, "its_tf"] = []

    return ncbi_df


def get_tf_target(tf_file: str, filename: str, filter_deg: boolean, bridgedb_df: pd.DataFrame):
    """Add tfs and targets from tflink.

    :param tf_file: The TF-Target dataset to download. Human "TFLink_Homo_sapiens_interactions_All_simpleFormat_v1.0.tsv.gz"
    :param filename: The local file path to save the downloaded dataset.
    :param filter_deg: Filter the data based on DEA output, if true, makes sure the column to filter and thereshold should be checked
    :param bridgedb_df: BridgeDb output for creating the list of gene ids to query.
    :returns: A TFLink DataFrame and dictionary of the TFLink metadata.
    """
    # Dowanload TFLink dataset and metadata
    tflink_df, tflink_metadata = download_tflink_dataset(tf_file, filename)

    # Extract the "target" values in bridgedb_df
    data_df = get_identifier_of_interest(bridgedb_df, TFLINK_GENE_INPUT_ID)

    # Filter rows where both TF and target exist in bridgedb_df
    tflink_df = tflink_df[tflink_df["NCBI.GeneID.TF"].isin(data_df["target"])]
    tflink_df = tflink_df[tflink_df["NCBI.GeneID.Target"].isin(data_df["target"])]

    # Extract the TF and the TF targets
    tf_list = list(tflink_df["NCBI.GeneID.TF"])
    target_list = list(tflink_df["NCBI.GeneID.Target"])

    # Add 'is_tf' and 'is_target' columns to bridgedb_df
    data_df["is_tf"] = data_df["target"].isin(tf_list)
    data_df["is_target"] = data_df["target"].isin(target_list)

    # kepp only rows where the target is a DEG
    if filter_deg:
        tflink_df = tflink_df[
            tflink_df["NCBI.GeneID.Target"].isin(data_df[data_df["padj_dea"] <= 0.01]["target"])
        ]

    # add ensembl gene id for the TF and target
    ncbi_to_ensembl = dict(zip(data_df["target"], data_df["identifier"]))
    tflink_df["Ensembl.GeneID.TF"] = tflink_df["NCBI.GeneID.TF"].map(ncbi_to_ensembl)
    tflink_df["Ensembl.GeneID.Target"] = tflink_df["NCBI.GeneID.Target"].map(ncbi_to_ensembl)

    merged_df = add_target_and_tf_interaction(data_df, tflink_df)

    return merged_df, tflink_metadata
