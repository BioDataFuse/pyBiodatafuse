# coding: utf-8

"""Python file for extracting patent data from PubChem.

The module contains special functions that are server expensive and can only be performed for smaller datasets.
"""
from typing import Dict

import pandas as pd
import requests
from tqdm import tqdm

from pyBiodatafuse.utils import get_identifier_of_interest


def get_patent_data(data_input: pd.DataFrame) -> Dict[str, dict]:
    """Get patent data summary from PubChem.

    :param data_input: A dataframe with the BridgeDb or Pubchem harmonized output
    :returns: A dictionary with the PubChem Compound ID as key and the patent counts as value
    The output is the following: {CID: ["US: X", "EP: X", "WO: X", "Others: X"]}
    """
    # Get column of interest
    data_df = get_identifier_of_interest(data_input, "PubChem Compound")

    cid_pat_dict = {}

    # Neeed some pre-processing:
    # 1. Remove duplicates (WO-03078408-A1 against WO03078408A1, Not classic Patent offices such as AR, AU againts WO)
    # 2. Adding Granted, non-granted patent counts

    for cid in tqdm(data_df["target"]):
        patent_counter_dict = {"US": set(), "EP": set(), "WO": set(), "Others": set()}

        patent_dict = requests.get(
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/xrefs/PatentID/JSON"
        ).json()

        if "Fault" in patent_dict.keys():  # No patents found
            cid_pat_dict[cid] = ""
        else:
            patents = patent_dict["InformationList"]["Information"][0]["PatentID"]

            for patent in patents:
                p = patent.replace("-", "")
                if p.startswith("US"):
                    patent_counter_dict["US"].add(p)
                elif p.startswith("EP"):
                    patent_counter_dict["EP"].add(p)
                elif p.startswith("WO"):
                    patent_counter_dict["WO"].add(p)
                else:
                    patent_counter_dict["Others"].add(p)

        patent_counter_dict = [
            f"{k}: {len(v)}" for k, v in patent_counter_dict.items() if len(v) > 0
        ]

        cid_pat_dict[cid] = patent_counter_dict

    return cid_pat_dict


def _process_data_for_plot(data_dict: dict) -> pd.DataFrame:
    """Process data to map to plotting template."""
    data = []

    for cidx, pat_counter in data_dict.items():
        for pat in pat_counter:
            p, c = pat.split(": ")
            data.append(
                {
                    "cid": cidx,
                    "label": p,
                    "value": c,
                }
            )

    data_df = pd.DataFrame(data).set_index("cid")

    return data_df
