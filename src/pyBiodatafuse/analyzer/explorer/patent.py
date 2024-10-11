# coding: utf-8

"""Python file for extracting patent data.

The module contains special functions that are server expensive and can only be performed for smaller datasets.
"""

import pandas as pd
import requests
from tqdm import tqdm

from pyBiodatafuse.utils import get_identifier_of_interest


def get_pubchem_data(bridgedb_df: pd.DataFrame) -> dict:
    """Get patent data summary from PubChem compounds.

    The output is the following: {CID: ["US: X", "EP: X", "WO: X", "Others: X"]}
    :param bridgedb_df: A dataframe with the BridgeDb or Pubchem harmonized output
    :returns: A dictionary with the PubChem Compound ID as key and the patent counts as value
    """
    # Get column of interest
    data_df = get_identifier_of_interest(bridgedb_df, "PubChem Compound")

    cid_pat_dict = {}

    # Neeed some pre-processing:
    # 1. Remove duplicates (WO-03078408-A1 against WO03078408A1, Not classic Patent offices such as AR, AU againts WO)
    # 2. Adding Granted, non-granted patent counts

    for cid in tqdm(data_df["target"]):
        patent_detail_dict = {"US": set(), "EP": set(), "WO": set(), "Others": set()}  # type: ignore[var-annotated]

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
                    patent_detail_dict["US"].add(p)
                elif p.startswith("EP"):
                    patent_detail_dict["EP"].add(p)
                elif p.startswith("WO"):
                    patent_detail_dict["WO"].add(p)
                else:
                    patent_detail_dict["Others"].add(p)

        # cid_pat_dict[cid] = [
        #     f"{k}: {len(v)}" for k, v in patent_detail_dict.items() if len(v) > 0
        # ]  # type: ignore[assignment]

    return cid_pat_dict
