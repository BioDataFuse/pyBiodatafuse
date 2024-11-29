# coding: utf-8

"""Python file for mapping identifiers using BridgeDb."""

import csv
import datetime
import json
import logging
import os
import time
from importlib import resources
from typing import List, Optional, Tuple

import pandas as pd
import requests
from pubchempy import BadRequestError, PubChemHTTPError, get_compounds, get_synonyms
from rdkit.Chem import CanonSmiles
from tqdm import tqdm

from pyBiodatafuse.constants import BRIDGEDB_ENDPOINT

logger = logging.getLogger(__name__)


def read_resource_files() -> pd.DataFrame:
    """Read the datasource file.

    :returns: a DataFrame containing the data from the datasource file
    """
    with resources.path("pyBiodatafuse.resources", "datasources.csv") as df:
        identifier_options = pd.read_csv(df)

    return identifier_options


def get_version_webservice_bridgedb() -> dict:
    """Get version of BridgeDb web service.

    :returns: a dictionary containing the version information
    :raises ValueError: if failed to retrieve data
    """
    # Set the BridgeDb API
    url = "https://webservice.bridgedb.org"
    version_response = requests.get(url=f"{url}/config")

    # Check if the request was successful (status code 200)
    if version_response.status_code == 200:
        # Initialize an empty dictionary to store the data
        bridgedb_version = {}
        # Split the response content into lines and create a CSV reader
        lines = version_response.text.strip().split("\n")
        csv_reader = csv.reader(lines, delimiter="\t")

        # Iterate over the rows in the CSV and populate the dictionary
        for row in csv_reader:
            if len(row) == 2:
                key, value = row
                bridgedb_version[key] = value
    else:
        raise ValueError(f"Failed to retrieve data. Status code: {version_response.status_code}")

    return bridgedb_version


def get_version_datasource_bridgedb(input_species: Optional[str] = None) -> List[str]:
    """Get version of BridgeDb datasource.

    :param input_species: specify the species, for now only human would be supported
    :returns: a list containing the version information
    :raises ValueError: if failed to retrieve data
    """
    if input_species is None:
        input_species = "Human"
    # Set the BridgeDb API
    url = "https://webservice.bridgedb.org"

    # Add datasource version to metadata file
    datasource_response = requests.get(url=f"{url}/{input_species}/properties")

    # Check if the request was successful (status code 200)
    if datasource_response.status_code == 200:
        datasource_version = datasource_response.text.strip().split("\n")
        datasource_version = [line.replace("\t", ": ") for line in datasource_version]
    else:
        raise ValueError(f"Failed to retrieve data. Status code: {datasource_response.status_code}")

    return datasource_version


def bridgedb_xref(
    identifiers: pd.DataFrame,
    input_species: Optional[str] = None,
    input_datasource: str = "HGNC",
    output_datasource: Optional[list] = None,
) -> Tuple[pd.DataFrame, dict]:
    """Map input list using BridgeDb.

    :param identifiers: a dataframe with one column called identifier (the output of data_loader.py)
    :param input_species: specify the species, for now only human would be supported
    :param input_datasource: type of input identifier. More details at https://www.bridgedb.org/pages/system-codes.html
    :param output_datasource: specify which type of identifiers you want to map your input identifiers.
    :returns: a DataFrame containing the mapped identifiers and dictionary of the data resource metadata.
    :raises ValueError: if the input_datasource is not provided or if the request fails
    """
    if input_species is None:
        input_species = "Human"

    if not input_datasource:
        raise ValueError("Please provide the identifier datasource, e.g. HGNC")

    if output_datasource is None or "All":
        output_datasource = [
            "Uniprot-TrEMBL",
            "NCBI Gene",
            "Ensembl",
            "HGNC Accession Number",
            "HGNC",
        ]
    else:
        assert isinstance(output_datasource, list), "output_datasource must be a list"

    data_sources = read_resource_files()
    input_source = data_sources.loc[data_sources["source"] == input_datasource, "systemCode"].iloc[
        0
    ]

    if len(identifiers) < 1:
        raise ValueError("Please provide at least one identifier datasource, e.g. HGNC")

    post_con = (
        "\n".join([f"{identifier}\t{input_source}" for identifier in identifiers["identifier"]])
        + "\n"
    )

    # Setting up the query url
    query_link = f"{BRIDGEDB_ENDPOINT}/{input_species}/xrefsBatch"

    # Record the start time
    start_time = datetime.datetime.now()

    # Getting the response to the query
    try:
        s = requests.post(url=query_link, data=post_con.encode())
        s.raise_for_status()
    except Exception as e:
        raise ValueError("Error:", e)

    # Extracting the content in the raw text format
    out = s.content.decode()
    lines = out.split("\n")

    # Record the end time
    end_time = datetime.datetime.now()

    # Processing each line and splitting values
    parsed_results = []

    for line in lines:
        if line:
            parts = line.split("\t")
            identifier = parts[0]
            identifier_source = parts[1]
            targets = parts[2].split(",")

        for target in targets:
            target_parts = target.split(":")
            target_source = target_parts[0]
            target_id = ":".join(target_parts[1:])

            parsed_results.append([identifier, identifier_source, target_id, target_source])

    # Create a DataFrame
    bridgedb = pd.DataFrame(
        parsed_results,
        columns=["identifier", "identifier.source", "target", "target.source"],
    )

    # Replace 'target.source' values with complete source names from 'data_sources'
    bridgedb["target.source"] = bridgedb["target.source"].map(
        data_sources.set_index("systemCode")["source"]
    )

    # Drop not mapped ids
    bridgedb = bridgedb.dropna(subset=["target.source"])

    # Subset based on the output_datasource
    bridgedb_subset = bridgedb[bridgedb["target.source"].isin(output_datasource)]

    bridgedb_subset = bridgedb_subset.drop_duplicates()
    identifiers.columns = [
        "{}{}".format(c, "" if c in "identifier" else "_dea") for c in identifiers.columns
    ]
    bridgedb_subset = bridgedb_subset.merge(identifiers, on="identifier")

    """Metadata details"""
    # Get the current date and time
    current_date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    # Calculate the time elapsed
    time_elapsed = str(end_time - start_time)
    # Add BridgeDb version to metadata file
    bridgedb_version = get_version_webservice_bridgedb()
    datasource_version = get_version_datasource_bridgedb()

    # Add the datasource, query, query time, and the date to metadata
    bridgedb_metadata = {
        "datasource": "BridgeDb",
        "metadata": {
            "source_version": bridgedb_version,
            "data_version": datasource_version,
        },
        "query": {
            "size": len(identifiers),
            "input_type": input_datasource,
            "time": time_elapsed,
            "date": current_date,
            "url": s.url,
            "request_string": f"{post_con.encode().decode('utf-8')}",
        },
    }

    return bridgedb_subset, bridgedb_metadata


"""PubChem helper functions."""


def check_smiles(smile: Optional[str]) -> Optional[str]:
    """Canonicalize the smiles of a compound.

    :param smile: smiles string
    :returns: canonicalized smiles string
    """
    try:
        return CanonSmiles(smile)
    except Exception:
        logger.info(f"Cannot canonicalize {smile}")
        return None


def get_cid_from_data(idx: Optional[str], idx_type: str) -> Optional[str]:
    """Get PubChem ID from any query using PubChempy.

    :param idx: identifier to query
    :param idx_type: type of identifier to query. Potential curies include : smiles, inchikey, inchi, name
    :returns: PubChem ID
    """
    if idx_type.lower() == "smiles":
        idx = check_smiles(idx)

    if not idx:
        return None

    try:
        return get_compounds(idx, idx_type.lower())[0].cid
    except BadRequestError:
        logger.info(f"Issue with {idx}")
        return None

    except IndexError:
        logger.info(f"Issue with {idx}")
        return None


def get_cid_from_pugrest(idx: Optional[str], idx_type: str) -> Optional[str]:
    """Get PubChem ID from any query throung Pubchem PUGREST.

    :param idx: identifier to query
    :param idx_type: type of identifier to query. Potential curies include : smiles, inchikey, inchi, name
    :returns: PubChem ID
    """
    if idx_type.lower() == "smiles":
        idx = check_smiles(idx)

    if not idx:
        return None

    cid_data = requests.get(
        f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/{idx_type}/{idx}/property/Title/JSON"
    ).json()

    if "Fault" in cid_data:
        logger.info(f"Issue with {idx}")
        return None

    cidx = cid_data["PropertyTable"]["Properties"][0]["CID"]
    if "." in str(cidx):
        return str(cidx).split(".")[0]
    return str(cidx)


def pubchem_xref(
    identifiers: list, identifier_type: str = "name", cache_res: bool = False
) -> Tuple[pd.DataFrame, dict]:
    """Map chemical names or smiles or inchikeys to PubChem identifier.

    :param identifiers: a list of identifiers to query
    :param identifier_type: type of identifier to query. Potential curies include : smiles, inchikey, inchi, name
    :param cache_res: whether to cache the results
    :raises ValueError: if the input_datasource is not provided or if the request fails
    :returns: a DataFrame containing the mapped identifiers and dictionary of the data resource metadata.
    """
    if len(identifiers) < 1:
        raise ValueError("Please provide at least one input.")

    # Record the start time
    start_time = datetime.datetime.now()

    # Getting the response to the query
    cid_data = []
    c = 0

    if cache_res:
        if os.path.exists("pubchem_cache_results.json"):
            with open("pubchem_cache_results.json", "r") as f:
                cache_results = json.load(f)
        else:
            cache_results = {}
    else:
        cache_results = {}

    c = 0
    for idx in tqdm(identifiers, desc="Mapping PubChem"):
        if idx in cache_results:
            cid = cache_results[idx]
        else:
            c += 1
            if c == 100:
                if cache_res:
                    with open("pubchem_cache_results.json", "w") as f:
                        json.dump(cache_results, f)
                time.sleep(5)
                c = 0

            cid = get_cid_from_pugrest(idx, identifier_type)
            cache_results[idx] = cid

        cid_data.append(
            {
                "identifier": idx,
                "identifier.source": identifier_type,
                "target": f"pubchem.compound:{cid}" if cid is not None else None,
                "target.source": "PubChem Compound",
            }
        )

    # Record the end time
    end_time = datetime.datetime.now()

    pubchem_df = pd.DataFrame(cid_data)
    pubchem_df = pubchem_df.drop_duplicates()

    """Metadata details"""
    # Get the current date and time
    current_date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    # Calculate the time elapsed
    time_elapsed = str(end_time - start_time)
    # Add package version to metadata file
    stable_package_version = "1.0.4"  # Stable version for PubChemPy

    # Add the datasource, query, query time, and the date to metadata
    pubchem_metadata = {
        "datasource": "Pubchem python client",
        "metadata": {
            "package": "PubChemPy",
            "data_version": stable_package_version,
        },
        "query": {
            "size": len(identifiers),
            "input_type": identifier_type,
            "time": time_elapsed,
            "date": current_date,
        },
    }

    return pubchem_df, pubchem_metadata


def cid2chembl(cids: list) -> dict:
    """Map Pubchem CIDs to ChEMBL identifier.

    :param cids: a list of CIDs identifiers to query
    :raises ValueError: if the input_datasource is not provided or if the request fails
    :returns: a dictonary of ChEMBL mapped to CID identifiers and dictionary of the data resource metadata.
    """
    if len(cids) < 1:
        raise ValueError("Please provide at least one input.")

    # Getting the response to the query
    chembl_data = {}  # ChEMBL ids as keys and PubChem ids as values
    for pubchem_idx in cids:
        try:
            other_idenfitiers = get_synonyms(identifier=pubchem_idx)
        except (PubChemHTTPError, BadRequestError):  # too many request
            time.sleep(3)
            try:
                other_idenfitiers = get_synonyms(identifier=pubchem_idx)
            except BadRequestError:  # incorrect pubchem id
                continue

        if len(other_idenfitiers) < 1:
            continue

        other_idenfitiers = other_idenfitiers[0]

        for idx in other_idenfitiers["Synonym"]:
            if idx.startswith("CHEMBL"):
                chembl_data[idx] = pubchem_idx
                break

    return chembl_data
