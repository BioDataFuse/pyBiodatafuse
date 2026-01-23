# coding: utf-8

"""Python file for mapping identifiers using BridgeDb."""

import csv
import datetime
import json
import logging
import os
import re
import time
from importlib import resources
from typing import List, Literal, Optional, Tuple

import pandas as pd
import requests
from pubchempy import BadRequestError, PubChemHTTPError, get_compounds, get_synonyms
from rdkit.Chem import CanonSmiles
from tqdm import tqdm

import pyBiodatafuse.constants as Cons

logger = logging.getLogger(__name__)


def read_datasource_file() -> pd.DataFrame:
    """Read the datasource file.

    :returns: a DataFrame containing the data from the datasource file
    """
    with resources.path("pyBiodatafuse.resources", "datasources.csv") as df:
        identifier_options = pd.read_csv(df)

    return identifier_options


def match_input_datasource(identifiers) -> str:
    """Check if the input identifiers match the datasource.

    This function attempts to match the provided identifiers against known patterns
    in the datasource file and returns the corresponding data source.

    :param identifiers: a pandas DataFrame containing the identifiers to be matched
    :returns: data source
    :raises ValueError: if the identifiers series is empty, no match is found, or multiple matches are found
    """
    if identifiers.empty:
        raise ValueError("The identifiers series is empty.")

    with resources.path("pyBiodatafuse.resources", "datasources.csv") as df_file:
        datasources = pd.read_csv(df_file)

    matched_sources = set()
    for identifier in identifiers:
        match_found = False
        for _, row in datasources.iterrows():
            pattern = (
                str(row["pattern"]) if pd.notna(row["pattern"]) else None
            )  # Handle NaN patterns
            if not pattern:
                continue  # Skip rows with invalid patterns
            if "ENS" in identifier:
                return "Ensembl"
            try:
                if re.fullmatch(pattern, identifier):
                    if pattern not in [r"^\d+$", r"^\S+$"]:
                        matched_sources.add(row["source"])
                        match_found = True
            except re.error as e:
                logger.warning(f"Invalid regex pattern '{pattern}': {e}")
                continue  # Skip invalid regex patterns
        if not match_found:
            raise ValueError(f"Identifier '{identifier}' does not match any known pattern.")

    if len(matched_sources) > 1:
        logger.info(f"Matched data sources: {', '.join(matched_sources)}")
        raise ValueError(
            f"Multiple data sources match the provided identifiers (e.g., {identifier}): {', '.join(matched_sources)}. "
            "Please specify the datasource explicitly using `input_datasource`."
        )

    return matched_sources.pop()


def get_version_webservice_bridgedb() -> dict:
    """Get version of BridgeDb web service.

    :returns: a dictionary containing the version information
    :raises ValueError: if failed to retrieve data
    """
    # Set the BridgeDb API
    version_response = requests.get(url=f"{Cons.BRIDGEDB_ENDPOINT}/config")

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

    # Add datasource version to metadata file
    datasource_response = requests.get(url=f"{Cons.BRIDGEDB_ENDPOINT}/{input_species}/properties")

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
    output_datasource: Optional[list] = None,
    input_datasource: Literal[
        "Ensembl",
        "NCBI Gene",
        "HGNC",
        "HGNC Accession Number",
        "MGI",
        "miRBase mature sequence",
        "miRBase Sequence",
        "OMIM",
        "RefSeq",
        "Rfam",
        "RGD",
        "SGD",
        "UCSC Genome Browser",
        "NCBI Protein",
        "PDB",
        "Pfam",
        "Uniprot-TrEMBL",
        "Uniprot-SwissProt",
        "Affy",
        "Agilent",
        "Illumina",
        "Gene Ontology",
        "CAS",
        "ChEBI",
        "ChemSpider",
        "ChEMBL compound",
        "DrugBank",
        "HMDB",
        "Guide to Pharmacology Ligand ID",
        "InChIKey",
        "KEGG Compound",
        "KEGG Drug",
        "KEGG Glycan",
        "LIPID MAPS",
        "LipidBank",
        "PharmGKB Drug",
        "PubChem Compound",
        "PubChem Substance",
        "SwissLipids",
        "TTD Drug",
        "Wikidata",
        "Wikipedia",
    ] = "HGNC",
) -> Tuple[pd.DataFrame, dict]:
    """
    Map input identifiers using BridgeDb.

    :param identifiers: A pandas DataFrame with one column named 'identifier'.
    :param input_species: Optional species name. Only 'Homo sapiens' is currently supported.
    :param input_datasource: The type of identifier in the input DataFrame. Expected formats by datasource:
        - "HGNC": e.g. "TP53"
        - "HGNC Accession Number": e.g. "HGNC:11998"
        - "Ensembl": e.g. "ENSG00000141510"
        - "NCBI Gene": e.g. "7157"
        - "MGI": e.g. "MGI:104874"
        - "miRBase mature sequence": e.g. "hsa-miR-21-5p"
        - "miRBase Sequence": e.g. "MI0000077"
        - "OMIM": e.g. "191170"
        - "RefSeq": e.g. "NM_000546"
        - "Rfam": e.g. "RF00001"
        - "RGD": e.g. "RGD:620474"
        - "SGD": e.g. "YAL001C"
        - "UCSC Genome Browser": e.g. "uc001aaa.3"
        - "NCBI Protein": e.g. "NP_000537"
        - "PDB": e.g. "1TUP"
        - "Pfam": e.g. "PF00069"
        - "Uniprot-SwissProt": e.g. "P04637"
        - "Uniprot-TrEMBL": e.g. "Q9H0H5"
        - "Affy": e.g. "202763_at"
        - "Agilent": e.g. "A_23_P61180"
        - "Illumina": e.g. "ILMN_1803030"
        - "Gene Ontology": e.g. "GO:0006915"
        - "CAS": e.g. "50-00-0"
        - "ChEBI": e.g. "CHEBI:15377"
        - "ChemSpider": e.g. "5798"
        - "ChEMBL compound": e.g. "CHEMBL25"
        - "DrugBank": e.g. "DB01050"
        - "HMDB": e.g. "HMDB0000122"
        - "Guide to Pharmacology Ligand ID": e.g. "1234"
        - "InChIKey": e.g. "BSYNRYMUTXBXSQ-UHFFFAOYSA-N"
        - "KEGG Compound": e.g. "C00031"
        - "KEGG Drug": e.g. "D00001"
        - "KEGG Glycan": e.g. "G00001"
        - "LIPID MAPS": e.g. "LMFA01010001"
        - "LipidBank": e.g. "LBID0001"
        - "PharmGKB Drug": e.g. "PA449053"
        - "PubChem Compound": e.g. "2244"
        - "PubChem Substance": e.g. "12345678"
        - "SwissLipids": e.g. "SLM:000000001"
        - "TTD Drug": e.g. "D000001"
        - "Wikidata": e.g. "Q18216"
        - "Wikipedia": e.g. "Aspirin"
    :param output_datasource: Optional list of identifier types to map to.
    :returns: Tuple of:
        - DataFrame with mapped identifiers.
        - Dictionary of data resource metadata.
    :raises ValueError: If required inputs are missing or the mapping fails.
    """
    if input_species is None:
        input_species = "Human"
    data_sources = read_datasource_file()
    input_source = data_sources.loc[
        data_sources[Cons.SOURCE_COL] == input_datasource, "systemCode"
    ].iloc[0]
    input_type = data_sources.loc[data_sources[Cons.SOURCE_COL] == input_datasource, "type"].iloc[0]
    if output_datasource is None or "All":
        output_datasource = data_sources[data_sources["type"] == input_type]["source"].tolist()
    else:
        assert isinstance(output_datasource, list), "output_datasource must be a list"

    if len(identifiers) < 1:
        raise ValueError("Please provide at least one identifier datasource, e.g. HGNC")

    post_con = (
        "\n".join([f"{identifier}\t{input_source}" for identifier in identifiers["identifier"]])
        + "\n"
    )

    # Setting up the query url
    query_link = f"{Cons.BRIDGEDB_ENDPOINT}/{input_species}/xrefsBatch"

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

            parsed_results.append(
                [identifier, identifier_source, target_id, target_source]
            )  # Create a DataFrame
    bridgedb = pd.DataFrame(
        parsed_results,
        columns=[
            Cons.IDENTIFIER_COL,
            Cons.IDENTIFIER_SOURCE_COL,
            Cons.TARGET_COL,
            Cons.TARGET_SOURCE_COL,
        ],
    )

    # Replace 'target.source' values with complete source names from 'data_sources'
    bridgedb[Cons.TARGET_SOURCE_COL] = bridgedb[Cons.TARGET_SOURCE_COL].map(
        data_sources.set_index("systemCode")[Cons.SOURCE_COL]
    )
    # Drop not mapped ids
    bridgedb = bridgedb.dropna(subset=[Cons.TARGET_SOURCE_COL])
    # Subset based on the output_datasource
    bridgedb_subset = bridgedb[bridgedb[Cons.TARGET_SOURCE_COL].isin(output_datasource)]
    bridgedb_subset = bridgedb_subset.drop_duplicates()
    identifiers.columns = [
        "{}{}".format(c, "" if c in "identifier" else "_dea") for c in identifiers.columns
    ]
    bridgedb_subset = bridgedb_subset.merge(identifiers, on=Cons.IDENTIFIER_COL)
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
        "datasource": Cons.BRIDGEDB,
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
    if idx_type.lower() == Cons.SMILES.lower():
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
    if idx_type.lower() == Cons.SMILES.lower():
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
                Cons.IDENTIFIER_COL: idx,
                Cons.IDENTIFIER_SOURCE_COL: identifier_type,
                Cons.TARGET_COL: f"{Cons.PUBCHEM_COMPOUND_CID}:{cid}" if cid is not None else None,
                Cons.TARGET_SOURCE_COL: Cons.PUBCHEM_COMPOUND,
            }
        )

    if cache_res:
        with open("pubchem_cache_results.json", "w") as f:
            json.dump(cache_results, f)

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
        "datasource": Cons.PUBCHEM,
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
            if idx.startswith(Cons.CHEMBL):
                chembl_data[idx] = pubchem_idx
                break

    return chembl_data
