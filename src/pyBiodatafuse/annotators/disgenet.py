# coding: utf-8

"""Python file for queriying DisGeNet database (https://www.disgenet.org/home/)."""

import datetime
import json
import logging
import time
import warnings
from typing import Any, Dict, List, Set, Tuple

import pandas as pd
import requests
from tqdm import tqdm
from urllib3 import disable_warnings

from pyBiodatafuse.constants import (
    DISGENET,
    DISGENET_DISEASE_COL,
    DISGENET_DISEASE_OUTPUT_DICT,
    DISGENET_ENDPOINT,
    DISGENET_GENE_INPUT_ID,
)
from pyBiodatafuse.utils import (
    check_columns_against_constants,
    collapse_data_sources,
    get_identifier_of_interest,
)

logger = logging.getLogger("disgenet")
disable_warnings()


def check_endpoint_disgenet(api_key: str) -> bool:
    """Check the availability of the DisGeNET API.

    :param api_key: DisGeNET API key (more details can be found at https://disgenet.com/plans)
    :returns: True if the endpoint is available, False otherwise.
    """
    # Set HTTP headers
    httpheadersdict = {}
    httpheadersdict["Authorization"] = api_key
    httpheadersdict["accept"] = "application/json"
    # Set the DisGeNET API
    s = requests.Session()
    # Get version
    response = s.get("https://api.disgenet.com/api/v1/public/version", headers=httpheadersdict)
    # Check if API is down
    if response.json()["status"] == "OK":
        return True
    else:
        return False


def get_version_disgenet(api_key: str) -> dict:
    """Get version of DisGeNET API.

    :param api_key: DisGeNET API key (more details can be found at https://disgenet.com/plans)
    :returns: a dictionary containing the version information
    """
    # Set HTTP headers
    httpheadersdict = {}
    httpheadersdict["Authorization"] = api_key
    httpheadersdict["accept"] = "application/json"
    # Set the DisGeNET API
    s = requests.Session()
    # Get version
    version_response = s.get(
        "https://api.disgenet.com/api/v1/public/version", headers=httpheadersdict
    )
    disgenet_version = version_response.json()["payload"]

    return disgenet_version


def _format_dis_identifiers(row, namespace: str) -> List[str]:
    """Format the disease identifiers.

    :param row: List of disease identifiers
    :param namespace: Namespace to be added to the identifiers
    :returns: a list of formatted disease identifiers
    """
    new_vals = []  # type: List[Any]

    for val in row:
        if pd.isna(val):
            new_vals.append(None)
            continue

        if val == "":
            new_vals.append(None)
            continue

        t = []
        for v in val.split(", "):
            p = v.split("_")[-1]
            t.append(f"{namespace}:{p}")
        new_vals.append(", ".join(t))
    return new_vals


def _format_disgenet_output(intermediate_df: pd.DataFrame) -> pd.DataFrame:
    """Format the DisGeNET output.

    :param intermediate_df: DataFrame containing the DisGeNET output
    :returns: a DataFrame containing the formatted DisGeNET output
    """
    # extract disease identifiers from diseaseVocabularies column
    # Initialize dictionaries to store the columns
    source_types: Set[str] = set()
    # Process the 'diseaseVocabularies' column
    for entry in intermediate_df["diseaseVocabularies"]:
        for item in entry:
            if isinstance(item, str):
                # Remove everything after '_'
                prefix = item.split("_")[0]
                # Add to the set
                source_types.add(prefix)
    # Convert set to list
    source_type_list: List[str] = list(source_types)

    # Add new columns for each identifier type and initialize with empty lists
    for source in source_type_list:
        intermediate_df[source] = None

    # Populate the new columns with identifiers
    for index, entry in intermediate_df.iterrows():
        vocab_list = entry["diseaseVocabularies"]
        # Create a dictionary to hold identifiers by type
        identifiers_by_type: Dict[str, List[str]] = {source: [] for source in source_type_list}
        for item in vocab_list:
            if isinstance(item, str):
                # Extract the type and identifier
                parts = item.split("_")
                if len(parts) > 1:
                    source_type = parts[0]
                    if source_type in identifiers_by_type:
                        identifiers_by_type[source_type].append(item)
        # Populate the DataFrame with the collected identifiers
        for source in source_type_list:
            # Join the identifiers with comma and format as a list
            intermediate_df.at[index, source] = ", ".join(identifiers_by_type[source])

    intermediate_df.rename(
        columns={
            "geneNcbiID": "target",
            "diseaseName": "disease_name",
            "diseaseType": "disease_type",
            "diseaseUMLSCUI": "disease_umlscui",
        },
        inplace=True,
    )
    intermediate_df["target"] = intermediate_df["target"].values.astype(str)

    missing_cols = [
        col for col in DISGENET_DISEASE_OUTPUT_DICT.keys() if col not in intermediate_df.columns
    ]
    for col in missing_cols:
        intermediate_df[col] = None

    selected_columns = [
        # "geneDSI",
        # "geneDPI",
        # "genepLI",
        # "geneNcbiType",
        # "geneProteinClassIDs",
        # "geneProteinClassNames",
        # "diseaseVocabularies",
        "target",
        *DISGENET_DISEASE_OUTPUT_DICT.keys(),
    ]
    intermediate_df = intermediate_df[selected_columns]

    # Adding namespace prefixes to the identifiers
    intermediate_df["HPO"] = _format_dis_identifiers(intermediate_df["HPO"], namespace="HPO")
    intermediate_df["NCI"] = _format_dis_identifiers(intermediate_df["NCI"], namespace="NCI")
    intermediate_df["OMIM"] = _format_dis_identifiers(intermediate_df["OMIM"], namespace="MIM")
    intermediate_df["MONDO"] = _format_dis_identifiers(intermediate_df["MONDO"], namespace="MONDO")
    intermediate_df["ORDO"] = _format_dis_identifiers(intermediate_df["ORDO"], namespace="ORDO")
    intermediate_df["EFO"] = _format_dis_identifiers(intermediate_df["EFO"], namespace="EFO")
    intermediate_df["DO"] = _format_dis_identifiers(intermediate_df["DO"], namespace="DOID")
    intermediate_df["MESH"] = _format_dis_identifiers(intermediate_df["MESH"], namespace="MESH")
    intermediate_df["UMLS"] = _format_dis_identifiers(intermediate_df["UMLS"], namespace="UMLS")
    return intermediate_df


def get_gene_disease(api_key: str, bridgedb_df: pd.DataFrame) -> Tuple[pd.DataFrame, dict]:
    """Query gene-disease associations from DisGeNET.

    :param api_key: DisGeNET API key (more details can be found at https://disgenet.com/plans)
    :param bridgedb_df: BridgeDb output for creating the list of gene ids to query.
    :returns: a DataFrame containing the DisGeNET output and dictionary of the DisGeNET metadata.
    """
    # Check if the DisGeNET API is available
    api_available = check_endpoint_disgenet(api_key)

    if not api_available:
        warnings.warn(
            f"{DISGENET} endpoint is not available. Unable to retrieve data.", stacklevel=2
        )
        return pd.DataFrame(), {}

    # Add the API key to the requests headers
    httpheadersdict = {}
    httpheadersdict["Authorization"] = api_key
    httpheadersdict["accept"] = "application/json"
    # Set the DisGeNET API
    s = requests.Session()

    # Extract the "target" values
    data_df = get_identifier_of_interest(bridgedb_df, DISGENET_GENE_INPUT_ID)

    disgenet_output = []

    # Specify query parameters by means of a dictionary
    params = {}
    params["format"] = "json"
    params["source"] = "CURATED"

    # Record the start time
    disgenet_version = get_version_disgenet(api_key)
    start_time = datetime.datetime.now()

    c = 0
    for gene in tqdm(data_df["target"], desc="Querying DisGeNET"):
        # Retrieve disease associated to gene with NCBI ID
        params["gene_ncbi_id"] = gene
        params["page_number"] = str(0)

        c += 1
        # Get all the diseases associated with genes for the current chunk
        gda_response = s.get(
            DISGENET_ENDPOINT, params=params, headers=httpheadersdict, verify=False
        )

        # If the status code of gda_response is 429, it means you have reached one of your query limits
        # You can retrieve the time you need to wait until doing a new query in the response headers
        if gda_response.ok:
            # Parse response content in JSON format since we set 'accept:application/json' as HTTP header
            response_parsed = json.loads(gda_response.text)
            disgenet_output.extend(response_parsed["payload"])
        elif gda_response.status_code == 429:
            while gda_response.ok is False:
                time.sleep(int(gda_response.headers["x-rate-limit-retry-after-seconds"]))

                # Repeat your query
                gda_response = s.get(
                    DISGENET_ENDPOINT,
                    params=params,
                    headers=httpheadersdict,
                    verify=False,
                    timeout=None,
                )
                if gda_response.ok is True:
                    break

        if c == 100:
            time.sleep(20)
            c = 0
    # Record the end time
    end_time = datetime.datetime.now()

    """Metadata details"""
    # Get the current date and time
    current_date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    # Calculate the time elapsed
    time_elapsed = str(end_time - start_time)

    # Add version, datasource, query, query time, and the date to metadata
    disgenet_metadata: Dict[str, Any] = {
        "datasource": DISGENET,
        "metadata": disgenet_version,
        "query": {
            "size": len(data_df["target"].drop_duplicates()),
            "input_type": DISGENET_GENE_INPUT_ID,
            "time": time_elapsed,
            "date": current_date,
            "url": DISGENET_ENDPOINT,
        },
    }

    # Organize the annotation results as an array of dictionaries
    intermediate_df = pd.DataFrame(disgenet_output)
    if "geneNcbiID" not in intermediate_df:
        warnings.warn(
            f"There is no annotation for your input list in {DISGENET}.",
            stacklevel=2,
        )
        return pd.DataFrame(), disgenet_metadata

    # Format the DisGeNET output
    intermediate_df = _format_disgenet_output(intermediate_df)

    # Check if all keys in df match the keys in OUTPUT_DICT
    check_columns_against_constants(
        data_df=intermediate_df,
        output_dict=DISGENET_DISEASE_OUTPUT_DICT,
        check_values_in=["HPO", "NCI", "OMIM", "MONDO", "ORDO", "EFO", "DO", "MESH", "UMLS"],
    )

    merged_df = collapse_data_sources(
        data_df=bridgedb_df,
        source_namespace=DISGENET_GENE_INPUT_ID,
        target_df=intermediate_df,
        common_cols=["target"],
        target_specific_cols=list(DISGENET_DISEASE_OUTPUT_DICT.keys()),
        col_name=DISGENET_DISEASE_COL,
    )

    """Update metadata"""
    # Calculate the number of new nodes
    num_new_nodes = intermediate_df["disease_name"].nunique()
    # Calculate the number of new edges
    num_new_edges = intermediate_df.drop_duplicates(subset=["target", "disease_name"]).shape[0]

    # Check the intermediate_df
    if num_new_edges != len(intermediate_df):
        warnings.warn(
            f"The intermediate_df in {DISGENET} annotatur should be checked, please create an issue on https://github.com/BioDataFuse/pyBiodatafuse/issues/.",
            stacklevel=2,
        )

    # Add the number of new nodes and edges to metadata
    disgenet_metadata["query"]["number_of_added_nodes"] = num_new_nodes
    disgenet_metadata["query"]["number_of_added_edges"] = num_new_edges

    return merged_df, disgenet_metadata
