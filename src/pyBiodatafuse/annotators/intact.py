import datetime
import warnings
import requests
import pandas as pd
from pyBiodatafuse.constants import INTACT, INTACT_GENE_INPUT_ID
from pyBiodatafuse.utils import get_identifier_of_interest

INTACT_ENDPOINT = "https://www.ebi.ac.uk/intact/ws/interaction"

def check_endpoint_intact() -> bool:
    """Check if the IntAct API is reachable by making a test request."""
    test_gene = "P53"  # Example known gene ID
    response = requests.get(f"{INTACT_ENDPOINT}/findInteractions/{test_gene}")
    return response.status_code == 200

def get_intact_interactions(gene_id: str):
    """Retrieve protein interactions for a given gene from IntAct."""
    response = requests.get(f"{INTACT_ENDPOINT}/findInteractions/{gene_id}")

    if response.status_code == 200:
        data = response.json()
        if not data.get("content"):  # Check if response is empty
            return []
        
        interactions = [
            {
                "interaction_id": item.get("ac", "N/A"),
                "interactor_id_A": item.get("acA", ""),
                "interactor_id_B": item.get("acB", ""),
                "confidence_values": item.get("confidenceValues", []),
                "biological_role_A": item.get("biologicalRoleA", ""),
                "biological_role_B": item.get("biologicalRoleB", ""),
                "type": item.get("type", ""),
                "stoichiometry_A": item.get("stoichiometryA", ""),
                "stoichiometry_B": item.get("stoichiometryB", ""),
                "detection_method": item.get("detectionMethod", ""),
                "host_organism": item.get("hostOrganism", ""),
                "interactor_A_name": item.get("intactNameA", ""),
                "interactor_B_name": item.get("intactNameB", ""),
                "interactor_A_species": item.get("speciesA", ""),
                "interactor_B_species": item.get("speciesB", ""),
            }
            for item in data["content"]
        ]
        return interactions
    
    return []  # Return empty list if API call fails

def get_interactions(bridgedb_df):
    """Annotate genes with interaction data from IntAct."""
    api_available = check_endpoint_intact()
    if not api_available:
        warnings.warn("IntAct API endpoint is unavailable. Cannot retrieve data.", stacklevel=2)
        return pd.DataFrame(), {}

    # Record the start time
    start_time = datetime.datetime.now()

    # Get identifiers of interest
    data_df = get_identifier_of_interest(bridgedb_df, INTACT_GENE_INPUT_ID)

    if isinstance(data_df, tuple):  # Ensure it's a DataFrame
        data_df = data_df[0]

    data_df = data_df.reset_index(drop=True)
    gene_list = list(set(data_df["target"].tolist()))

    # Retrieve interactions for each gene
    data_df["IntAct_interactions"] = data_df["target"].apply(get_intact_interactions)

    # Record the end time
    end_time = datetime.datetime.now()

    """Metadata details"""
    # Get the current date and time
    current_date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    # Calculate the time elapsed
    time_elapsed = str(end_time - start_time)
    # Calculate new edges
    num_new_edges = data_df.shape[0]

    # Metadata details
    intact_metadata = {
        "datasource": INTACT,
    #    "metadata": {"source_version": intact_version},
        "query": {
            "size": len(gene_list),
            "input_type": INTACT_GENE_INPUT_ID,
            "number_of_added_edges": num_new_edges,
            "time": time_elapsed,
            "date": current_date,
            "url": INTACT_ENDPOINT,
        },
    }

    return data_df, intact_metadata
