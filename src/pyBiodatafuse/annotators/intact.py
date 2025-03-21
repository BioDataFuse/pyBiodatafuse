import datetime
import warnings
import numpy as np
import pandas as pd
import requests

INTACT_ENDPOINT = "www.ebi.ac.uk/intact/ws/interaction"
INTACT_INPUT_GENE_ID = "Ensembl"

def check_endpoint_intact() -> bool:
    """Check if the IntAct API endpoint is available."""
    response = requests.get(f"{INTACT_ENDPOINT}/ping")
    return response.status_code == 200

def get_intact_interactions(gene_id: str):
    """Retrieve protein interactions for a given gene from IntAct."""
    response = requests.get(f"{INTACT_ENDPOINT}/search/{gene_id}")
    
    if response.status_code == 200:
        data = response.json()
        interactions = [
            {"interaction_id": item["id"], "participants": item["participants"]}
            for item in data.get("interactions", [])
        ]
        return interactions
    
    return []

def get_interactions(bridgedb_df):
    """Annotate genes with interaction data from IntAct."""
    api_available = check_endpoint_intact()
    if not api_available:
        warnings.warn("IntAct API endpoint is unavailable. Cannot retrieve data.", stacklevel=2)
        return pd.DataFrame(), {}

    # Record the start time
    start_time = datetime.datetime.now()

    # Filter for genes where target.source is "Ensembl"
    data_df = bridgedb_df[bridgedb_df["target.source"] == INTACT_INPUT_GENE_ID].copy()
    data_df = data_df.reset_index(drop=True)
    gene_list = list(set(data_df["target"].tolist()))

    # Retrieve interactions for each gene
    data_df["IntAct_interactions"] = data_df["target"].apply(
        lambda x: get_intact_interactions(x)
    )

    # Record the end time
    end_time = datetime.datetime.now()

    # Metadata details
    intact_metadata = {
        "datasource": KEGG,
        "metadata": {"source_version": intact_version},
        "query": {
            "size": len(gene_list),
            "input_type": INTACT_GENE_INPUT_ID,
            "number_of_added_edges": num_new_edges,
            "time": time_elapsed,
            "date": current_date,
            "url": INTACT_ENDPOINT,
        },
    }

    return data_df, metadata
