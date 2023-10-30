# coding: utf-8

"""Python constant file."""

import os

HERE = os.path.dirname(os.path.realpath(__file__))  # pyBiodatafuse path
DATA_DIR = os.path.join(HERE, "..", "..", "data")  # Data path
RESOURCES_DIR = os.path.join(HERE, "..", "..", "resources")  # Resources path

BRIDGEDB_DIR = os.path.join(DATA_DIR, "bridgedb")  # BridgeDb path
os.makedirs(BRIDGEDB_DIR, exist_ok=True)

DISGENET_DIR = os.path.join(DATA_DIR, "disgenet")  # DisGeNET path
os.makedirs(DISGENET_DIR, exist_ok=True)

COMBINED_DIR = os.path.join(DATA_DIR, "combined_annotated_data")  # Cobimbined data path
os.makedirs(COMBINED_DIR, exist_ok=True)
