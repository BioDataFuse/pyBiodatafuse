# coding: utf-8

"""Python constant file."""

import os

HERE = os.path.dirname(os.path.realpath(__file__))  # Src path
DATA_DIR = os.path.join(HERE, "..", "data")  # Data path
RESOURCES_DIR = os.path.join(os.path.dirname(HERE), "..", "resources")  # Resources path

BRIDGEDB_DIR = os.path.join(DATA_DIR, "bridgedb")  # BridgeDb path
DISGENET_DIR = os.path.join(DATA_DIR, "disgenet")  # DisGeNET path
COMBINED_DIR = os.path.join(DATA_DIR, "combined_annotated_data")  # Cobimbined data path
