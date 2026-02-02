Core Utilities
==============

Here you can find the core utility functions for data loading, ID mapping, and annotation orchestration.

.. module:: pyBiodatafuse

Data Loaders
~~~~~~~~~~~~

.. autofunction:: pyBiodatafuse.data_loader.create_df_from_file

.. autofunction:: pyBiodatafuse.data_loader.create_df_from_text

.. autofunction:: pyBiodatafuse.data_loader.create_df_from_dea

.. autofunction:: pyBiodatafuse.data_loader.filter_dea

ID Mapping
~~~~~~~~~~

.. autofunction:: pyBiodatafuse.id_mapper.bridgedb_xref

.. autofunction:: pyBiodatafuse.id_mapper.pubchem_xref

Homologs
~~~~~~~~

.. autofunction:: pyBiodatafuse.human_homologs.get_homologs

Annotation Orchestrator
~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: pyBiodatafuse.id_annotator.run_gene_selected_sources

Utility Functions
~~~~~~~~~~~~~~~~~

.. autofunction:: pyBiodatafuse.utils.get_identifier_of_interest

.. autofunction:: pyBiodatafuse.utils.collapse_data_sources

.. autofunction:: pyBiodatafuse.utils.combine_sources

.. autofunction:: pyBiodatafuse.utils.combine_with_homologs
