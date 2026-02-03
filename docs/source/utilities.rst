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

.. autofunction:: pyBiodatafuse.id_mapper.read_datasource_file

.. autofunction:: pyBiodatafuse.id_mapper.match_input_datasource

.. autofunction:: pyBiodatafuse.id_mapper.get_version_webservice_bridgedb

.. autofunction:: pyBiodatafuse.id_mapper.get_version_datasource_bridgedb

.. autofunction:: pyBiodatafuse.id_mapper.bridgedb_xref

.. autofunction:: pyBiodatafuse.id_mapper.check_smiles

.. autofunction:: pyBiodatafuse.id_mapper.get_cid_from_data

.. autofunction:: pyBiodatafuse.id_mapper.get_cid_from_pugrest

.. autofunction:: pyBiodatafuse.id_mapper.pubchem_xref

.. autofunction:: pyBiodatafuse.id_mapper.cid2chembl

Homologs
~~~~~~~~

.. autofunction:: pyBiodatafuse.human_homologs.check_endpoint_ensembl

.. autofunction:: pyBiodatafuse.human_homologs.check_version_ensembl

.. autofunction:: pyBiodatafuse.human_homologs.get_human_homologs

.. autofunction:: pyBiodatafuse.human_homologs.get_homologs

Annotation Orchestrator
~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: pyBiodatafuse.id_annotator.run_gene_selected_sources

Utility Functions
~~~~~~~~~~~~~~~~~

.. autofunction:: pyBiodatafuse.utils.get_identifier_of_interest

.. autofunction:: pyBiodatafuse.utils.collapse_data_sources

.. autofunction:: pyBiodatafuse.utils.create_or_append_to_metadata

.. autofunction:: pyBiodatafuse.utils.combine_sources

.. autofunction:: pyBiodatafuse.utils.combine_with_homologs

.. autofunction:: pyBiodatafuse.utils.check_columns_against_constants

.. autofunction:: pyBiodatafuse.utils.create_harmonized_input_file

.. autofunction:: pyBiodatafuse.utils.give_annotator_warning
