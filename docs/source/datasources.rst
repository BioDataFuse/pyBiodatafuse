Data sources
============

Here you can find the different data resources that you can successfully query with the pyBiodatafuse package.

.. module:: pyBiodatafuse.annotators

Bgee
~~~~~~~~~~~~~~~

`Bgee <https://www.bgee.org/>` is a database for retrieval and comparison of gene expression patterns across multiple animal species.

.. autofunction:: pyBiodatafuse.annotators.bgee.get_gene_expression

DisGeNET
~~~~~~~~~~~~~~~

`DisGeNET <https://www.disgenet.org/home/>` is a discovery platform containing one of the largest publicly available collections of genes and variants associated to human diseases.

.. autofunction:: pyBiodatafuse.annotators.disgenet.get_gene_disease

MolMeDB
~~~~~~~~~~~~~~~

`MolMeDB <https://molmedb.upol.cz/detail/intro>` is an open chemistry database about interactions of molecules with membranes.

.. autofunction:: pyBiodatafuse.annotators.molmedb.get_gene_compound_inhibitor

.. autofunction:: pyBiodatafuse.annotators.molmedb.get_compound_gene_inhibitor

OpenTargets
~~~~~~~~~~~~~~~

`OpenTargets database <https://www.opentargets.org/>`_ uses human genetics and genomics data for systematic drug target identification and prioritisation.

.. autofunction:: pyBiodatafuse.annotators.opentargets.get_gene_go_process

.. autofunction:: pyBiodatafuse.annotators.opentargets.get_gene_reactome_pathways

.. autofunction:: pyBiodatafuse.annotators.opentargets.get_gene_compound_interactions

.. autofunction:: pyBiodatafuse.annotators.opentargets.get_disease_compound_interactions

StringDB
~~~~~~~~~~~~~~~

`StringDB <https://string-db.org/>` aims to collect, score and integrate all publicly available sources of protein-protein interaction information, and to complement these with computational predictions.

.. autofunction:: pyBiodatafuse.annotators.stringdb.get_ppi

Wikidata
~~~~~~~~~~~~~~~

`Wikidata <https://www.wikidata.org/>` acts as central storage for the structured data of its Wikimedia sister projects including Wikipedia, Wikivoyage, Wiktionary, Wikisource, and others.

.. autofunction:: pyBiodatafuse.annotators.wikidata.get_gene_cellular_component

WikiPathways
~~~~~~~~~~~~~~~

`Wikipathways <https://www.wikipathways.org/>` is an open science platform for biological pathways contributed, updated, and used by the research community.

.. autofunction:: pyBiodatafuse.annotators.wikipathways.get_gene_wikipathways

MINERVA
~~~~~~~~~~~~~~~

`MINERVA <https://minerva.pages.uni.lu/doc/>`_ is a standalone webserver for visual exploration, analysis and management of molecular networks encoded in following systems biology formats.

.. autofunction:: pyBiodatafuse.annotators.minerva.get_minerva_components
.. autofunction:: pyBiodatafuse.annotators.minerva.get_gene_pathways

AOP-Wiki
~~~~~~~~~~~~~~~

`AOP-Wiki <https://aopwiki.org/>`_ is a database for Adverse Outcome Pathways (AOPs), which describe mechanistic information about the linkage between a molecular initiating event and an adverse outcome.

.. autofunction:: pyBiodatafuse.annotators.aopwiki.get_aops_gene

.. autofunction:: pyBiodatafuse.annotators.aopwiki.get_aops_compound

.. autofunction:: pyBiodatafuse.annotators.aopwiki.get_aops

IntAct
~~~~~~~~~~~~~~~

`IntAct <https://www.ebi.ac.uk/intact/>`_ provides a freely available, open source database system and analysis tools for molecular interaction data.

.. autofunction:: pyBiodatafuse.annotators.intact.get_gene_interactions

.. autofunction:: pyBiodatafuse.annotators.intact.get_compound_interactions

KEGG
~~~~~~~~~~~~~~~

`KEGG <https://www.kegg.jp/>`_ is a database resource for understanding high-level functions and utilities of biological systems from genomic and molecular-level information.

.. autofunction:: pyBiodatafuse.annotators.kegg.get_pathways

.. autofunction:: pyBiodatafuse.annotators.kegg.get_compounds

PubChem
~~~~~~~~~~~~~~~

`PubChem <https://pubchem.ncbi.nlm.nih.gov/>`_ is a database of chemical molecules and their activities against biological assays.

.. autofunction:: pyBiodatafuse.annotators.pubchem.get_protein_compound_screened

CompoundWiki
~~~~~~~~~~~~~~~

`CompoundWiki <https://compoundwiki.org/>`_ is a crowd-sourced database that provides comprehensive information about chemical compounds.

.. autofunction:: pyBiodatafuse.annotators.compoundwiki.get_compound_annotations

gProfiler
~~~~~~~~~~~~~~~

`gProfiler <https://biit.cs.ut.ee/gprofiler/>`_ is a web server for functional enrichment analysis and conversions of gene lists.

.. autofunction:: pyBiodatafuse.annotators.gprofiler.get_gene_enrichment

TFLink
~~~~~~~~~~~~~~~

`TFLink <https://tflink.net/>`_ is a comprehensive database of transcription factor-target interactions.

.. autofunction:: pyBiodatafuse.annotators.tflink.get_tf_target

MitoCarta
~~~~~~~~~~~~~~~

`MitoCarta <https://www.broadinstitute.org/mitocarta>`_ is an inventory of mammalian mitochondrial genes.

.. autofunction:: pyBiodatafuse.annotators.mitocarta.get_gene_mito_pathways
