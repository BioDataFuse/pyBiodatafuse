Data sources
============

Here you can find the different data resources that you can successfully query with the pyBiodatafuse package.

.. module:: pyBiodatafuse.annotators


Bgee
~~~~~~~~~~~~~~~

`Bgee <https://www.bgee.org/>` is a database for retrieval and comparison of gene expression patterns across multiple animal species.

.. autofunction:: pyBiodatafuse.annotators.bgee.get_gene_expression()
    :members:

DisGeNET
~~~~~~~~~~~~~~~

`DisGeNET <https://www.disgenet.org/home/>` is a discovery platform containing one of the largest publicly available collections of genes and variants associated to human diseases.

.. autofunction:: pyBiodatafuse.annotators.disgenet.get_gene_disease()
    :members:

MolMeDB
~~~~~~~~~~~~~~~

`MolMeDB <https://molmedb.upol.cz/detail/intro>` is an open chemistry database about interactions of molecules with membranes.

.. autofunction:: pyBiodatafuse.annotators.molmedb.get_gene_mol_inhibitor()
    :members:

.. autofunction:: pyBiodatafuse.annotators.molmedb.get_mol_gene_inhibitor()
    :members:

OpenTargets
~~~~~~~~~~~~~~~

`OpenTargets database <https://www.opentargets.org/>` uses human genetics and genomics data for systematic drug target identification and prioritisation

.. autofunction:: pyBiodatafuse.annotators.opentargets.get_gene_location()
    :members:

.. autofunction:: pyBiodatafuse.annotators.opentargets.get_gene_go_process()
    :members:

.. autofunction:: pyBiodatafuse.annotators.opentargets.get_gene_reactome_pathways()
    :members:

.. autofunction:: pyBiodatafuse.annotators.opentargets.get_gene_tractability()
    :members:

.. autofunction:: pyBiodatafuse.annotators.opentargets.get_gene_drug_interactions()
    :members:

.. autofunction:: pyBiodatafuse.annotators.opentargets.get_gene_disease_associations()
    :members:

StringDB
~~~~~~~~~~~~~~~

`StringDB <https://string-db.org/>` aims to collect, score and integrate all publicly available sources of protein-protein interaction information, and to complement these with computational predictions.

.. autofunction:: pyBiodatafuse.annotators.stringdb.get_ppi()
    :members:

Wikidata
~~~~~~~~~~~~~~~

`Wikidata <https://www.wikidata.org/>` acts as central storage for the structured data of its Wikimedia sister projects including Wikipedia, Wikivoyage, Wiktionary, Wikisource, and others.

.. autofunction:: pyBiodatafuse.annotators.wikidata.get_gene_literature()
    :members:

.. autofunction:: pyBiodatafuse.annotators.wikidata.get_gene_cellular_component()
    :members:

WikiPathways
~~~~~~~~~~~~~~~

`Wikipathways <https://www.wikipathways.org/>` is an open science platform for biological pathways contributed, updated, and used by the research community.

.. autofunction:: pyBiodatafuse.annotators.wikipathways.get_gene_wikipathway()
    :members:

