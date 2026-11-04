Welcome to Nanomotif's documentation!
=====================================

**Nanomotif** finds DNA methylation motifs in prokaryotic genomes and metagenomes from Nanopore
sequencing data. It finds the motifs of every bin (or genome) de novo from a modkit pileup, and its
MTase-linker module links them to the methyltransferase genes that produce them.

Nanomotif is fast and runs on large metagenomes: a 12-bin mock community with a 4 GB pileup takes about
20 seconds, and a 490-bin fecal metagenome with a 330 GB pileup took 21 minutes.

Start with :doc:`installation` and :doc:`required_files`, then :doc:`usage`. :doc:`how_it_works`
explains how motifs are found.

.. note::

   Version 1.2.0 removed the bin contamination and inclusion commands (``detect_contamination``,
   ``include_contigs``). For bin refinement from methylation, see
   `EpiRefine <https://github.com/SorenHeidelbach/epirefine>`_.

.. toctree::
   :maxdepth: 2

   installation
   required_files
   usage
   output
   how_it_works
   examples
   known_issues
   troubleshooting
   changelog
