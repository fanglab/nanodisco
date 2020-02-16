================
Data requirement
================

* Base called native nanopore sequencing dataset (\*.fast5)
* Base called WGA nanopore sequencing dataset (\*.fast5)
* (Optional) Reference genome (can be replace by *de novo* assembly with proper polishing)

Nanopore sequencing reads (\*.fast5) in “single” format (one read per file) are supported with Albacore and Guppy base calling but the ones in “multi” format (multiple reads per file) are only supported with Guppy base calling.

:ref:`Methylation typing and fine mapping<characterize_motifs>` is supported by a model we have trained using nanopore sequencing data base called with Albacore v2.3.4. Meanwhile, we are actively training new models for Guppy base caller and they will be released soon. We obtained good results with 75-100x coverage (see Supplementary Figure 8a in the preprint).

:ref:`Methylation binning of metagenomic contigs<microbiome>` can be performed on nanopore data base called with any software (e.g. Albacore, Guppy) but note that we expect better results with more accurate base calling, therefore the latest base caller version is recommended. Sequencing depth needed depends on the complexity of specific microbiome samples. In the samples described in the preprint, we obtained metagenomic bins with coverage from ~15x to ~370x.
