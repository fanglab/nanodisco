=================
Commands overview
=================

``nanodisco <subtask> [options]``, ``<subtask>`` include:

* :ref:`preprocess`: Extract reads (.fasta) from base called fast5 files and map reads on reference (meta)genome
* :ref:`chunk_info`: Display chunks information regarding supplied reference (meta)genome
* :ref:`difference`: Compute nanopore signal difference between a native and a WGA dataset
* :ref:`merge`: Combine nanopore signal difference for all processed chunks in directory
* :ref:`motif`: *De novo* discovery of methylation motifs
* :ref:`characterize`: Predict the methylation type and fine-map the modification within de novo discovered methylation motifs file
* :ref:`coverage`: Compute average coverage for each contig in a reference genome (uses bedtools genomecov)
* :ref:`profile`: Compute the methylation profile matrix for a metagenome sample (methylation feature at common or expected methylation motifs)
* :ref:`select_feature`: Select informative feature from a methylation profile matrix
* :ref:`filter_profile`: Compute the methylation profile matrix for selected features for a metagenome sample
* :ref:`binning`: Perform methylation binning, cluster metagenomic contigs according to methylation feature similarities using t-SNE
* :ref:`plot_binning`: Plot results of methylation binning
* :ref:`help`: Print help
