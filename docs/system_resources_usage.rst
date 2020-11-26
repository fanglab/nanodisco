======================
System resources usage
======================

Generating current differences
==============================

**nanodisco preprocess**

* This command performs sequencing reads extraction from the fast5, and then mapping fasta on the reference genome.
* A typical runtime is ~40 min with 20 threads for a 5.5 Mb genome sequenced at 140x (~3 GB memory used).

**nanodisco difference**

* A typical runtime is ~16 min with 4 threads (2 jobs with 2 threads each) for 10 chunks at sequencing depth of ~85x (~10 GB memory used).
* The computation of current differences can be optimized by using ``-nj`` to use maximum available computing resources.
* HPC infrastructure can also be leverage by starting jobs on multiple nodes using ``-f`` and ``-l`` to specify unique chunk ranges. 

Individual bacteria, 4.5 Mb (or a metagenomic bin)
==================================================

**nanodisco motif**

* A typical runtime for the detection of 4 motifs from a 4.5 Mb genome is ~20 min with 4 threads.
* Runtime generally depends on the genome size and the number of methylation motifs.

**nanodisco characterize**

* A typical runtime for the characterization of 4 motifs is ~7 min with 4 threads (~10 GB memory used).
* Runtime generally depends on the genome size, the number of input methylation motifs, and the number of occurrences.

Microbiome, average size of 94 Mb (for methylation binning)
===========================================================

**nanodisco profile (motif-guided)**

* A typical runtime for the profiling of 60 motifs in a 90 Mb metagenome is ~3 min with 20 threads (~140 GB of memory used).
* Runtime generally depends on the number of input motifs and the metagenome size.

**nanodisco profile (automated)**

* A typical runtime for the profiling of all common motifs in a 90 Mb metagenome is ~41 hours with 20 threads (>240 GB of memory used).
* Runtime generally depends on the metagenome size (contigs longer than 100 kb by default).

**nanodisco select_feature**

* A typical runtime is ~45 min with 4 threads (~140 GB of memory used).
* Runtime generally depends on the diversity of methylation motifs in the microbiome and the metagenome size.

**nanodisco filter_profile**

* A typical runtime for the profiling of 3769 methylation features in a 90 Mb metagenome is ~25 min with 20 threads (>240 GB of memory used).
* Runtime generally depends on the diversity of methylation motifs in the microbiome and the metagenome size.

.. admonition:: Note

All tests were performed on an Intel(R) Xeon(R) Platinum 8168 CPU @ 2.70GHz.
