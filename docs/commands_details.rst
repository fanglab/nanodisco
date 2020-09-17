.. _commands-details:

================
Commands details
================

``nanodisco <subtask> [options]``, ``<subtask>`` include:

* `preprocess`_: Extract reads (.fasta) from base called fast5 files and map reads on reference (meta)genome
* `chunk_info`_: Display chunks information regarding supplied reference (meta)genome
* `difference`_: Compute nanopore signal difference between a native and a whole genome amplified (WGA) dataset
* `merge`_: Combine nanopore signal difference for all processed chunks in directory
* `motif`_: *De novo* discovery of methylation motifs
* `refine`_: Generate refine plot for candidate methylation motifs
* `characterize`_: Predict the methylation type and fine-map the modification within de novo discovered methylation motifs file
* `coverage`_: Compute average coverage for each contig in a reference genome (uses bedtools genomecov)
* `profile`_: Compute the methylation profile matrix for a metagenome sample (methylation feature at common or expected methylation motifs)
* `select_feature`_: Select informative feature from a methylation profile matrix
* `filter_profile`_: Compute the methylation profile matrix for selected features for a metagenome sample
* `binning`_: Perform methylation binning, cluster metagenomic contigs according to methylation feature similarities using t-SNE
* `plot_binning`_: Plot results of methylation binning
* `version`_: Print version
* `help`_: Print help

.. _preprocess:

preprocess
==========

Extract reads (.fasta) from base called fast5 files and map reads on a reference (meta)genome.

**Usage:**

.. code-block:: none

   nanodisco preprocess -p <nb_threads> -f <path_fast5> -s <sample_name> -o <path_output> -r <path_reference_genome>
     -p : Number of threads to use.
     -f : Path to the directory containing .fast5 files (nanopore sequencing dataset). Note that fast5 files are searched recursively within the directory.
     -s : Name to the sample processed (e.g. Ecoli_native).
     -o : Path to output directory (e.g. ./analysis/Ecoli).
     -r : Path to a reference genome (i.e. fasta). Necessary index will be generated at runtime.
     --basecall_version : (Optional) Specify basecalling version if multiple ones available (<basecaller:version>, e.g. Guppy:3.2.4).
     -h : Print help.

**Output:**

* Fasta file (``<path_output>/<sample_name>.fasta``)
* Bam file (``<path_output>/<sample_name>.sorted.bam``)
* Bam index (``<path_output>/<sample_name>.sorted.bam.bai``)
* (Optional) Create reference index

.. _chunk_info:

chunk_info
==========

Display chunks information regarding the supplied reference (meta)genome.

**Usage:**

.. code-block:: none

   nanodisco chunk_info -r <path_reference_genome> [-t <target_region> -s <chunk_size>]
     -r : Path to a reference genome (i.e. fasta).
     -t : Specify a genomic region whose chunks need to be processed (e.g. chr1:2500-85000).
     -s : Size of chunk in bp to use when dividing the reference genome (default is 5000).
     -h : Print help.

**Output:**

* Default is the number of chunks in the reference genome
* With ``-t``, index of the first and last chunk to process for the targeted region

.. _difference:

difference
==========

Compute the current difference between a native and a WGA dataset.

**Usage:**

.. code-block:: none

   nanodisco difference -nj <nb_jobs> -nc <nb_chunks> -p <nb_threads> -i <path_input> -o <path_output> -w <name_WGA> -n <name_native> -r <path_genome> [-f <first_chunk> -l <last_chunk> + advanced parameters]
     -nj : Number of jobs to run in parallel (affect CPU and memory usage).
     -nc : Number of chunks to process in a row (affect memory usage).
     -p  : Number of threads to use.
     -i  : Path to input data folder containing .fasta, .bam, and .bam.bai generated with nanodisco preprocess.
     -o  : Path to output directory for current difference file and logs.
     -w  : Specify the name of WGA sample (same than for nanodisco preprocess -s <sample_name>).
     -n  : Specify the name of native sample (same than for nanodisco preprocess -s <sample_name>).
     -r  : Path to a reference genome (i.e. fasta).
     -h : Print help.

Advanced parameters. We recommend leaving them set to default values:

.. code-block:: none

     -f : First chunk to process. -l needs to be set. All chunks between -f and -l will be processed. All genome processed if not provided.
     -l : Last chunk to process. -f needs to be set. All chunks between -f and -l will be processed. All genome processed if not provided.
     -x : Execution type between seq or batch. Default is batch and seq is for development only.
     -a : IQR factor for outliers removal (0 to skip; smaller is harsher). Default is 1.5.
     -z : Type of additional signal normalization (0 is none, 1 is lm, and 2 is rlm). Default is 2.
     -b : Correct for strand bias (ori is no and revc is yes). Default is revc.
     -e : Minimum number of events per position. Default is 5.
     -j : Type of filtering for mapping. Default is noAddSupp.
     -k : Minimum mapped read length. Default is 0 (no filtering).
     --basecall_version : (Optional) Specify basecalling version if multiple ones available, for tracking only (<basecaller:version>, e.g. Guppy:3.2.4).

**Output:**

* Current difference files (``<path_output>/chunk.*.difference.rds``), one per chunk:

.. code-block:: none

   columns:
     contig       name of contig
     position     genomic position
     dir          genomic strand, fwd or rev
     strand       read strand, used when 2D nanopore reads
     N_wga        number of current values at this position and strand in WGA dataset
     N_nat        number of current values at this position and strand in native dataset
     mean_diff    current difference in pA
     t_test_pval  p-values from t-test
     u_test_pval  p-values from Mann-Whitney u-test

.. _merge:

merge
=====

Combine nanopore signal difference for all processed chunks in directory.

**Usage:**

.. code-block:: none

   nanodisco merge -d <path_difference> -o <path_output> -b <name_output>
     -d : Path to current differences directory (*.rds produced from nanodisco difference).
     -o : Path to output directory. Default is current directory.
     -b : Base name for outputting results (e.g. Ecoli_K12). Default is 'results'.
     -h : Print help.

**Output:**

* Current difference file (``<path_output>/<name_output>_difference.RDS``; same format as ``nanodisco difference`` output)

.. _motif:

motif
=====

*De novo* discovery of methylation motifs from current differences file.

**Usage:**

.. code-block:: none

   nanodisco motif -p <nb_threads> -b <name_output> -d <path_difference> -o <path_output> -r <path_genome> [+ advanced parameters]
     -p : Number of threads to use.
     -b : Base name for outputting results (e.g. Ecoli_K12). Default is 'results'.
     -d : Path to current differences file (*.RDS produced from nanodisco difference).
     -o : Path to output directory. Default is current directory.
     -r : Path to a reference genome (i.e. fasta).
     -h : Print help.

Advanced parameters. We recommend leaving them set to default values:

.. code-block:: none

     -c                : (Optional) Comma separated list of contigs (e.g. contig_1,contig_3).
     -m                : (Optional) Comma separated list of motifs (<motif_1,motif_2>, e.g. GATC,CCWGG).
     --contigs_file    : (Optional) Path to file with list of contigs (one per line).
     -a                : Disable manual motif discovery procedure (not recommended). Default is FALSE.
     -t                : Smoothed peaks p-values threshold for sequence selection (if double: peaks > <threshold> or if NA: top <nb_peaks> only). Default is NA.
     --nb_peaks        : Number of sequence with p-value peaks to keep for each round. Default is 2000.
     --stat_type       : Select which type of p-value sources used. Default is u_test_pval.
     --smooth_func     : Function to use for p-values smoothing. Default is sumlog.
     --smooth_win_size : Window size used for smoothing p-values. Default is 5.
     --peak_win_size   : Window size used for p-values peaks detection. Default is 2.

**Output:**

* Comma separated list of *de novo* discovered methylation motifs.
* Intermediate meme files (``<path_output>/motif_detection/``)
* Refinement plots for each motif without ``-a`` option.

.. _refine:

refine
======

Generate refine plot for candidate methylation motifs.

**Usage:**

.. code-block:: none

   nanodisco refine -p <nb_threads> -b <name_output> -d <path_difference> -o <path_output> -m <motif_1,motif_2> -M <motif_3,motif_4|all> -r <path_genome>
     -p : Number of threads to use.
     -b : Base name for outputting results (e.g. Ecoli_K12). Default is 'results'.
     -d : Path to current differences file (*.RDS produced from nanodisco difference).
     -o : Path to output directory. Default is current directory.
     -m : Comma separated list of discovered motifs (e.g. GATC,CCWGG).
     -M : Comma separated list of candidate motifs or 'all' to analyze '-m' motifs individually (e.g. GATC,CCWGG).
     -r : Path to a reference genome (i.e. fasta).
     -h : Print help.

**Output:**

* Refinement plots for the candidate motif(s) with ``-M <motif_3,motif_4>`` option, or each motif with ``-M all`` option.

.. _characterize:

characterize
============

Predict the methylation type and fine map the modification within *de novo* discovered methylation motifs file.

**Usage:**

.. code-block:: none

   nanodisco characterize -p <nb_threads> -b <name_output> -d <path_difference> -o <path_output> -m <motif1,motif2,...> -t <models> -r <path_genome>
     -p : Number of threads to use.
     -b : Base name for outputting results (e.g. Ecoli_K12). Default is 'results'.
     -d : Path to current differences file (*.RDS produced from nanodisco difference).
     -o : Path to output directory. Default is current directory.
     -m : Comma separated list of motifs following IUPAC nucleotide code (e.g. GATC,CCWGG).
     -t : Comma separated list of model to apply (nn: neural network, rf: random forest, or knn: k-nearest neighbor; e.g. nn,rf)
     -r : Path to a reference genome (i.e. fasta).
     -c : (Optional) Comma separated list of contigs (e.g. contig_1,contig_3).
     --contigs_file : (Optional) Path to file with list of contigs (one per line).
     -h : Print help.

**Output:**

* Identified methylation type and methylated position summarized in a heatmap (``Motifs_classification_Ecoli_<model_name>_model.pdf``) as presented in the preprint Figure 4d.
* Figure representing the data used to define the motif signature center as presented in the preprint Supplementary Figure 5a.

.. _coverage:

coverage
========

Compute average coverage for each contig in a reference genome (uses ``bedtools genomecov``).

**Usage:**

.. code-block:: none

   nanodisco coverage -b <path_mapping> -r <path_metagenome> -o <path_output>
     -b : Path of mapping data (.sorted.bam)
     -r : Path to a reference metagenome (i.e. fasta).
     -o : Path to output directory (.sorted.bam suffix replaced by .cov).

**Output:**

* Genomic coverage for each contig (``<path_output>/<bam_file_name>.cov``)

.. _profile:

profile
=======

Compute the methylation profile matrix for a metagenome sample (methylation feature at common or expected methylation motifs).

**Usage:**

.. code-block:: none

   nanodisco profile -p <nb_threads> -r <path_fasta> -d <path_difference> -w <path_WGA_cov> -n <path_NAT_cov> -b <analysis_name> -o <path_output> (-a || -m <motif1,motif2,...> || --motifs_file <path_motif>) [+ advanced parameters]
     -p : Number of threads to use.
     -r : Path to reference metagenome (.fasta).
     -d : Path to current differences file (*.RDS produced from nanodisco difference).
     -w : Path to WGA sample coverage (*.cov).
     -n : Path to native sample coverage (*.cov).
     -b : Base name for outputting results (e.g. Ecoli_K12). Default is 'results'.
     -o : Path to output directory. Default is current directory.
     -a : Compute methylation profile from predefined common motifs followed by filtering (automated binning; all|4mer|5mer|6mer|noBi). -a & -m & --motifs_file are exclusive.
     -m : Comma separated list of motifs following IUPAC nucleotide code (e.g. GATC,CCWGG). -a & -m & --motifs_file are exclusive.
     --motifs_file : Path to file with list of motifs (one per line) following IUPAC nucleotide code. -a & -m & --motifs_file are exclusive.

Advanced parameters. We recommend leaving them set to default values:

.. code-block:: none

     -c : Minimum coverage/number of current values needed at given position for methylation feature computation. Default is 10.
     --min_contig_len : Minimum length to consider a contig for feature selection. Default is 100000.

**Output:**

* Methylation profile matrix (``<path_output>/methylation_profile_<base_name>.RDS``):

.. code-block:: none

   columns:
     contig          name of contig
     motif           motif sequence (e.g. CCWGG)
     distance_motif  relative distance to first base of motif occurrence (0-based)
     signal_ratio    for development only. Expected strength of signal if motif known.
     dist_score      methylation feature value at relative distance (absolute average current difference across all motif occurrences)
     nb_occurrence   number of motif occurrence in the contig
   attribute:
     contig_coverage (data.frame):
         chr                 name of contig
         contig_length       contig length
         avg_cov.dataset_A   average contig coverage in -w dataset (WGA)
         avg_cov.dataset_B   average contig coverage in -n dataset (native)
         diff                coverage difference (A - B)
         ratio               coverage difference (A + 0.001)/(B + 0.001)

* With ``-a``\ , additional attribute ``min_contig_len`` for minimum length to consider a contig for feature selection.

.. _select_feature:

select_feature
==============

Select informative feature from a methylation profile matrix.

**Usage:**

.. code-block:: none

   nanodisco select_feature -p <nb_threads> -r <path_fasta> -s <path_profile> -b <analysis_name> -o <path_output> [+ advanced parameters]
     -p : Number of threads to use.
     -r : Path to reference metagenome (.fasta).
     -s : Path to methylation profile file (*.RDS produced from nanodisco profile).
     -b : Base name for outputting results (e.g. Ecoli_K12). Default is 'results'.
     -o : Path to output directory. Default is current directory.

Advanced parameters. We recommend leaving them set to default values:

.. code-block:: none

     --fsel_min_contig_len : Minimum length to consider a contig for feature selection. Default is 100000.
     --fsel_min_cov        : Minimum average coverage to consider a contig for feature selection. Default is 10.
     --fsel_min_motif_occ  : Minimum number of motif occurrences in a contig for feature selection. Default is 20.
     --fsel_min_signal     : Absolute threshold for considering a feature informative. Default is 1.5.

**Output:**

* Selected methylation features (``<path_output>/selected_features_<base_name>.RDS``):

.. code-block:: none

   columns:
     feature_name     feature identification (<motif>_<relative_distance>)
     motif            motif sequence (e.g. CCWGG)
     contigs_origin   list of contigs with significant feature (e.g. contig1|contig2)
   attribute:
     contig_coverage (data.frame), same format as in nanodisco profile

.. _filter_profile:

filter_profile
==============

Compute the methylation profile matrix for selected features for a metagenome sample.

**Usage:**

.. code-block:: none

   nanodisco filter_profile -p <nb_threads> -r <path_fasta> -d <path_difference> -f <path_feature> -b <analysis_name> -o <path_output> [+ advanced parameters]
     -p : Number of threads to use.
     -r : Path to reference metagenome (.fasta).
     -d : Path to current differences file (*.RDS produced from nanodisco difference).
     -f : Path to selected features file (*.RDS produced from nanodisco selected_feature).
     -b : Base name for outputting results (e.g. Ecoli_K12). Default is 'results'.
     -o : Path to output directory. Default is current directory.

Advanced parameters. We recommend leaving them set to default values:

.. code-block:: none

     -c : Minimum coverage/number of current values needed at given position for methylation feature computation. Default is 10.

**Output:**

* Filtered methylation profile matrix (``<path_output>/methylation_profile_<base_name>.RDS``):

.. code-block:: none

   columns:
     contig          name of contig
     motif           motif sequence (e.g. CCWGG)
     distance_motif  relative distance to first base of motif occurrence (0-based)
     signal_ratio    for development only. Expected strength of signal if motif known.
     dist_score      methylation feature value at relative distance (absolute average current difference across all motif occurrences)
     nb_occurrence   number of motif occurrence in the contig
   attribute:
     contig_coverage (data.frame), same format as in nanodisco profile

.. _binning:

binning
=======

Perform methylation binning, cluster metagenomic contigs according to methylation feature similarities using t-SNE.

**Usage:**

.. code-block:: none

   nanodisco binning -r <path_fasta> -s <path_profile> -b <analysis_name> -o <path_output> [+ advanced parameters]
     -r : Path to reference metagenome (.fasta).
     -s : Path to methylation profile file (*.RDS produced from nanodisco profile).
     -b : Base name for outputting results (e.g. Ecoli_K12). Default is 'results'.
     -o : Path to output directory. Default is current directory.

Advanced parameters. We recommend leaving them set to default values:

.. code-block:: none

     --min_motif_occ : Minimum number of motif occurrence to conserve entry in the methylation profile matrix. Default is 5.
     --min_contig_len : Minimum contig length to conserve entry in the methylation profile matrix. Default is 25000.
     --contig_weight_unit : Weight unit (bp) used for additional exaggeration in binning. Default is 50000.
     --max_relative_weight : Maximum relative weight a contig can have, weighting ceiling. Default is 0.05.
     --tsne_perplexity : t-SNE perplexity parameter. Default is 30.
     --tsne_max_iter : t-SNE maximum iteration parameter. Default is 2500.
     --tsne_seed : Seed set before t-SNE processing using set.seed function. Default is 101.
     --rdm_seed : Seed used for random number generation in missing value filling using set.seed function. Default is 42.

**Output:**

* Methylation binning results from t-SNE dimensionality reduction (``<path_output>/methylation_binning_<base_name>.RDS``):

.. code-block:: none

   columns:
     tSNE_1          x coordinate from t-SNE dimensionality reduction
     tSNE_2          y coordinate from t-SNE dimensionality reduction
     contig          name of contig
     contig_length   contig length
     id              contig identifier (e.g. species), is NA by default

.. _plot_binning:

plot_binning
============

Plot results of methylation binning.

**Usage:**

.. code-block:: none

   nanodisco plot_binning -r <path_fasta> -u <path_methylation_binning> -b <analysis_name> -o <path_output> [+ advanced parameters]
     -r : Path to reference metagenome (.fasta).
     -u : Path to methylation binning file (*.RDS produced from nanodisco binning).
     -b : Base name for outputting results (e.g. Ecoli_K12). Default is 'results'.
     -o : Path to output directory. Default is current directory.

Advanced parameters. We recommend leaving them set to default values:

.. code-block:: none

     -a                : Path to contig annotation. We expect two columns .txt or .RDS file with contig_name<tab>custom_name.
     -c                : Comma separated list of MGE contigs (e.g. contig_1,contig_3).
     --list_MGE_contig : Comma separated list of MGE contigs (e.g. contig_1,contig_3).
     --MGEs_file       : Path to file with list of MGE contigs (one per line).
     --xlim            : Optional x-axis zooming (e.g. -5:10).
     --ylim            : Optional y-axis zooming (e.g. -10:9).
     --min_contig_len  : Minimum length for plotting contigs. Default is 25000 bp.

**Output:**

* Methylation binning figure (``Contigs_methylation_tsne_<base_name>.pdf``) similar to Figure 5a-b in the preprint

.. _version:

version
=======

Print version.

**Usage:**

.. code-block:: none

   nanodisco version

.. _help:

help
====

Print help.

**Usage:**

.. code-block:: none

   nanodisco help
