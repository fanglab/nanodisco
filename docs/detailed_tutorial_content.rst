.. sectnum::
  :suffix: .

.. _detailed-tutorial-content:

To showcase the toolbox applications and facilitate the further understanding of the methods, we provide more detailed examples for the analysis of two datasets presented in the preprint. These datasets can be downloaded with the following commands from within a nanodisco container: ``get_data_bacteria`` and ``get_data_microbiome``.

.. _generate_differences:

Generating current differences
==============================

**Goals:**

* `Preprocessing fast5 datasets`_
* `Computing current differences`_ between native and whole genome amplified (WGA) data

Preprocessing fast5 datasets
----------------------------

**Description:**
For *de novo* methylation motifs detection and characterization (typing and fine mapping), ``nanodisco`` only needs the directory containing fast5 reads files and a reference genome as initial input. The following commands directly extract indexed read sequences (needed for ``nanopolish`` event alignment), and map them to the reference genome (missing reference genome index will be generated at runtime).  

**Inputs:**

#. Nanopore sequencing base called datatset
#. Reference (meta)genome file (.fasta)

**Outputs:** Extracted and indexed fasta files, as well as read mapping results (sorted.bam and sorted.bam.bai files).

.. code-block:: sh

   get_data_bacteria # If not already retrieved
   nanodisco preprocess -p 4 -f dataset/EC_WGA -s EC_WGA -o analysis/preprocessed_subset -r reference/Ecoli_K12_MG1655_ATCC47076.fasta
   nanodisco preprocess -p 4 -f dataset/EC_NAT -s EC_NAT -o analysis/preprocessed_subset -r reference/Ecoli_K12_MG1655_ATCC47076.fasta

Parameters ``-p`` number of threads to use, ``-f`` path to fast5 directory, ``-o`` path to output directory, and ``-r`` path to reference genome. See parameters detailed and advanced parameters in :ref:`preprocess` section.

Genomic chunks information
--------------------------

**Description:**
For current differences computation, we choose a processing approach based on genomic chunks. This make it easy to process large dataset in parallel by dividing processing across available resources (e.g. computing nodes). For a given reference genome, one can use the following command to display the number of genomic chunks that need to be processed to compute current differences for the whole genome.

**Inputs:**

#. Reference (meta)genome file (.fasta)

**Outputs:** Number of genomic chunks in the supplied reference (meta)genome file. Alternatively, chunks indexes for query genomic region.

.. code-block:: sh

   nanodisco chunk_info -r reference/Ecoli_K12_MG1655_ATCC47076.fasta

Parameter ``-r`` path to reference genome.

A subset of the reference genome can also be specifically targeted using the ``-t <chr:start-end>`` parameter which will display the subset of chunks that need to be processed.

.. code-block:: sh

   nanodisco chunk_info -r reference/Ecoli_K12_MG1655_ATCC47076.fasta -t CP014225.1:25050-84950

See parameters detailed and advanced parameters in :ref:`chunk_info` section.

Computing current differences
-----------------------------

**Description:**
Current differences are computed by first mapping nanopore sequencing events (plateaus in measured current) to genomic position using ``nanopolish eventalign`` (`nanopolish GitHub <https://github.com/jts/nanopolish>`_). Events current values are then aggregated for the native and the WGA datasets independently after removing of outliers. Current differences between the native and WGA datasets are then computed as well as additional statistics (e.g. p-values from Mann–Whitney U test).

**Inputs:**

#. Preprocessed native and WGA datasets directory (see `Preprocessing fast5 datasets`_)
#. Reference genome file (.fasta)

**Outputs:** Current differences files, one per chunk.

.. code-block:: sh

   nanodisco difference -nj 2 -nc 1 -p 2 -f 281 -l 290 -i analysis/preprocessed_subset -o analysis/difference_subset -w EC_WGA -n EC_NAT -r reference/Ecoli_K12_MG1655_ATCC47076.fasta

Parameters ``-nj`` number of jobs in parallel, ``-nc`` number of chunks to process in a row, ``-p`` number of threads per jobs, ``-f`` first chunk and ``-l`` last to process, ``-i`` path to input directory (used as output in ``nanodisco preprocess``), ``-o`` path to output directory, ``-w`` name of WGA sample and ``-n`` native sample (used in ``nanodisco preprocess``), and ``-r`` path to reference genome. See parameters detailed and advanced parameters in :ref:`difference` section. See how genomic chunks are defined in `Genomic chunks information`_ section.

Output file description (chunk.*.difference.rds):

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

Merge current differences files
-------------------------------

**Description:**
Individual current differences files can be combined using the following command.

**Inputs:**

#. Path to directory with a set of current differences files

**Outputs:** Combined current differences file.

.. code-block:: sh

   nanodisco merge -d analysis/difference_subset -o analysis -b EC_subset

Parameters ``-d`` path to input directory with individual current difference files (output directory used in ``nanodisco difference``), ``-o`` path to output directory, and ``-b`` name of the combined current differences file. Note that all chunk.*.difference.rds files will be combined in numerical order. Eventual gaps in chunks' coverage of the genome are not reported. See parameters detailed in :ref:`merge` section.

.. _bacteria:

Individual bacteria (or metagenomic bins)
=========================================

**Goals:**

* :ref:`De novo discovery of methylation motifs<de_novo_motif_detection>`
* :ref:`Characterize methylation motifs: typing and fine mapping<characterize_motifs>`

Examples dataset can be retrieve by executing ``get_data_bacteria`` within the container. This include a subset of fast5 reads from a native and WGA *E. coli* sample, *E. coli* reference genome, and current difference file for the complete *E. coli* genome.

.. _de_novo_motif_detection:

*De novo* discovery of methylation motifs
-----------------------------------------

**Description:**
When current differences are computed (see :ref:`Generating current differences<generate_differences>`), we also compute p-values for each genomic position (using Mann–Whitney U test by default). These p-values are combined locally with a sliding window-based approach using sumlog followed by peak detection. Flanking sequences around the center of peaks are then used as input for MEME motif discovery analysis. See Methods section in the preprint.

**Inputs:**

#. Current differences file (see :ref:`Generating current differences<generate_differences>`)
#. Reference genome file (.fasta)

**Outputs:** A list of *de novo* discovered methylation motifs

.. code-block:: sh

   get_data_bacteria # If not already retrieved
   nanodisco motif -p 4 -b test_EC -d /home/nanodisco/dataset/EC_difference.RDS -o /home/nanodisco/analysis -r /home/nanodisco/reference/Ecoli_K12_MG1655_ATCC47076.fasta -a

Parameters ``-p`` number of threads, ``-b`` base sample name, ``-d`` path to current differences file, ``-o`` path to output directory, ``-a`` automated processing without user input (Optional, not recommended),  and ``-r`` path to reference genome. We highly recommend not using the ``-a`` option and refining the potential motif before proceeding. A subset of contigs can be processed using ``-c`` or ``--contigs_file``. See parameters detailed and advanced parameters in :ref:`motif` section. **Runtime is ~20 min with 4 threads**.

.. _characterize_motifs:

Characterize methylation motifs: typing and fine mapping
--------------------------------------------------------

**Description:**
Identify the specific type (6mA, 5mC or 4mC, namely typing) of a methylation motif, and identify specific position within the motif is methylated (namely fine mapping). Detailed method is described in the preprint.

**Inputs:**

#. Current differences file (see :ref:`Generating current differences<generate_differences>`)
#. Reference genome file (.fasta)
#. Methylation motifs for which one wants to perform typing and fine mapping (see `de_novo_motif_detection`_)

**Outputs:** For each queried methylation motif, nanodisco identifies the methylation type and the methylated position summarized in a heatmap (``Motifs_classification_Ecoli_nn_model.pdf``). See Figure 4d in the preprint as an example. Filling colors correspond to percentage of occurrences classified to a specific class ranging from blue (0%) to red (100%). Grey columns correspond to positions without prediction. Motif position corresponds to position relative to motif start 0-based, e.g. position 1 for GATC motif corresponds to GATC. We include models for three classifiers that we found to be more accurate (presented in the preprint), while Supplementary Figure 6-7 in the preprint compared more classifiers with a leave-one-out method. In addition, plots representing the data used to find the motif signature center are made for each motif (similar to Supplementary Figure 5a in the preprint).

.. code-block:: sh

   get_data_bacteria # If not already retrieved
   nanodisco characterize -p 4 -b Ecoli -d dataset/EC_difference.RDS -o analysis/Ecoli_motifs -m GATC,CCWGG,GCACNNNNNNGTT,AACNNNNNNGTGC -t nn,rf,knn -r reference/Ecoli_K12_MG1655_ATCC47076.fasta

Parameters ``-p`` number of threads, ``-b`` base sample name, ``-d`` path to current differences file, ``-o`` path to output directory, ``-m`` comma separated list of motifs, ``-t`` comma separated list of models (nn/rf/knn), and ``-r`` path to reference genome. See parameters detailed and advanced parameters in :ref:`characterize` section. In this example, the current differences file (``EC_difference.RDS``) was generated on a whole E. coli nanopore sequencing dataset, from the preprint, using ``nanodisco difference`` (see :ref:`Generating current differences<generate_differences>`). **Runtime is ~7 min with 4 threads** (~10GiB memory used).

.. _microbiome:

Microbiome (methylation binning)
================================

Examples dataset can be retrieve by executing ``get_data_microbiome`` within the container. This includes current differences computed for a microbiome sample presented in our preprint (MGM1, see :ref:`Generating current differences<generate_differences>`), as well as the *de novo* metagenome assembly, and the intermediate methylation binning files (methylation profile from automated binning and motif driven binning).

**Goals:**

* `Compute metagenomic contigs coverage`_
* `Compute methylation profiles`_ (two approaches as described in the preprint):

  * :ref:`Motif driven <motif-driven>` (from *de novo* discovered motifs)
  * :ref:`Automated <automated>` (including methylation features selection)

* `Perform methylation binning`_

Compute metagenomic contigs coverage
------------------------------------

**Description:**
The metagenomic contigs binning from microbiome sample rely on methylation information measured from comparison between a native and a WGA dataset (current differences). The accuracy of those measures depends in part on the depth of sequence on each contigs (see :ref:`FAQ<faq>`). We used ``bedtools genomecov`` to compute contigs coverage for each dataset (`bedtools GitHub <https://github.com/arq5x/bedtools2>`_). Note that we do not provide the mapping files for the microbiome example but the coverage files are directly provided with ``get_data_microbiome``.

**Inputs:**

#. Mapped reads (see `Preprocessing fast5 datasets`_)
#. Reference genome file (.fasta)

**Outputs:** Metagenomic contigs coverage (``<sample_name>.cov``).

.. code-block:: sh

   get_data_bacteria # If not already retrieved. We do not provide reads for the microbiome sample
   # First run example commands for "nanodisco preprocess", see above
   nanodisco coverage -b analysis/preprocessed_subset/EC_NAT.sorted.bam -r reference/Ecoli_K12_MG1655_ATCC47076.fasta -o analysis/preprocessed_subset
   nanodisco coverage -b analysis/preprocessed_subset/EC_WGA.sorted.bam -r reference/Ecoli_K12_MG1655_ATCC47076.fasta -o analysis/preprocessed_subset

Parameters ``-b`` path to aligned reads (sorted.bam generated with ``nanodisco preprocess``), ``-r`` path to reference genome, ``-o`` path to output directory (output files will be analysis/EC_NAT.cov in this example). See parameters detailed in :ref:`coverage` section.

Compute methylation profiles
----------------------------

To perform methylation binning of metagenome contigs, we first need to extract methylation signal from the dataset by computing current difference (see :ref:`difference` section). Then, for each metagenomic contig, we create a methylation profile, which corresponds to the collection of signature averages from motifs of interest. Meanwhile, the motif signature corresponds to the ensemble of current differences distributions near a given motif, which are constructed by aggregating current differences for all motif occurrences in a contig. This signify that the methylation profile correspond to current difference averaged at relative position from a given motif, which are, during the methylation binning, consider as methylation features. 

Methylation profiles can be constructed by two approaches:

#. :ref:`Motif driven <motif-driven>`: a methylation profile matrix is computed for a set of specific motifs of interest on all contigs and all resulting methylation features will be conserved. The set of motif of interest can be define from prior knowledge regarding the microbiome sample (e.g. identified species, identified MTases, etc.) or from *de novo* discovered methylation motifs obtained with ``nanodisco motif`` analysis (see :ref:`motif` section). Methylation motifs can be found from the analysis of individual contigs of interest, from the analysis of bins found with automated methylation binning (``nanodisco profile -a``, see below), or from the analysis of bins found with conventional binning (e.g. with coverage and composition).
#. :ref:`Automated <automated>`: a methylation profile matrix is computed for a large set of commonly methylated motifs (n=200,000+) on large contigs only. The resulting methylation features are then filtered to conserve only informative features (the ones with strong signal). Those features are then computed for the remaining contigs in the metagenome.

.. _motif-driven:

**1. Methylation binning based on methylation profile of specific motifs:**

**Inputs:**

#. Current differences file (see :ref:`Generating current differences<generate_differences>`)
#. Metagenomic *de novo* assembly (.fasta)
#. Metagenomic contigs coverage files (see `Compute metagenomic contigs coverage`_)
#. Specific methylation motifs (see :ref:`De novo motif detection<de_novo_motif_detection>`)

**Outputs:** Methylation profile matrix

**Description:** Compute methylation profile for methylation motifs of interest for all contigs with sufficient coverage (>10x in native and WGA dataset).

.. code-block:: sh

   get_data_microbiome # If not already retrieved
   nanodisco profile -p 4 -r reference/metagenome.fasta -d dataset/metagenome_subset_difference.RDS -w dataset/metagenome_WGA.cov -n dataset/metagenome_NAT.cov -b MGM1_motif -o analysis/binning --motifs_file dataset/list_de_novo_discovered_motifs.txt

Parameters ``-p`` number of threads, ``-r`` path to reference genome, ``-d`` path to current differences file, ``-w`` and ``-n`` path to WGA and native coverage files (generated with ``nanodisco coverage``), ``-b`` name of sample/analysis, ``-o`` path to output directory, and ``--motifs_file`` path to file with list of motifs following IUPAC nucleotide code (one per line). See parameters detailed and advanced parameters in :ref:`profile` section. **Runtime is ~1 min with 4 threads** and ~4 Gb of memory used.

You can generate a similar methylation profile than in the preprint by applying the following command. Note that the operation need more memory therefore we also provided its output file (``dataset/methylation_profile_MGM1_motif.RDS``) within the example dataset retrieved with ``get_data_microbiome``. This allows you to skip this step and directly proceed to ``nanodisco binning``.

.. code-block:: sh

   get_data_microbiome # If not already retrieved
   nanodisco profile -p 20 -r reference/metagenome.fasta -d dataset/metagenome_difference.RDS -w dataset/metagenome_WGA.cov -n dataset/metagenome_NAT.cov -b MGM1_motif -o analysis/binning --motifs_file dataset/list_de_novo_discovered_motifs.txt

With this command the **runtime is ~3 min with 20 threads** and ~140 Gb of memory used.

.. _automated:

**2. Methylation binning based on methylation profile without specific de novo discovered motifs:**


For this section, all inputs files are already available within the example dataset retrieved with ``get_data_microbiome`` allowing you to start from any step.

**Description:** Compute methylation profile for predefined common methylation motifs (n=210,176) for a subset of long contigs (by default >100kbp) with sufficient coverage (>10x in native and WGA dataset).

**Inputs:**

#. Current differences file (see :ref:`Generating current differences<generate_differences>`)
#. Metagenomic *de novo* assembly (.fasta)
#. Metagenomic contigs coverage files (see `Compute metagenomic contigs coverage`_)

**Output:** A filtered methylation profile matrix

.. code-block:: sh

   get_data_microbiome # If not already retrieved
   nanodisco profile -p 4 -r reference/metagenome.fasta -d dataset/metagenome_subset_difference.RDS -w dataset/metagenome_WGA.cov -n dataset/metagenome_NAT.cov -b MGM1_auto_subset -o analysis/binning -a 4mer

Parameters ``-p`` number of threads, ``-r`` path to reference genome, ``-d`` path to current differences file, ``-w`` and ``-n`` path to WGA and native coverage files (generated with ``nanodisco coverage``), ``-b`` name of sample/analysis, ``-o`` path to output directory, and ``-a`` request automated methylation binning (4mer, 5mer, 6mer, noBi, or all). See parameters detailed and advanced parameters in :ref:`profile` section. This procedure **runtime is ~2 min with 4 threads**.

You can generate methylation profile with the same depth than in the preprint by applying the following command. Note that the operation is computationally intensive therefore we also provided its output file (``dataset/methylation_profile_MGM1_auto.RDS``) within the example dataset retrieved with ``get_data_microbiome``. This allows you to skip this step and directly proceed to ``nanodisco select_feature``.

.. code-block:: sh

   nanodisco profile -p 20 -r reference/metagenome.fasta -d dataset/metagenome_difference.RDS -w dataset/metagenome_WGA.cov -n dataset/metagenome_NAT.cov -b MGM1_auto -o analysis/binning -a all

This is a long procedure (>24 hours with 20+ threads depending on the metagenome characteristics) because many potential methylation motifs are considered before methylation features filtering. This procedure **runtime is ~41 hours with 20 threads** and more than 240 GB of memory used.

**Description:** Select informative methylation features from the previously generated methylation profile. Informative methylation features are the ones with more than 20 occurrences and an absolute feature values >=1.5 in at least one contig.

.. code-block:: sh

   get_data_microbiome # If not already retrieved
   nanodisco select_feature -p 4 -r reference/metagenome.fasta -s dataset/methylation_profile_MGM1_auto.RDS -b MGM1_auto -o analysis/binning

Parameters ``-p`` number of threads, ``-r`` path to reference genome, ``-d`` path to current differences file, ``-w`` and ``-n`` path to WGA and native coverage files (generated with ``nanodisco coverage``), ``-b`` name of sample/analysis, ``-o`` path to output directory, and ``-a`` request automated methylation binning. See parameters detailed and advanced parameters in :ref:`select_feature` section. **Runtime is ~45 min with 4 threads** and ~140 Gb of memory used.

**Description:** Compute methylation profile for select informative methylation features for all contigs with sufficient coverage (>10x in native and WGA dataset).

.. code-block:: sh

   get_data_microbiome # If not already retrieved
   nanodisco filter_profile -p 20 -r reference/metagenome.fasta -d dataset/metagenome_difference.RDS -f dataset/selected_features_MGM1_auto.RDS -b MGM1_auto_filtered -o analysis/binning

Parameters ``-p`` number of threads, ``-r`` path to reference genome, ``-d`` path to current differences file, ``-w`` and ``-n`` path to WGA and native coverage files (generated with ``nanodisco coverage``), ``-b`` name of sample/analysis, ``-o`` path to output directory, and ``-a`` request automated methylation binning. See parameters detailed and advanced parameters in :ref:`filter_profile` section. **Runtime is ~25 min with 20 threads** and more than 240 GB of memory used.

Perform methylation binning
---------------------------

**Description:** Compute methylation profile for methylation motifs of interest for all contigs with sufficient coverage (>10x in native and WGA dataset).
The methylation profile generated with ``nanodisco filter_profile`` or ``nanodisco profile`` (1. Automated or 2. Motif driven) can then be further processed to reveal bins of contigs with similar methylation profile using t-SNE dimensionality reduction.

**Inputs:**

#. Methylation profile matrix (see `Compute methylation profiles`_)
#. Metagenomic *de novo* assembly (.fasta)
#. (Optional) Annotation for metagenome contigs (e.g. species of origin)
#. (Optional) List of contigs from Mobile Genetic Elements (MGEs)

**Outputs:** t-SNE scatter plots that demonstrates the species level clustering of metagenomic contigs as presented in the preprint Figure 5a-b or Supplementary Figure 11-12.

.. code-block:: sh

   nanodisco binning -r reference/metagenome.fasta -s dataset/methylation_profile_MGM1_motif.RDS -b MGM1_motif -o analysis/binning
   # OR
   nanodisco binning -r reference/metagenome.fasta -s dataset/methylation_profile_MGM1_auto_filtered.RDS -b MGM1_auto -o analysis/binning

Parameters ``-r`` path to reference genome, ``-s`` path to methylation profile file (generated with ``nanodisco profile``), ``-b`` name of sample/analysis, and ``-o`` path to output directory.  See parameters detailed and advanced parameters in :ref:`binning` section.

Results of the methylation binning can then be plot using the following function:

.. code-block:: sh

   nanodisco plot_binning -r reference/metagenome.fasta -u analysis/binning/methylation_binning_MGM1_motif.RDS -b MGM1_motif -o analysis/binning -a reference/motif_binning_annotation.RDS --MGEs_file dataset/list_MGE_contigs.txt
   # OR
   nanodisco plot_binning -r reference/metagenome.fasta -u analysis/binning/methylation_binning_MGM1_auto.RDS -b MGM1_auto -o analysis/binning -a reference/motif_binning_annotation.RDS --MGEs_file dataset/list_MGE_contigs.txt

Parameters ``-r`` path to reference genome, ``-u`` path to methylation binning file (generated with ``nanodisco binning``), ``-b`` name of sample/analysis, ``-o`` path to output directory, ``-a`` (optional) path to contigs annotation, and ``--MGEs_file`` (optional) path to file with list of MGE contigs (one per line). See parameters detailed and advanced parameters in :ref:`plot_binning` section. Output file from ``nanodisco binning`` (``methylation_binning_<sample_name>.RDS``) can also be open and used within R in interactive mode if you want to create your own figures or extract specific cluster details. 

**Output files:** The above commands generate a motif driven methylation binning of mouse gut microbiome metagenome contigs as shown in Figure 5a of the preprint. Methylation status of common motifs (n=210,176) were screened across large contigs (>=100 kb) through computation of methylation feature vectors. Informative features were selected and their status evaluated across remaining contigs. Resulting methylation features are projected on two dimensions using t-SNE and bins were identified (see preprint Supplementary Figure 11a). For each bin, we performed :ref:`de novo motif detection<de_novo_motif_detection>` and generate a combined list of methylation motifs. Motif driven methylation binning was then performed for two rounds to further expose additional bins and methylation motifs (see preprint Supplementary Figure 11b and c). A final list of methylation motifs was then created (``list_de_novo_discovered_motifs.txt``) and used in this ultimate round of methylation binning. Contigs are colored based on bin identities with point sizes matching contig length according to legend.
