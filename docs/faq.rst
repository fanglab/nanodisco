.. _faq:

===
FAQ
===
 
* Q: Is ``nanodisco`` available for other systems than Singularity in Linux?
    Sylabs also offers solutions to `install Singularity on Windows and Mac <https://sylabs.io/guides/3.5/admin-guide/installation.html#installation-on-windows-or-mac>`_. While we recommend using Singularity for running analysis with ``nanodisco``, we also provide a solution to use `Docker <https://www.docker.com/>`_ as image management software (:ref:`see instructions <using-docker>`).

* Q: Why don't container times match the host machine?
     ``nanodisco`` image uses a default timezone. You can change it by exporting TZ variable as follow: 1. Use ``tzselect`` to find your time zone. 2. Enter ``export TZ='America/New_York'`` and add it to the ``/home/nanodisco/.bashrc`` file.

* Q: What reference genome to use?
     You can use individual bacterial genomes or metagenomes *de novo* assembled from native nanopore sequencing data but the resulting assembly is better polished with either WGA nanopore reads or with Illumina/PacBio reads. This additional polishing is very helpful because systematic errors can exist at methylation sites with only native nanopore reads even after polishing the assembly with `nanopolish <https://github.com/jts/nanopolish>`_ or `medaka <https://github.com/nanoporetech/medaka>`_.

* Q: Do I need nanopore sequencing data from whole genome amplified (methylation free)? 
     Currently, yes, for two considerations. First, it is better to polish nanopore based *de novo* assembly of individual bacteria and metagenome with whole genome amplified (WGA) nanopore reads (or with Illumina/PacBio reads). This additional polishing is very helpful because systematic errors can exist at methylation sites with only native nanopore reads even after polishing the assembly with nanopolish or medaka. Second, although in theory, different modifications can be called directly from native nanopore sequencing reads, our observation, as described in the preprint, is that some methylation motifs have weak signal to noise ratio, and it is best to discover them by comparing native and WGA data.

* Q: Why can methylation binning from the same data result in slightly different figures?
     t-SNE dimensionality reduction uses a random seed in its computation. Although we have used a fixed seed before t-SNE dimensionality reduction to facilitate reproducibility, we found that the output may slightly differ depending on users’ platform/environment. This is a known behavior (`Rtsne issue <https://github.com/jkrijthe/Rtsne/issues/45>`_). Nevertheless, clustering results should be essentially the same.

* Q: My nanopore data isn't base called with a supported software version, can I still use ``nanodisco``?
     Our methylation binning procedure is base caller agnostic as long as all nanopore datasets were base called with the same software version as well as the methylated motif detection. It's not the case for motif typing and fine mapping procedure because it relies on provided models that we trained with data from a specific base caller version. Substantial changes in base calling software would affect the methylation signal and could reduce model accuracy. Unfortunately, with the fast pace of base caller update, it's unlikely that we can keep up and create a model for every base caller version. However, most software update should have little effect on the signal and using the most recent model should provide good performance but you should proceed with cautions.

* Q: Can I perform individual methylation event calling with ``nanodisco``?
     ``nanodisco`` is currently focused on *de novo* discovering, typing and fine mapping methylation motifs, the key challenge that prevents broad use of nanopore sequencing for the *de novo* study of bacterial epigenomes. ``nanodisco`` achieve very high accuracy for methylation motif typing and fine mapping. For a *de novo* discovered, typed and fine mapped methylation motif, other tools such as Tombo may be used to detect individual methylation motif sites and estimate partial methylation. However, it is worth noting that, based on our analysis in the preprint, the signal to noise ratio and hence detectability at individual methylation motif sites can vary across different methylation motifs, e.g. it is relatively easier to determine 6mA methylation states at individual GATC sites and much harder at GAGG sites.

* Q: Which software should I use to base call my data?
     Methylation typing and fine mapping is supported by a model we have trained using nanopore sequencing data basecalled with Albacore v2.3.4. In the meantime, we are actively training new models for Guppy base caller and they will be released soon. Methylation binning of metagenomic contigs can be performed on nanopore data basecalled with any software (e.g. Albacore, Guppy) but note that we expect better results with more accurate base calling, therefore the latest base caller version is recommended.

* Q: How much coverage is needed?
     We obtained good results for methylation typing and fine mapping with 75-100x coverage (see Supplementary Figure 8a in the preprint). For methylation binning of metagenomic contigs, sequencing depth needed depends on the complexity of specific microbiome samples. In the samples described in the preprint, we obtained metagenomic bins with coverage from ~15x to ~370x.

* Q: What are .RDS and .rds files?
     Those files are ``R``'s own data file format, which conserved all object properties. Use ``readsRDS`` function to read a R data file. 
