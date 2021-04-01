.. _faq:

===
FAQ
===

.. _q_system:

* Q1: Is ``nanodisco`` available for other systems than Singularity in Linux?
     Sylabs also offers solutions to `install Singularity on Windows and Mac <https://sylabs.io/guides/3.5/admin-guide/installation.html#installation-on-windows-or-mac>`_. While we recommend using Singularity for running analysis with ``nanodisco``, we also provide a solution to use `Docker <https://www.docker.com/>`_ as image management software (:ref:`see instructions <using-docker>`).

.. _q_timezone:

* Q2: Why don't container times match the host machine?
     ``nanodisco`` image uses a default timezone. You can change it by exporting TZ variable as follow: 1. Use ``tzselect`` to find your time zone. 2. Enter ``export TZ='America/New_York'`` and add it to the ``/home/nanodisco/.bashrc`` file.

.. _q_ref_genome:

* Q3: What reference genome to use?
     You can use individual bacterial genomes or metagenomes *de novo* assembled from native nanopore sequencing data but the resulting assembly is better polished with either WGA nanopore reads or with Illumina/PacBio reads. This additional polishing is very helpful because systematic errors can exist at methylation sites with only native nanopore reads even after polishing the assembly with `nanopolish <https://github.com/jts/nanopolish>`_ or `medaka <https://github.com/nanoporetech/medaka>`_.

.. _q_wga:

* Q4: Do I need nanopore sequencing data from whole genome amplified (methylation free)? 
     Currently, yes, for two considerations. First, it is better to polish nanopore based *de novo* assembly of individual bacteria and metagenome with whole genome amplified (WGA) nanopore reads (or with Illumina/PacBio reads). This additional polishing is very helpful because systematic errors can exist at methylation sites with only native nanopore reads even after polishing the assembly with nanopolish or medaka. Second, although in theory, different modifications can be called directly from native nanopore sequencing reads, our observation, as described in the preprint, is that some methylation motifs have weak signal to noise ratio, and it is best to discover them by comparing native and WGA data.

.. _q_tsne_rep:

* Q5: Why can methylation binning from the same data result in slightly different figures?
     t-SNE dimensionality reduction uses a random seed in its computation. Although we have used a fixed seed before t-SNE dimensionality reduction to facilitate reproducibility, we found that the output may slightly differ depending on users’ platform/environment. This is a known behavior (`Rtsne issue <https://github.com/jkrijthe/Rtsne/issues/45>`_). Nevertheless, clustering results should be essentially the same.

.. _q_basecall_version:

* Q6: My nanopore data isn't base called with a supported software version, can I still use ``nanodisco``?
     Our methylation binning procedure is base caller agnostic as long as all nanopore datasets were base called with the same software version as well as the methylated motif detection. It's not the case for motif typing and fine mapping procedure because it relies on provided models that we trained with data from a specific base caller version. Substantial changes in base calling software would affect the methylation signal and could reduce model accuracy. Unfortunately, with the fast pace of base caller update, it's unlikely that we can keep up and create a model for every base caller version. However, most software update should have little effect on the signal and using the most recent model should provide good performance but you should proceed with cautions.

.. _q_methylation_event:

* Q7: Can I perform individual methylation event calling with ``nanodisco``?
     ``nanodisco`` is currently focused on *de novo* discovering, typing and fine mapping methylation motifs, the key challenge that prevents broad use of nanopore sequencing for the *de novo* study of bacterial epigenomes. ``nanodisco`` achieve very high accuracy for methylation motif typing and fine mapping. For a *de novo* discovered, typed and fine mapped methylation motif, other tools such as `Tombo <https://github.com/nanoporetech/tombo>`_ and `Nanopolish <https://github.com/jts/nanopolish>`_ may be used to detect individual methylation motif sites and estimate partial methylation. However, it is worth noting that, based on our analysis in the preprint, the signal to noise ratio and hence detectability at individual methylation motif sites can vary across different methylation motifs, e.g. it is relatively easier to determine 6mA methylation states at individual GATC sites and much harder at GAGG sites.

.. _q_basecall_version_req:

* Q8: Which software should I use to base call my data?
     Methylation typing and fine mapping is supported by a model we have trained using nanopore sequencing data basecalled with Albacore v2.3.4. In the meantime, we are actively training new models for Guppy base caller and they will be released soon. Methylation binning of metagenomic contigs can be performed on nanopore data basecalled with any software (e.g. Albacore, Guppy) but note that we expect better results with more accurate base calling, therefore the latest base caller version is recommended.

.. _q_flowcell:

* Q9: Which Nanopore instruments and flow cells are compatible?
     Datasets generated on the MinION and GridION with R9.4 flow cells are fully supported and all corresponding sequencing kits are compatible. For users with PromethION data (not tested yet), current tool is expected to work for methylation binning, and we are happy to help evaluate methylation typing and fine mapping. We also plan to support R10 data when it becomes more widely used.

.. _q_coverage:

* Q10: How much coverage is needed?
     We obtained good results for methylation typing and fine mapping with 75-100x coverage (see Supplementary Figure 8a in the preprint). For methylation binning of metagenomic contigs, sequencing depth needed depends on the complexity of specific microbiome samples. In the samples described in the preprint, we obtained metagenomic bins with coverage from ~15x to ~370x.

.. _q_rds:

* Q11: What are .RDS and .rds files?
     Those files are ``R``'s own data file format, which conserved all object properties. Use ``readRDS`` function to read a R data file. 

.. _q_eukaryote:

* Q12: Is it recommended to use ``nanodisco`` for methylation detection in eukaryotic species, such as human?
     ``nanodisco`` is currently focused on *de novo* discovering, typing and fine mapping methylation motifs, the key challenge that prevents broad use of nanopore sequencing for the *de novo* study of bacterial epigenomes. We do not recommend using ``nanodisco`` to look for methylation in human because 1) in the human genome, 5mC is mostly at CpG sites, for which multiple existing tools, such as `Nanopolish <https://github.com/jts/nanopolish>`_ and `Tombo <https://github.com/nanoporetech/tombo>`_, have been specifically trained; 2) for 6mA/4mC, if they exists in the human genome, it is at much lower level than 5mC, and would probably require different strategies for detection. Please see `Q7 <https://nanodisco.readthedocs.io/en/latest/faq.html#q-methylation-event>`_ in the FAQ page for more information. 

.. _q_hpc:

* Q13: How to best use an HPC infrastructure?
     For users with access to an HPC infrastructure, we recommend to use the job scheduler instead of relying on the built-in parallel approach. A loop can be used to generate all chunk start/end combinations (e.g. 1 to 5, then 6 to 10, etc.) for which a single job can be spawned by updating `-f` and `-l` (e.g. `-nj 1 -p 5 -nc 5 -f 1 -l 5`). This would makes resource management easier, and allows for a better use of available HPC computing power.
