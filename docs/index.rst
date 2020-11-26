.. _nanodisco-index:

nanodisco documentation
=======================

``nanodisco`` is a toolbox for *de novo* discovery of all the three types (6mA, 5mC and 4mC) of DNA methylation from individual bacteria and microbiomes using nanopore sequencing. For microbiomes, nanodisco also supports the use of DNA methylation patterns as natural epigenetic barcodes to facilitate high resolution metagenomic binning.

Features
--------

* *De novo* discover DNA methylation motifs, identify specific type (6mA, 5mC or 4mC, namely *typing*) of a methylation motif, and identify which specific position within the motif is methylated (namely *fine mapping*). 
* Perform metagenomic binning based on microbial DNA methylation pattern by constructing and clustering a methylation profile matrix. 
* Integrate the two functionalities above together for *de novo* methylation motif discovery from microbiomes, and metagenomic analysis.

Authors' notes
--------------

We are actively developing ``nanodisco`` to facilitate usage and broaden features. All feedback is more than welcome. You can reach us on twitter (`@iamfanggang <https://twitter.com/iamfanggang>`_ and `@AlanTourancheau <https://twitter.com/AlanTourancheau>`_) or directly through the `GitHub issues system <https://github.com/fanglab/nanodisco/issues>`_.

Installation
------------

``nanodisco`` is distributed as a fully functional image, bypassing the need to install any dependencies others than the virtualization software. We currently recommend using Singularity (v3.2.1 and above), which can be installed on Linux systems and is often the preferred solution by HPC administrators (`Quick Start <https://sylabs.io/guides/3.5/user-guide/quick_start.html>`_). ``nanodisco`` was tested extensively with Singularity v3.2.1 and v3.5.2.

.. code-block:: sh

   singularity pull --name nanodisco.sif shub://fanglab/nanodisco # Download the image from singularity-hub.org
   singularity build nd_env nanodisco.sif                         # Create a container named nd_env

Contents
--------

.. toctree::
   :maxdepth: 2

   overview
   tool_showcase
   detailed_tutorial
   commands_overview
   commands_details
   data_requirement
   analyze_dataset
   system_resources_usage
   faq
   citation

Contribute
----------

* `Download source code <https://github.com/fanglab/nanodisco>`__
* `Report issue <https://github.com/fanglab/nanodisco/issues>`__

Search
------

* :ref:`search`
