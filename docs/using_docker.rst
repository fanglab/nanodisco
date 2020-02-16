:orphan:

.. _using-docker:

Install Docker
==============

``nanodisco`` is distributed as a fully functional image bypassing the need to install any dependencies others than the virtualization software. We currently recommend using Singularity (v3.2.1 and above), which can be installed on Linux systems and is often the preferred solution by HPC administrators (`Quick Start <https://sylabs.io/guides/3.5/user-guide/quick_start.html>`_). ``nanodisco`` was tested extensively in Singularity v3.2.1 and v3.5.2.

For the analysis of individual bacteria on desktops and laptops, an alternative to Singularity is Docker, which can be install on multiple Linux platforms (\ `CentOS <https://docs.docker.com/install/linux/docker-ce/centos/>`_\ , `Debian <https://docs.docker.com/install/linux/docker-ce/debian/>`_\ , `Fedora <https://docs.docker.com/install/linux/docker-ce/fedora/>`_\ , `Ubuntu <https://docs.docker.com/install/linux/docker-ce/ubuntu/>`_\ ), as well as on `Windows <https://docs.docker.com/docker-for-windows/>`_\ , and on `Mac <https://docs.docker.com/docker-for-mac/>`_ following the documentation. For the analysis of microbiome data, we currently recommend using Singularity. Notes that allowed container resources (the number of CPUs and amount of memory) need to be configured for Windows and Mac from "Preferences" or "Setting" > "Advances". We recommend setting the memory limit as high as possible for your system. Note that processing one chunk with nbCPU=2 uses ~ 5.1GiB of memory, which will increase if using more than 2 CPUs (~ 5.9GiB for nbCPU=3) as well as processing more than one chunk in a row (~ 5.8GiB for two chunks).

Install nanodisco with Docker
=============================

.. code-block:: sh

   docker pull touraa01/smtm # Download the image from hub.docker.com

   # Create and start an interactive shell to use nanodisco
   docker run -it --name my_analysis touraa01/smtm bash # For running example analysis
   # Type `exit` to leave the container

The image retrieved from `Docker Hub <https://hub.docker.com/>`_ with ``docker pull`` is already built and can be reused at will. The command ``docker run`` creates and run a container named ``my_analysis`` from the image. This command will start an interactive shell within ``my_analysis`` container.

Using nanodisco with Docker
===========================

By binding a directory with nanopore sequencing datasets using ``-v /full_path/to/my_datasets:/home/nanodisco/dataset:ro``\ , you can find the content from ``/full_path/to/my_datasets`` in ``/home/nanodisco/dataset`` for processing. Note that ``:ro`` forces the binding in read-only mode to secure your data but will prevent any writting in the directory.

.. code-block:: sh

   # Start an interactive shell to use nanodisco and bind your nanopore sequencing dataset (*.fast5 files) to /home/nanodisco/dataset
   docker run -it --name my_analysis -v /full_path/to/my_datasets:/home/nanodisco/dataset:ro touraa01/smtm bash # For new analysis
   # Type `exit` to leave the container

With docker, you cannot directly reuse the same container and bind another dataset (\ ``-v``\ ) but you can create multiple containers to compartmentalize analysis of different datasets (e.g. my_analysis and my_analysis2). Containers build with those instructions are writable but not accessible with a file explorer or terminal. You can retrieve files or directories from a running or stopped container using ``docker cp my_analysis:$source_path $destination_path``. Alternatively you can directly create the container with another bind mount directory (\ ``-v``\ ) where analysis results can be saved and shared between the container and the host meaning that results from nanodisco analysis can be access with a file explorer or terminal. Note that, to maintain an organized container you can bind another directory containing a genome sequence file (.fasta) by adding ``-v /full_path/to/my_genomes:/home/nanodisco/reference`` or directly use ``docker cp $genome_path my_analysis:/home/nanodisco/reference/``. However, the container being writable also means that any modification made to the image's files, including system files and files from binded directories with ``-v`` (without ``:ro``\ ) will be permanent.

Useful Docker commands
======================

.. code-block:: sh

   # Reusing a container
   docker start my_analysis
   docker exec -it my_analysis bash

   # Stopping a container
   # Type `exit` to leave the container
   docker stop my_analysis

   docker container ls -a # Show all existing containers (without -a only show running ones)
