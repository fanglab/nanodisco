======================
Analyzing your dataset
======================

To analyze your dataset, ``nanodisco`` container can be used solely as an processing environment. This approach facilitate the processing of local datasets because no binding is required. Although, please note that the examples datasets (e.g. used in :ref:`Tool showcase<tool_showcase>`) cannot be retrieved because the container is not writable.

.. code-block:: sh

   singularity build nd_env nanodisco.sif # Create a default container named nd_env (not a sandbox)
   singularity exec nd_env nanodisco # Show general help
   # singularity exec nd_env nanodisco <subtask> [parameters]

Alternatively, by binding a directory with nanopore sequencing datasets using ``-B ./path/to/my_datasets:/home/nanodisco/dataset:ro``, you can find the content from ``./path/to/my_datasets`` in ``/home/nanodisco/dataset`` for processing. Note that ``:ro`` forces the binding in read-only mode to secure your data but will prevent any writing in the directory.

.. code-block:: sh

   # Start an interactive shell to use nanodisco and bind your nanopore sequencing dataset (*.fast5 files) to /home/nanodisco/dataset
   singularity shell --pwd /home/nanodisco --no-home -e -s /bin/bash -B ./path/to/my_datasets:/home/nanodisco/dataset:ro -w my_analysis # For new analysis
   # source .bashrc # If Singularity version >=3.3 for a pretty prompt
   # Type `exit` to leave the container

However, the container being writable also means that any modification made to the image's files, including system files and files from bound directories with ``-B`` (without ``:ro``) will be permanent. Note that, to maintain an organized container one can bind another directory containing a genome sequence file (.fasta) by adding ``-B ./path/to/my_genome_directory:/home/nanodisco/reference`` or directly use ``cp $genome_path ./my_analysis/home/nanodisco/reference/``.
