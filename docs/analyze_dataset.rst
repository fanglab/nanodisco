======================
Analyzing your dataset
======================

We highly recommend users to consult the :ref:`detailed tutorial <detailed-tutorial>` before analyzing your own dataset.

To analyze your dataset, ``nanodisco`` container can be used solely as a processing environment. This approach facilitates the processing of local datasets because no binding is required. Although, please note that the examples datasets used in :ref:`Tool showcase<tool-showcase>` and the :ref:`detailed tutorial <detailed-tutorial>` cannot be retrieved because the container is not writable with the following command.

.. code-block:: sh

   singularity build nd_env nanodisco.sif # Create a default container named nd_env (not a sandbox)
   singularity exec nd_env nanodisco # Show general help
   # singularity exec nd_env nanodisco <subtask> [parameters]

Alternatively, by binding a directory with nanopore sequencing datasets using ``-B ./path/to/my_datasets:/home/nanodisco/dataset:ro``, you can find the content from ``./path/to/my_datasets`` in ``/home/nanodisco/dataset`` for processing. Note that ``:ro`` forces the binding in read-only mode to secure your data but will prevent any writing in the directory.

.. code-block:: sh

   singularity build --sandbox nd_analysis nanodisco.sif # Create a writable container (directory) named nd_analysis
   # Start an interactive shell to use nanodisco and bind your nanopore sequencing dataset (*.fast5 files) to /home/nanodisco/dataset
   singularity run --no-home -B ./path/to/my_datasets:/home/nanodisco/dataset:ro -w nd_analysis
   # Type `exit` to leave the container

However, the container being writable also means that any modification made to the image's files, including system files and files from bound directories with ``-B`` (without ``:ro``) will be permanent. Note that, to maintain an organized container one can bind another directory containing a genome sequence file (.fasta) by adding ``-B ./path/to/my_genome_directory:/home/nanodisco/reference`` or directly use ``cp $genome_path ./my_analysis/home/nanodisco/reference/``.
