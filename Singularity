Bootstrap: docker
From: rocker/r-ver:4.1.0

%help
For more information, please consult https://github.com/fanglab/nanodisco

# Install dependencies
%post
    # Install basic dependencies
    apt-get update && apt-get install -y --no-install-recommends \
      "vim=2:8.1.2269-1ubuntu5" \
      "git=1:2.25.1-1ubuntu3.2" \
      "wget=1.20.3-1ubuntu1" \
      "bzip2=1.0.8-2" \
      "parallel=20161222-1.1" \
      "python" \
      "ghostscript=9.50~dfsg-5ubuntu4.3" \
      "libcurl4-openssl-dev=7.68.0-1ubuntu2.7" \
      "libssl-dev=1.1.1f-1ubuntu2.8" \
      "libxml2-dev=2.9.10+dfsg-5ubuntu0.20.04.1" \
      "libxslt1-dev=1.1.34-4" \
      "zlib1g-dev=1:1.2.11.dfsg-2ubuntu1.2" \
      "libncurses5-dev=6.2-0ubuntu2" \
      "libncursesw5-dev=6.2-0ubuntu2" \
      "libexpat1-dev=2.2.9-1build1" \
      "libjson-perl=4.02000-2" \
      "libhtml-tree-perl=5.07-2" \
      "libbz2-dev=1.0.8-2" \
      "liblzma-dev=5.2.4-1ubuntu1" \
      "openmpi-bin=4.0.3-0ubuntu1" \
      "libopenmpi-dev=4.0.3-0ubuntu1" \
      "openssh-client=1:8.2p1-4ubuntu0.3" \
      "openssh-server=1:8.2p1-4ubuntu0.3" \
    && rm -rf /var/lib/apt/lists/*

    # Prepare for devtools dependencies
    apt-get update && apt-get install -y --no-install-recommends "cmake=3.16.3-1ubuntu1" "libgit2-dev=0.28.4+dfsg.1-2" && rm -rf /var/lib/apt/lists/*

    # Include nanodisco toolbox
    mkdir /home/nanodisco
    git clone --depth 1 --branch dev_latest https://github.com/fanglab/nanodisco
    cp -r /nanodisco/code /home/nanodisco/code

    # Install remaining dependencies from sources (nanopolish, minimap2, samtools, R packages, MEME, bedtools)
    # mv /nanodisco/postInstall /postInstall
    bash /nanodisco/postInstall

    # Define working directory
    cd /home/nanodisco

    # Create folders for analysis
    mkdir /home/nanodisco/analysis
    mkdir /home/nanodisco/dataset
    
    # Set default behavior
    cat > /.singularity.d/env/99-custom.sh <<EOF
export PS1="\[\033[36m\]\u\[\033[m\]@\[\033[32m\]nanodisco:\[\033[33;1m\]\w\[\033[m\]$ "
export HDF5_PLUGIN_PATH=/nanopolish-0.13.3/ont-vbz-hdf-plugin-1.0.1-Linux/usr/local/hdf5/lib/plugin
SINGULARITY_SHELL=/bin/bash
EOF

%environment
    export HOME=/home/nanodisco

%runscript
    cd /home/nanodisco
    exec /bin/bash
