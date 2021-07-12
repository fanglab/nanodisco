Bootstrap: docker
From: rocker/r-ver:3.5.3

%help
For more information, please consult https://github.com/fanglab/nanodisco

# Install dependencies
%post
    # Install basic dependencies
    apt-get update && apt-get install -y --no-install-recommends \
    "vim" \
    "git" \
    "wget" \
    "bzip2" \
    "parallel" \
    "python" \
    "ghostscript" \
    "libcurl4-openssl-dev" \
    "libssl-dev" \
    "libxml2-dev" \
    "libxslt1-dev" \
    "zlib1g-dev" \
    "libncurses5-dev" \
    "libncursesw5-dev" \
    "libexpat1-dev" \
    "libjson-perl" \
    "libhtml-tree-perl" \
    "libbz2-dev" \
    "liblzma-dev" \
    "openmpi-bin" \
    "libopenmpi-dev" \
    "openssh-client" \
    "openssh-server" \
    && rm -rf /var/lib/apt/lists/*

    # Prepare for devtools dependencies
    apt-get update && apt-get install -y --no-install-recommends "cmake" "libgit2-dev" && rm -rf /var/lib/apt/lists/*

    # Include nanodisco toolbox
    mkdir /home/nanodisco
    git clone --depth 1 --branch dev_minimap2 https://github.com/fanglab/nanodisco
    cp -r /nanodisco/code /home/nanodisco/code

    # Install remaining dependencies from sources (nanopolish, bwa, samtools, R packages, MEME, bedtools)
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
SINGULARITY_SHELL=/bin/bash
EOF

%environment
    export HOME=/home/nanodisco

%runscript
    cd /home/nanodisco
    exec /bin/bash
