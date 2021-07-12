FROM rocker/r-ver:4.1.0

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
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

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    "cmake" \
    "libgit2-dev" \
    && rm -rf /var/lib/apt/lists/*

# Include nanodisco toolbox
RUN mkdir /home/nanodisco

COPY code /home/nanodisco/code
# Retrieve repository
# RUN git clone --depth 1 --branch dev_minimap2 https://github.com/fanglab/nanodisco
# RUN cp -r /nanodisco/code /home/nanodisco/

# Install remaining dependencies from sources (nanopolish, bwa, samtools, R packages, MEME, bedtools).
COPY postInstall /
RUN bash /postInstall
# RUN bash /nanodisco/postInstall

# Define working directory.
WORKDIR /home/nanodisco

# Set environment variables.
ENV HOME /home/nanodisco

# Create folders for analysis
RUN mkdir /home/nanodisco/analysis
RUN mkdir /home/nanodisco/dataset

# Display start-up information
CMD echo "Usage:  docker run -it --name <IMAGE_NAME> fanglab/nanodisco bash" && \
    echo "" && \
    echo "This will start a sub-shell within your current directory named /home/nanodisco." && \
    echo "You can then run analysis, including:" && \
    echo "   Typing and fine-mapping methylation motif." && \
    echo "   Bin metagenomic contigs using methylation signal." && \
    echo "" && \
    echo "For more information, please consult https://github.com/fanglab/nanodisco"
