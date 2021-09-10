FROM rocker/r-ver:4.1.0

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
      "vim=2:8.1.2269-1ubuntu5" \
      "git=1:2.25.1-1ubuntu3.1" \
      "wget=1.20.3-1ubuntu1" \
      "bzip2=1.0.8-2" \
      "parallel=20161222-1.1" \
      "python" \
      "ghostscript=9.50~dfsg-5ubuntu4.3" \
      "libcurl4-openssl-dev=7.68.0-1ubuntu2.6" \
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

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    "cmake=3.16.3-1ubuntu1" \
    "libgit2-dev=0.28.4+dfsg.1-2" \
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
