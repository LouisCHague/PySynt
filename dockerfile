# Base image
FROM ubuntu:20.04

# Prevents interative installations
ENV DEBIAN_FRONTEND=noninteractive

# Install Packages and Dependencies
RUN apt-get update \
    && apt-get install -y \
        build-essential \
        curl \
        wget \
    && rm -rf /var/lib/apt/lists/*

# Install MUMmer4 Dependencies
RUN apt-get update && \
    apt-get install -y \
    sudo \
    perl \
    sed \
    gawk \
    dash \
    fig2dev \
    gnuplot \
    xfig \
    g++ \
    && apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# MUMmer4 Install: Grab the tarball from GitHub
RUN mkdir /mummer && \
    curl -L https://github.com/mummer4/mummer/releases/download/v4.0.0rc1/mummer-4.0.0rc1.tar.gz -o /mummer/mummer-4.0.0rc1.tar.gz && \
    tar -xzvf /mummer/mummer-4.0.0rc1.tar.gz -C /mummer --strip-components=1

# Build, and install MUMmer
WORKDIR /mummer
RUN ./configure --prefix=/usr/local && \
    make && \
    make install

# https://github.com/mummer4/mummer/issues/48
RUN ldconfig

# Add MUMmer to PATH (Accessible everywhere) and remove mummer install directory
ENV PATH="/usr/local/bin:${PATH}" 
RUN rm -rf /mummer

# SAMTools Dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    libncurses5-dev \
    libncursesw5-dev \
    libbz2-dev \
    liblzma-dev \
    zlib1g-dev \
    libcurl4-gnutls-dev \
    libssl-dev \
    wget \
    nano \
    bzip2

# Install Samtools 1.15.1
RUN wget -q https://github.com/samtools/samtools/releases/download/1.15.1/samtools-1.15.1.tar.bz2 && \
    tar -xjf samtools-1.15.1.tar.bz2 && \
    cd samtools-1.15.1 && \
    ./configure --prefix=/usr/local && \
    make && \
    make install && \
    cd .. && \
    rm -rf samtools-1.15.1 samtools-1.15.1.tar.bz2

# Step down directory
WORKDIR /

# Set default command to open a shell
CMD ["/bin/bash"]