FROM ubuntu:26.04

# install necessary packages, tools. noninteractive stops apt from asking for
# geographic time zone (stops Docker)
RUN apt-get update && \
    DEBIAN_FRONTEND="noninteractive" apt-get -y install \
    samtools \
    ivar \
    picard-tools \
    bwa \
    fastqc \
    fastp \
    bedtools \
    git \
    python3 \
    python3-pip \
    python3-dev \
    python3-venv \
    build-essential \
    libparasail-dev \
    minimap2 \
    mafft \
    autoconf \
    automake \
    libtool libtool-bin \
    pkg-config \
    libsigsegv-dev \
    openjdk-21-jdk-headless

# PicardCommandLine is a wrapper script that resolves `java` from PATH, so no
# JAVA_HOME is needed.

COPY setup_snper-env.sh setup_snper-env.sh

RUN chmod +x /setup_snper-env.sh \
    && ./setup_snper-env.sh

RUN git clone https://github.com/jsh58/NGmerge.git && \
    cd NGmerge && \
    make && \
    mv NGmerge /usr/local/bin/

ENV PATH="/root/bin/:/snper-env/bin:$PATH"