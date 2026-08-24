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
    curl \
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

# Pinned: the aligner determines every published MSA, GFF and REF_POS, so it must not
# move under the snapshots on a rebuild. dpkg over BuildKit's TARGETARCH so this works
# under plain `docker build` and Apptainer conversion too.
ARG NEXTCLADE_VERSION=3.23.0
RUN set -eux; \
    arch="$(dpkg --print-architecture)"; \
    case "$arch" in \
        amd64) nc_arch='x86_64-unknown-linux-gnu' ;; \
        arm64) nc_arch='aarch64-unknown-linux-gnu' ;; \
        *) echo "unsupported architecture for nextclade: $arch" >&2; exit 1 ;; \
    esac; \
    curl -fsSL -o /usr/local/bin/nextclade \
        "https://github.com/nextstrain/nextclade/releases/download/${NEXTCLADE_VERSION}/nextclade-${nc_arch}"; \
    chmod +x /usr/local/bin/nextclade; \
    nextclade --version

ENV PATH="/root/bin/:/snper-env/bin:$PATH"