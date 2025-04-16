FROM ubuntu:25.04

# install necessary packages, tools. noninteractive stops apt from asking for
# geographic time zone (stops Docker)
RUN apt-get update && \
    DEBIAN_FRONTEND="noninteractive" apt-get -y install \
    samtools \
    ivar \
    picard-tools \
    bwa \
    fastqc \
    bedtools \
    git \
    python3 \
    python3-pip \
    python3-dev \
    python3-venv \
    build-essential \
    libparasail-dev \
    minimap2 \
    mafft

COPY setup_snper-env.sh setup_snper-env.sh

RUN chmod +x /setup_snper-env.sh \
    && ./setup_snper-env.sh

ENV PATH="/snper-env/bin:$PATH"