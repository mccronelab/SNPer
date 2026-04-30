FROM ubuntu:25.04
ARG TARGETARCH

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
    mafft \
    autoconf \
    automake \
    libtool libtool-bin \
    pkg-config \
    libsigsegv-dev \
    openjdk-21-jdk-headless

ENV JAVA_HOME=/usr/lib/jvm/java-21-openjdk-amd64
ENV PATH="$JAVA_HOME/bin:$PATH"

COPY setup_snper-env.sh setup_snper-env.sh

RUN chmod +x /setup_snper-env.sh \
    && ./setup_snper-env.sh

RUN git clone https://github.com/timflutre/trimmomatic.git && \
    cd trimmomatic && \
    make && \
    make check && \
    make INSTALL=/usr/local/ install

RUN git clone https://github.com/jsh58/NGmerge.git && \
    cd NGmerge && \
    make && \
    mv NGmerge /usr/local/bin/

FROM ubuntu:24.04
COPY --from=nextstrain/nextclade:latest /usr/bin/nextclade /usr/local/bin/nextclade

ENV PATH="root/bin/:/snper-env/bin:$PATH"
ENV TRIMMOMATIC_JAR="/usr/local/bin/trimmomatic.jar"