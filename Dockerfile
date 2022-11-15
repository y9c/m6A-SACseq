# x86_64
FROM debian:bullseye-slim
RUN echo "This is a analysis pipelie (v2.0) for SAC-seq data" > /README.md

ENV PATH="/opt/pipeline/bin:/bin:/opt/micromamba/bin:$PATH"
ENV MAMBA_ROOT_PREFIX="/opt/micromamba"
ENV MAMBA_EXE="/bin/micromamba"

# install system dependencies
RUN DEBIAN_FRONTEND=noninteractive apt-get update && DEBIAN_FRONTEND=noninteractive TZ=Etc/UTC apt-get -y --no-install-recommends install tzdata apt-utils wget git bzip2 make cmake xsltproc gcc g++ pkg-config zlib1g-dev libxml2-dev python3 python3-distutils default-jre && apt-get clean && rm -rf /var/lib/apt/lists/*
# install package by micromamba
RUN wget -qO- https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj -C /tmp bin/micromamba && cp /tmp/bin/micromamba $MAMBA_EXE && $MAMBA_EXE shell init -s bash -p $MAMBA_ROOT_PREFIX && eval "$($MAMBA_EXE shell hook --shell=posix)" && \
  micromamba install -n base -y -c conda-forge -c bioconda libzlib falco bowtie2 star samtools bedtools fastp seqtk blast subread pyfaidx && micromamba clean -afy
# install package by python-pip
RUN pip install --no-cache-dir snakemake==7.16 cutadapt==4.1 multiqc panoptes-ui numpy pandas pysam 
# install software directly
RUN wget -q https://github.com/getzlab/rnaseqc/releases/download/v2.4.2/rnaseqc.v2.4.2.linux.gz -O /tmp/rnaseqc.gz && gunzip /tmp/rnaseqc.gz && chmod +x /tmp/rnaseqc && cp /tmp/rnaseqc /bin/rnaseqc && rm -rf /tmp/rnaseqc.gz
# install blast2bam
RUN git clone --quiet --depth 1 https://github.com/guyduche/Blast2Bam /tmp/blast2bam && make -C /tmp/blast2bam/ -j && cp /tmp/blast2bam/bin/blast2bam /bin/blast2bam && rm -rf /tmp/blast2bam
# install rnaseqmut
RUN git clone --quiet --depth 1 https://github.com/davidliwei/rnaseqmut.git /tmp/rnaseqmut && mkdir -p /tmp/rnaseqmut/src/bamtools/build && cmake -B /tmp/rnaseqmut/src/bamtools/build /tmp/rnaseqmut/src/bamtools && make -C /tmp/rnaseqmut/src/bamtools/build -j && make -C /tmp/rnaseqmut/src -j && cp /tmp/rnaseqmut/bin/rnaseqmut /bin/rnaseqmut && rm -rf /tmp/rnaseqmut
# install umicollapse
RUN git clone --quiet --depth 1 https://github.com/Daniel-Liu-c0deb0t/UMICollapse.git /tmp/UMICollapse &&  cp /tmp/UMICollapse/umicollapse.jar /bin/umicollapse.jar && rm -rf /tmp/UMICollapse && \
  mkdir -p /bin/lib && wget -q -P /bin/lib https://repo1.maven.org/maven2/com/github/samtools/htsjdk/2.19.0/htsjdk-2.19.0.jar && wget -q -P /bin/lib https://repo1.maven.org/maven2/org/xerial/snappy/snappy-java/1.1.7.3/snappy-java-1.1.7.3.jar
# install cpup
RUN git clone --quiet --depth 1 https://github.com/y9c/cpup.git /tmp/cpup && make -C /tmp/cpup/ -j && cp /tmp/cpup/cpup /bin/cpup && rm -rf /tmp/cpup
# clean up and reduce size
RUN apt-get purge -y wget git bzip2 make cmake xsltproc gcc g++ pkg-config && apt-get clean && rm -rf /var/lib/apt/lists/*

ADD ./bin /opt/pipeline/bin
COPY ./Snakefile /opt/pipeline/Snakefile
COPY ./config.yaml /opt/pipeline/config.yaml
COPY ./sacseq /opt/pipeline/bin/sacseq
