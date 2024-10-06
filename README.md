[![Docker](https://img.shields.io/docker/pulls/y9ch/sacseq.svg)](https://hub.docker.com/r/y9ch/sacseq)
[![DOI:10.1038/s41587-022-01243-z](https://zenodo.org/badge/DOI/10.1038/s41596-022-00765-9.svg)](https://doi.org/10.1038/s41596-022-00765-9)
[![Citation Badge](https://api.juleskreuer.eu/citation-badge.php?doi=10.1038/s41596-022-00765-9)](https://www.nature.com/articles/s41596-022-00765-9/metrics)

# m<sup>6</sup>A-SAC-seq

## Overview of the workflow

<p align="center">
  <a href="https://y9c.github.io/m6A-SACseq/Introduction#gh-light-mode-only">
    <img src="./docs/scheme.svg" />
  </a>
  <a href="https://y9c.github.io/m6A-SACseq/Introduction#gh-dark-mode-only">
    <img src="./docs/scheme_dark.svg" />
  </a>
</p>

## How to use?

A [docker image](https://hub.docker.com/r/y9ch/sacseq) containing the source code and dependencies has been published for reproducibility. You can run it using the [singularity](https://sylabs.io/singularity) container runtime.

The entire analysis can be completed in just three steps:

1. **Specific the path (with label) of both rawdata and references for your project in a YAML format.**

   <details>
     <summary><code>data.yaml</code> for example<sup>(Click to expand)</sup></summary>

   ```yaml
   samples:
     HeLa-WT:
       input:
         rep1:
           - R1: ./rawdata/HeLa-WT-polyA-input-rep1-run1_R1.fq.gz
             R2: ./rawdata/HeLa-WT-polyA-input-rep1-run1_R2.fq.gz
           - R1: ./rawdata/HeLa-WT-polyA-input-rep1-run2_R1.fq.gz
             R2: ./rawdata/HeLa-WT-polyA-input-rep1-run2_R2.fq.gz
         rep2:
           - R1: ./rawdata/HeLa-WT-polyA-input-rep2-run1_R1.fq.gz
             R2: ./rawdata/HeLa-WT-polyA-input-rep2-run1_R2.fq.gz
           - R1: ./rawdata/HeLa-WT-polyA-input-rep2-run2_R1.fq.gz
             R2: ./rawdata/HeLa-WT-polyA-input-rep2-run2_R2.fq.gz
       treated:
         rep1:
           - R1: ./rawdata/HeLa-WT-polyA-treated-rep1-run1_R1.fq.gz
             R2: ./rawdata/HeLa-WT-polyA-treated-rep1-run1_R2.fq.gz
           - R1: ./rawdata/HeLa-WT-polyA-treated-rep1-run2_R1.fq.gz
             R2: ./rawdata/HeLa-WT-polyA-treated-rep1-run2_R2.fq.gz
         rep2:
           - R1: ./rawdata/HeLa-WT-polyA-treated-rep2-run1_R1.fq.gz
             R2: ./rawdata/HeLa-WT-polyA-treated-rep2-run1_R2.fq.gz
           - R1: ./rawdata/HeLa-WT-polyA-treated-rep2-run2_R1.fq.gz
             R2: ./rawdata/HeLa-WT-polyA-treated-rep2-run2_R2.fq.gz
   references:
     spike:
       fa: ./ref/spike_expand.fa
       bt2: ./ref/spike_expand
     spikeN:
       fa: ./ref/spike_degenerate.fa
       blast: ./ref/spike_degenerate
     rRNA:
       fa: ./ref/Homo_sapiens.GRCh38.rRNA.fa
       bt2: ./ref/Homo_sapiens.GRCh38.rRNA
     smallRNA:
       fa: ./ref/Homo_sapiens.GRCh38.smallRNA.fa
       bt2: ./ref/Homo_sapiens.GRCh38.smallRNA
     genome:
       fa: ./ref/Homo_sapiens.GRCh38.genome.fa
       star: ./ref/Homo_sapiens.GRCh38.genome
       gtf: ./ref/Homo_sapiens.GRCh38.genome.gtf
       gtf_collapse: ./ref/Homo_sapiens.GRCh38.genome.collapse.gtf
     contamination:
       fa: ./ref/contamination.fa
       bt2: ./ref/contamination
   ```

   _Read the [documentation](https://y9c.github.io/m6A-SACseq/Run-the-pipeline.html#refer-rawdata-and-references-in-the-configuration-file) on how to customize._

   </details>

2. **Run all the analysis by one command**:

   ```bash
   apptainer run docker://y9ch/sacseq:latest
   ```

   Note that when you storge your input file in a mounted partition, don't forget to add `--bind / -B` command to mount the partition.
   For example, using `apptainer run -B /data docker://sacseq:latest`...

    <details>
      <summary>default settings<sup>(Click to expand)</sup></summary>

   - default config file: `data.yaml`
   - default output dir: `./results`
   - default jobs in parallel: `48`

   _Read the [documentation](https://y9c.github.io/m6A-SACseq/Run-the-pipeline.html#customized-analysis-parameters) on how to customize._

   </details>

3. **View the analytics report and use the m<sup>6</sup>A sites for downstream analysis**.

   The output of all the steps will be in one folder (`./results`) under the current path. A webpage report of all the analysis will be in `./results/report.html` <sup>([example](https://y9c.github.io/m6A-SACseq/demo_output.html))</sup>.

## Documentation

https://y9c.github.io/m6A-SACseq/

## Citation

- Ge, R., Ye, C., Peng, Y. et al. m6A-SAC-seq for quantitative whole transcriptome m6A profiling. Nat Protoc (2022). https://doi.org/10.1038/s41596-022-00765-9

&nbsp;

<p align="center">
  <img
    src="https://raw.githubusercontent.com/y9c/y9c/master/resource/footer_line.svg?sanitize=true"
  />
</p>
<p align="center">
  Copyright &copy; 2021-present
  <a href="https://github.com/y9c" target="_blank">Chang Y</a>
</p>
<p align="center">
  <a href="https://github.com/y9c/m6A-SACseq/blob/master/LICENSE">
    <img src="https://img.shields.io/static/v1.svg?style=for-the-badge&label=License&message=GPLv3&logoColor=d9e0ee&colorA=282a36&colorB=c678dd" />
  </a>
</p>
