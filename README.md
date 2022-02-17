# m6A-SAC-seq

## Overview of the workflow

## How to use?

For reproducibility, source code and dependencies have been packaged into one [docker image](https://hub.docker.com/r/y9ch/sacseq).

You can run it thought [singularity](https://sylabs.io/singularity) container runtime.

- Specific the path and the label of your data a yaml file:

for example,

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
```

- Then run all the analysis by one command:

```bash
singularity exec docker://y9ch/sacseq:latest sacseq data.yaml
```

## Documentation

https://github.com/y9c/m6A-SACseq/wiki

## Citation

- (nbt)
