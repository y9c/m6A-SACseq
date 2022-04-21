---
title: Introduction
nav_exclude: false
nav_order: 1
---

<!-- prettier-ignore-start -->
# Introduction
{: .fs-9 }
<!-- prettier-ignore-end -->

## Mechanism

![scheme_mechanism](scheme_mechanism.svg)

## Protocol

![scheme_protocol](scheme_protocol.svg)

1. 2-50 ng of poly-A enriched or ribosome RNA-depleted RNAs were fragmented and ligated, then divided in a 2:1 ratio.
1. 2/3 of the starting materials are labeled by MjDim1, while the remaining 1/3 serve as the untreated control.
1. After RT the cyclic allyl m<sup>6</sup>A sites (red dot) are converted to mismatches (orange dot), while unconverted m<sup>6</sup>A sites (cyan dot) in the control group are read as A (blue dot).

## Pipeline

![scheme_pipeline](scheme_pipeline.svg)

### Quality control:

We will first perform read trimming to remove adapters, primer sequences, molecular barcode (UMI), and low-quality bases using the `cutadapt` software.

### Data Processing:

All trimmed reads will be mapped to E. coli and Mycoplasma genome to filter biological contamination RNA using the `bowtie2` tool, and then unmapped reads will be mapped to spike-in sequence to filter RNA spike-in.
Similarly, unmapped reads will be mapped to the ribosomal RNA, small RNA and mRNA (whole genome) reference sequence sequentially.
After mapping, reads with identical UMI and mapped to the same location will be treated as PCR duplicates and dropped from downstream analysis.

### m6A Analyses:

m6A sites (mutation signal) were detected simultaneously, and mutation number and sequencing depth for each mutation position among all the samples were recorded.
Then putative m6A sites will be detected based on the mutation ratio, mutation number, and sequencing depth of each mutation position.
