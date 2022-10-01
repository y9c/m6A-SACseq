#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2022 Ye Chang yech1990@gmail.com
# Distributed under terms of the GNU license.
#
# Created: 2022-09-30 18:45

import argparse

import numpy as np
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument(
    "-f", "--files", nargs="+", help="<Required> input files", required=True
)
parser.add_argument(
    "-n", "--names", nargs="+", help="<Required> input names", required=True
)
parser.add_argument(
    "-o", "--output", help="<Required> output file", required=True
)
args = parser.parse_args()

uniq_names = list(dict.fromkeys(args.names))


# chr     pos     ref_base        strand  HeLa-WT-input-rep1_depth        HeLa-WT-input-rep1_mut  HeLa-WT-input-rep1_ratio        HeLa-WT-input-rep2_depth        HeLa-WT-input-rep2_mut  HeLa-WT-input-rep2_ratio        HeLa-WT-treated-rep1_depth      HeLa-WT-treated-rep1_mut        HeLa-WT-treated-rep1_ratio      HeLa-WT-treated-rep2_depth      HeLa-WT-treated-rep2_mut  HeLa-WT-treated-rep2_ratio
dfs = []
for f, n in zip(args.files, args.names):
    df = pd.read_csv(
        f,
        sep="\t",
        names=["chr", "pos", "ref_base", "mut", "depth"],
        low_memory=False,
    ).assign(sample=n)
    dfs.append(df)

df = (
    pd.concat(dfs, axis=0)
    .assign(strand=lambda x: np.where(x.ref_base == "A", "+", "-"))
    .pivot(
        index=["chr", "pos", "strand"],
        columns=["sample"],
        values=["mut", "depth"],
    )
    .swaplevel(axis=1)
    .loc[:, uniq_names]
    .fillna(0)
    .astype(pd.Int64Dtype())
)


df.columns = ["_".join(c) for c in df.columns]

df.to_csv(
    args.output, sep="\t", compression={"method": "gzip", "compresslevel": 5}
)
