#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2021 Ye Chang yech1990@gmail.com
# Distributed under terms of the MIT license.
#
# Created: 2021-11-26 00:58

import sys

import numpy as np
import pandas as pd

df = pd.read_csv(sys.argv[1], sep="\t", low_memory=False)

col_names = list(df.columns)[:4]
sample_names = list(df.columns)[4:]

df["ref"] = df["ref_base"].str.upper()


for c in sample_names:
    df_acgt = df[c].str.split(",", expand=True).iloc[:, 1:5].astype(int)
    s_ref = np.where(
        df["ref"] == "A",
        df_acgt[1],
        np.where(
            df["ref"] == "C",
            df_acgt[2],
            np.where(
                df["ref"] == "G",
                df_acgt[3],
                np.where(df["ref"] == "T", df_acgt[4], 0),
            ),
        ),
    ).astype(int)
    s_depth = df_acgt.sum(axis=1)
    s_mut = s_depth - s_ref
    s_ratio = s_mut / s_depth
    df[c + "_depth"] = s_depth
    df[c + "_mut"] = s_mut
    df[c + "_ratio"] = s_ratio
    col_names += [c + "_depth", c + "_mut", c + "_ratio"]

df.loc[:, col_names].dropna().to_csv(sys.argv[2], sep="\t", index=False)
