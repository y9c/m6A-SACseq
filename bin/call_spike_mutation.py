#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2021 Ye Chang yech1990@gmail.com
# Distributed under terms of the MIT license.
#
# Created: 2021-09-05 17:08

"""picked mapped read within region."""

import sys

import pysam

infile = sys.argv[1]

samfile = pysam.AlignmentFile(infile, "rb")

for x in samfile.fetch():
    s = x.query_sequence
    if s:
        r = samfile.get_reference_name(x.reference_id)
        m = ""
        for i, j in x.get_aligned_pairs():
            if j and 15 <= j <= 19:
                m += s[i] if i else ""
        if len(m) == 5 and "N" not in m:
            print(r.replace("probe_", ""), m[:2] + "-" + m[3:], m[2], sep="\t")
