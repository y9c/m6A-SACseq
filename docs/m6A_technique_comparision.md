---
title: m<sup>6</sup>A Technique Comparision
nav_exclude: false
nav_order: 4
---

<!-- prettier-ignore-start -->
# m<sup>6</sup>A Technique Comparision
{: .fs-9 }
<!-- prettier-ignore-end -->

| Method                                  | m<sup>6</sup>A-seq and MeRIP-seq            | m<sup>6</sup>A-SEAL                         | miCLIP-seq                                  | DART-seq                                    | m<sup>6</sup>A-REF-seq and MAZTER-seq       | m<sup>6</sup>A-label-seq                    | m<sup>6</sup>A-SAC-seq                      |
| --------------------------------------- | ------------------------------------------- | ------------------------------------------- | ------------------------------------------- | ------------------------------------------- | ------------------------------------------- | ------------------------------------------- | ------------------------------------------- |
| Site-specific                           | No                                          | No                                          | <span style= 'background:AFE1AF'>Yes</span> | <span style= 'background:AFE1AF'>Yes</span> | <span style= 'background:AFE1AF'>Yes</span> | <span style= 'background:AFE1AF'>Yes</span> | <span style= 'background:AFE1AF'>Yes</span> |
| Quantitative                            | No                                          | No                                          | No                                          | No[^a]                                      | <span style= 'background:AFE1AF'>Yes</span> | No                                          | <span style= 'background:AFE1AF'>Yes</span> |
| Covering all motifs                     | <span style= 'background:AFE1AF'>Yes</span> | <span style= 'background:AFE1AF'>Yes</span> | <span style= 'background:AFE1AF'>Yes</span> | No                                          | No                                          | <span style= 'background:AFE1AF'>Yes</span> | <span style= 'background:AFE1AF'>Yes</span> |
| Directly detecting m<sup>6</sup>A sites | No                                          | No                                          | No[^b]                                      | No[^b]                                      | No[^c]                                      | <span style= 'background:AFE1AF'>Yes</span> | <span style= 'background:AFE1AF'>Yes</span> |
| Frozen sample compatible                | <span style= 'background:AFE1AF'>Yes</span> | <span style= 'background:AFE1AF'>Yes</span> | <span style= 'background:AFE1AF'>Yes</span> | No[^d]                                      | <span style= 'background:AFE1AF'>Yes</span> | No[^e]                                      | <span style= 'background:AFE1AF'>Yes</span> |
| Antibody-free                           | No                                          | <span style= 'background:AFE1AF'>Yes</span> | No                                          | <span style= 'background:AFE1AF'>Yes</span> | <span style= 'background:AFE1AF'>Yes</span> | No                                          | <span style= 'background:AFE1AF'>Yes</span> |
| Starting material                       | 2-400 μg mRNA                               | 5 μg mRNA                                   | 20 μg mRNA                                  | 10 ng-1 μg total RNA                        | 100 ng mRNA                                 | 5 μg total RNA                              | 2-50 ng mRNA                                |
| Identified sites                        | 10,841-12,558                               | 8,605                                       | 21,494-37,183                               | 12,672                                      | 4,231-17,007                                | 2,479                                       | 71,547                                      |

[^a]: Can estimate the modification fraction based on C to U mutation ratio, which is indirect and inaccurate.
[^b]: Mutate flanking U (miCLIP-seq) or C (DART-seq) sites, thereby inferring the position of m<sup>6</sup>A sites.
[^c]: Cut at A sites, thereby inferring the fraction of m<sup>6</sup>A modification.
[^d]: Need to express the APOBEC1-YTH fusion protein in cells.
[^e]: Need to culture cells in Se-allyl-L-selenohomocysteine medium for 16 h.
