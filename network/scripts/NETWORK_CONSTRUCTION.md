##Network Construction

This document describes how `edgelist_final.txt` was generated from raw
experimental data. 

---

##Step 1 — Compute log2 fold-changes from raw batch data

Each batch CSV contains metabolite concentrations for WT and amino-acid
treatment conditions across biological replicates. Log2 fold-changes are
computed relative to the WT mean within each batch.

###Method
| Step | Detail |
|---|---|
| Reference | Mean WT concentration per metabolite per batch |
| Replicate FC | `log2(replicate_value / wt_mean)` per metabolite per replicate |
| Summary FC | `log2(mean(treatment) / mean(WT))` per metabolite per treatment |
| Significance | Two-sample t-test (treatment vs WT) per metabolite |
| Multiple testing | Benjamini–Hochberg FDR across all metabolite × treatment pairs |

---

##Step 2 — Build the metabolite interaction network

The network combines two edge types: direct perturbation edges (treatment →
metabolite) and metabolite–metabolite correlation edges derived from
replicate-level co-variation.

###Parameters

| Parameter | Value | Description |
|---|---|---|
| `QVALUE_THRESH` | 0.05 | FDR threshold for including a perturbation edge |
| `COR_THRESH` | 0.75 | Minimum absolute Pearson *r* to include a correlation edge |
| `ALPHA` | 0.6 | Weight fraction assigned to perturbation edges; `1 − ALPHA` assigned to correlation edges |

###Perturbation edges

Edges from each treatment to metabolites with a significant fold-change
(`qvalue < 0.05`). Edge weight:

```
weight = |log2FC| × (1 − qvalue)
```

Self-loops (treatment == metabolite) are removed.

###Correlation edges

Replicate-level log2FC values are assembled into a sample × metabolite
matrix, median-imputed within treatment groups, and scaled. Pairwise
Pearson correlations are computed across all metabolites. Edges with
`|r| > 0.75` are retained.

Correlation edges originating from perturbed nodes (i.e. nodes that already
appear as a source in the perturbation edge set) are removed to avoid
redundancy.

###Weight normalisation and combining

Both edge sets are min-max normalised to [0, 1] independently, then scaled before merging into a single edge list.

Each edge carries:
- `treatment` — source node
- `metabolite` — target node
- `norm_weight` — normalised weight
- `sign` — direction: +1 (positive FC or positive correlation), −1 (negative)

