# Delineating Mutation Bias and Selection during Plant Development

This repository contains SLiM simulations and R code used to quantify genic hypomutation in plants and to separate **true mutation-rate reductions** from **apparent reductions caused by selection on nonsynonymous sites**.

## Contents

- `code/`
  - `functions.R` – shared helpers (loading TAIR10, feature overlap, plotting, utilities).
  - `du_estimate.R` (or inline in analysis scripts) – **core function** that takes observed genic and reference mutation counts, the CDS fraction, and a neutral NS/S ratio, and returns:
    - observed genic/reference multiplier (`duO = µ_G / µ_ref`);
    - selection-corrected multiplier (`duEst`);
    - the nonsynonymous correction term (`RNS`);
    - class-specific multipliers for NS, S, and noncoding genic mutations.
  - plotting helpers to reproduce the figures (genic vs intergenic, Δµ vs M, NS/S diagnostics).

- `data/`
  - `slim_out/*.csv` – mutation tables exported by SLiM for each parameter combination.
  - annotated Arabidopsis MA mutations (SBS-only, singletons, TAIR10-anchored).

- `slim/` or `R/meristems2.slim`
  - forward simulation of 3 “generations” of meristem-like growth;
  - tunable parameters: genic mutation multiplier (`M`), selection on nonsynonymous (`S`), dominance (`D`), and separate synonymous / noncoding genic classes;
  - writes `fixed.txt` to `data/slim_out/`, which the R scripts aggregate.

- `Figures/`
  - R scripts export the same panels shown in the manuscript (mutation-rate contrasts and component plots).

## Key idea

1. **Simulate** genic vs intergenic mutations under different (`M`, `S`, `D`) to see how much hypomutation can be created by selection alone.
2. **Re-analyze** Arabidopsis MA mutations against TAIR10, classify mutations (CDS / gene body / intergenic, NS vs S), and **apply the same estimator**.
3. **Use `du_estimate()`** to ask: *given the neutral NS/S and coding fraction, how much of the observed drop in genic mutation rate can be explained by missing nonsynonymous mutations?* The remainder is interpreted as a true mutation-rate effect.

## How to run

1. Run SLiM to populate `data/slim_out/` (or use the provided CSVs).
2. In R:

```r
source("code/functions.R")
muts <- data.table::fread("data/slim_out/SLiM_mutsM.csv")

est <- du_estimate(
  gene_muts    = sum(muts$genic),
  ref_muts     = sum(!muts$genic),
  gene_len     = 4000,
  ref_len      = 4000,
  nonsyn_muts  = sum(muts$Dn),
  syn_muts     = sum(muts$Ds),
  p_cds        = 0.75,
  neutral_ns_s = 2
)
print(est)
```

3. Use the plotting helpers to make PDFs in `Figures/`.

## Notes

- All Arabidopsis coordinates and features are based on **TAIR10**.
- Mutation tables were filtered to **single-base substitutions** and **singletons** to reduce calling/reporting artifacts.
- The selection-correction logic in `du_estimate()` is the reusable part of this repository.

