# Delineating Mutation Bias and Selection during Plant Development

This repository contains code and simulation material from **Monroe et al. 2025**. The goal is to quantify genic hypomutation in plants and to separate **true mutation-rate reductions** from **apparent reductions caused by selection on nonsynonymous sites**.

## Contents

- `code/`
  - `functions.R` — main R helpers for wrangling, loading, plotting, and basic stats.
    - `neutral_NSS_estimate()` — estimates a neutral nonsynonymous/synonymous (NS/S) ratio from:
      1. observed mutation spectrum (REF→ALT),
      2. CDS codon usage,
      3. genome base composition.
      Returns a list with (i) spectra and codon-weighted mutation probabilites, (ii) effect-level weights, and (iii) the neutral NS/S ratio.
    - `du_estimate()` Estimate genic mutation rate relative to a reference and correct for NS selection. Returns ∆µ~Est~, and component ∆µ~NS~, ∆µ~S~, ∆µ~NC

- `slim/`
  - `meristems.slim` — forward SLiM simulation of 3 “generations” of meristem-like growth and reproduction; parameters for genic mutation multiplier (`M`, true ∆µ), nonsynonymous selection strength (`S`), dominance (`D`), and separate synonymous / noncoding genic classes.

- `data/`
  - mutation datasets, SLiM outputs, and reference genomes (not all tracked in git due to size).

-`tables/`
  - `mutations_Athal_annotated.csv`
  - `SLiM_out.csv`
  - `mutations_Athal_annotated.csv`
  - `codon_mut_counts.csv`

## Requirements

- R (≥ 4.2 recommended)
- packages: `data.table`, `seqinr`, `ggplot2`, `ggrepel`, `polymorphology2`
- for SLiM runs: SLiM 5
