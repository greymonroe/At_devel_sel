# Forward Evolutionary Simulations

We simulated the accumulation and transmission of de novo mutations through a simplified plant-like developmental life cycle using SLiM (nucleotide-based mode). The goal was to model how variation in genic mutation rate, strength of purifying selection, and dominance affects the probability that mutations ultimately become fixed after repeated rounds of clonal expansion and extreme developmental bottlenecks. The baseline point mutation rate was chosen so that the realized number of mutations per genome per generation was close to empirical estimates for *Arabidopsis thaliana* (i.e. on the order of 10^-9^ to 10^-8^ per site (~1-2/genome) per generation, scaled here to our shorter simulated genome).

## Genome Architecture

We modeled a haploid genome of 8,000 bp. The genome was partitioned into:

- **0–3,999 bp**: non-genic / intergenic region  
- **4,000–7,999 bp**: genic region

This 50:50 split was chosen to approximate a small genome with a substantial functional fraction. Within the genic portion, we further modeled coding versus non-coding genic sequence by assigning different mutation types and sampling weights, such that coding sequence produced nonsynonymous mutations twice as often as synonymous mutations (2:1).

Two genomic element types were defined:

1. **Intergenic element** (`g1`): positions 0–3,999, using mutation type `m1` (neutral) and a baseline Jukes–Cantor mutation model with rate **3.5 × 10^{-6}** per site.
2. **Genic element** (`g2`): positions 4,000–7,999, using a mixture of genic mutation types (see below) and the same model but with its rate multiplied by a parameter **M** (0.25–1.0) to reduce or restore the genic mutation rate relative to intergenic sequence.

Thus, the genic mutation rate was M × 3.5 × 10^-6^. When **M = 1**, genic and intergenic mutated at the same rate; when **M = 0.25**, genic sequence mutated 4× less frequently.

## Mutation Types and Selection

We defined four nucleotide mutation types:

1. **m1 (intergenic)**: dominance = 0.5, fitness effect = 0 (neutral).
2. **m2 (nonsynonymous)**: dominance = **D**, fitness effect drawn from a **gamma distribution** with **mean −S** and **shape = 0.14** (shape taken from Plavskin *et al.* 2024 to capture a realistic right-skewed distribution of deleterious effects). This models deleterious nonsynonymous mutations.
3. **m21 (synonymous)**: dominance = **D2**, fitness effect from a gamma with **mean −S2** and **shape = 0.14**. This allows us to model slightly deleterious or effectively neutral synonymous changes.
4. **m3 (non-coding genic)**: dominance = **D3**, fitness effect from a gamma with **mean −S3** and **shape = 0.14**, representing regulatory or intronic sequence inside genes.

These three genic mutation types were combined in the genic element (`g2`) in the ratio `c(2, 1, 1)` for `(m2, m3, m21)`, so that nonsynonymous mutations were twice as frequent as each of the other two genic categories. All selected mutations used negative means (i.e. purifying selection).

## Parameters Explored

We systematically varied three key axes:

- **Nonsynonymous selection means (S):**

  ```r
  Ss <- c(0, 0.001, 0.01, 0.03, 0.05, 0.1, 0.5, 1)
  ```

  These values let us go from effectively neutral nonsynonymous mutations (S = 0) to quite strongly deleterious ones (S = 1), with the gamma still allowing a range of effects. Note that `0.001, 0.01` represent biologically realistic selection coefficients, but exploring extremely deleterious selection was necessary to induce sufficient selection to test efficacy of selection correction for re-estimating mutation rate heterogeneity.

- **Dominance of nonsynonymous mutations (D):**

  ```r
  Ds <- c(0, 0.1, 0.5, 1)
  ```

  This spans fully recessive (0), partially recessive (0.1), additive (0.5), and fully dominant (1) deleterious alleles."

- **Genic mutation multipliers (M):**

  ```r
  Ms <- seq(from = 0.25, to = 1, length.out = 5)
  # 0.25, 0.4375, 0.625, 0.8125, 1.0
  ```

  This provide true mutation rate reudctions against which observed simulated results can be compared, to directly measure confounding caused by selection along with ability to correct for selection to correctly esimate selection.

For each combination of **S × D × M**, we ran **4,000 independent SLiM simulations** (replicates) to obtain stable counts of fixed mutations, **8 × 4 × 5 = 160 parameter combinations**, **4,000 runs each**, for **640,000 SLiM runs** in total 


## Developmental / Population Schedule

We modeled a repeated plant-like life cycle with clonal expansion followed by severe bottlenecks approximating meristem dynamics and representing single seed descent. This models mutations arising and expand in a small number of meristem cells which can ultimately reach the gametes/zygote and be detected in subsequent generations as fixed within a genome. Each SLiM simulation proceeds for 3 whole plant generations, capturing somatic selection but also allowing for fixation of homozygous mutations that can be detected at the end of each run.

The schedule was:

1. **Generation 1**: found a population of 2 cells and force **clonal reproduction**. These represent two effective embryonic cells.
2. **Generations 2–4**: doubling of the population each generation (2 → 4 → 8 → 16), still clonal. This mimics early expansion to form a vegetative/meristematic zone; any mutation arising here can be propagated clonally and potentially contribute to the germline.
3. **Generations 10–18**: a second period of clonal expansion, again by doubling, to represent proliferation of **gametogenic** cells — i.e. making a large pool of potential reproductive cells that may carry mutations.
4. **Generation 19**: cloning is turned off and sexual reproduction is enabled. This is the "make gametes and zygote" phase.
5. **Generation 20**: the population is **reduced to a single individual** (size = 1). This bottleneck reflect single seed descent used in mutation accumulation experiments.
6. **Generations 21–24 and 30–38**: the same pattern is repeated — re-enable clonal growth, double up, expand from the zygote to build a new meristem. This corresponds to the **second organismal generation**.
7. **Generation 39**: sexual phase again (cloning off).
8. **Generation 40**: bottleneck to 1 again — start of **third organismal generation**.
9. **Generations 41–44 and 50–58**: clonal expansions from that single individual, forming the final meristematic population.
10. **Generation 60**: population set to 1 once more to finalize the lineage.
11. **Generation 61 (late)**: we sampled genomes and wrote only the fixed mutations to disk using:

    ```slim
    sim.outputFixedMutations();
    ```

Any mutation present in that individual at the end is, by definition, fixed. This allows us to count only those mutations that (i) arose under the specified mutation model, (ii) survived purifying selection with the chosen strength and dominance, and (iii) successfully passed through **three** cycles of "grow clonally → sexual phase → bottleneck."

Thus, each SLiM run represents **three organism-like generations of mutation accumulation and transmission through developmental bottlenecks.**

## Output and Postprocessing

Each run produced a file named `fixed.txt` in `data/slim_out/` that contained the list of fixed mutations for that run. Results were parsed and combined in R (see `run_meristems.R` ).

## Rationale

- **Mutation rate baseline**: the intergenic rate (3.5 × 10^-6^) and scaling of genic rates by M were chosen so that, given the number of cell/developmental generations in the schedule, the total number of mutations per simulated genome per organismal generation was similar to observed per-generation mutation rates in *A. thaliana*. This makes the resulting fixed-mutation counts interpretable in terms of realistic plant mutation supply.
- **Gamma with shape 0.14**: using a low shape parameter captures the empirically supported pattern that most new deleterious mutations are of small effect, but a non-trivial tail of larger-effect mutations exists; this was parameterized after Plavskin *et al.* (2024).
- Meristems and gamete cell numbers were based approximately on Burian  *et al.* (2016) and Kakui  *et al.* (2022). Our simulations use more conservative cell number estimates (within order of magnitude of experimental numbers) which should over influence the effect of selection (higher effective population size).
