
library(polymorphology2) #github.com/greymonroe/polymorphology
library(parallel)
library(pbapply)
library(seqinr)
library(ggrepel)


# SLiM --------------------------------------------------------------------


#' Run a single SLiM meristem simulation and parse fixed mutations
#'
#' This wrapper calls an external SLiM script (e.g. `meristems.slim`)
#' with the parameters passed from R. It then reads the resulting
#' `fixed.txt` file (SLiM output) and returns a data.table of fixed
#' mutations, with some helper columns added (genic / Dn / Ds).
#'
#' @param D  Dominance of deleterious non-synonymous mutations (m2)
#' @param N  Population size (currently *not* used in SLiM code, kept for compatibility)
#' @param S  Mean of gamma distribution for negative selection on nonsynonymous
#' @param S2 Mean of gamma distribution for negative selection on synonymous
#' @param S3 Mean of gamma distribution for negative selection on non-coding genic
#' @param M  Mutation-rate multiplier for genic regions
#' @param D2 Dominance of deleterious synonymous mutations
#' @param D3 Dominance of deleterious non-coding mutations
#' @param slim_file Path to the SLiM script to run
#' @param fixed_file Path to the SLiM output file with fixed mutations
#'
#' @return data.table with fixed mutations, or NULL if no mutations were fixed
run_slim_meristems <- function(
    D = 0.5,
    N = 100,       # NOTE: kept for compatibility; may be ignored by the .slim file
    S = 0.01,
    S2 = 0.01,
    S3 = 0.01,
    M = 0.5,
    D2 = 0.5,
    D3 = 0.5,
    slim_file  = "meristems.slim",
    fixed_file = "fixed.txt"
) {

  ## 1) build and run command
  # SLiM lets us pass variables with -d
  cmd <- paste0(
    "slim ",
    "-d D=",  D,  " ",
    "-d N=",  N,  " ",
    "-d S=",  S,  " ",
    "-d M=",  M,  " ",
    "-d S2=", S2, " ",
    "-d D2=", D2, " ",
    "-d S3=", S3, " ",
    "-d D3=", D3, " ",
    slim_file
  )

  # run SLiM quietly (you had ignore.stdout=TRUE)
  system(cmd, ignore.stdout = TRUE)

  ## 2) make sure the file exists and has content
  if (!file.exists(fixed_file)) {
    # SLiM didn't produce the expected file
    return(NULL)
  }

  lines <- readLines(fixed_file)
  # your original check: if <3 lines, skip
  if (length(lines) < 3) {
    return(NULL)
  }

  ## 3) read the fixed mutations table
  # original code used: skip = 2
  # that’s typical for SLiM fixed mutation log with 2 header lines
  fixed <- fread(fixed_file, skip = 2, verbose = FALSE)

  # you had: if(nrow(fixed)>0 & ncol(fixed)>1) ...
  if (nrow(fixed) == 0L || ncol(fixed) <= 1L) {
    return(NULL)
  }

  # name columns explicitly to something sensible
  # adjust if your slim output changes
  setnames(
    fixed,
    c(
      "Mutation_ID",  # id
      "X1",           # often generation or age, depends on SLiM logFormat
      "Mutation_Type",
      "POS",
      "Fitness",
      "Dominance",
      "Subpopulation",
      "Gen_mutated",
      "Gen_fixed",
      "ALT"
    )
  )

  ## 4) add helper columns
  # your original logic:
  # genic = POS > 4000
  fixed[, genic := POS > 4000]

  # classification by mutation type (your code used "m2" and "m21")
  fixed[, Dn := Mutation_Type == "m2"]   # deleterious nonsyn
  fixed[, Ds := Mutation_Type == "m21"]  # deleterious syn
  # NOTE: if your SLiM file uses different names, change here

  return(fixed)
}


#' Run many SLiM meristem simulations and bind results
#'
#' @param D,S,S2,S3,M,D2,D3 See `run_slim_meristems()`
#' @param N  Population size (passed through; may be unused by SLiM script)
#' @param n_reps Number of replicate SLiM runs to perform
#' @param slim_file Path to `.slim` script
#' @param fixed_file Name of output file to read (per run). If SLiM overwrites
#'   the same file every run, this is fine; we read it immediately each time.
#'
#' @return data.table of all fixed mutations across all runs, with parameters
#'         and replicate id attached; or NULL if none succeeded
bulk_slim <- function(
    D,
    N,
    S,
    S2,
    S3,
    M,
    D2,
    D3,
    n_reps = 10L,
    slim_file  = "meristems.slim",
    fixed_file = "fixed.txt"
) {

  sim_list <- lapply(seq_len(n_reps), function(rep_i) {
    out <- run_slim_meristems(
      D = D, N = N,
      S = S, S2 = S2, S3 = S3,
      M = M, D2 = D2, D3 = D3,
      slim_file  = slim_file,
      fixed_file = fixed_file
    )

    # if this run produced no fixed muts, skip it
    if (is.null(out)) return(NULL)

    # tag with replicate id
    out[, rep := rep_i]
    out
  })

  # drop NULLs
  sim_list <- Filter(Negate(is.null), sim_list)

  if (length(sim_list) == 0L) {
    return(NULL)
  }

  all_fixed <- rbindlist(sim_list, fill = TRUE)

  # attach parameter values to every row for later plotting / modeling
  all_fixed[, `:=`(
    D = D,
    S = S,
    S2 = S2,
    S3 = S3,
    N = N,
    M = M,
    D2 = D2,
    D3 = D3
  )]

  return(all_fixed)
}


## ------------------------------------------------------------
## 2. Plot helpers
## ------------------------------------------------------------
plot_est <- function(est, variable, yname) {
  est2 <- est[D %in% c(0, 1)]
  est2$variable <- est2[[variable]]

  ggplot(est2, aes(x = M, y = variable, col = factor(D), group = D)) +
    geom_point(shape = 20) +
    geom_line() +
    facet_grid(S2 ~ S) +
    scale_y_continuous(name = yname) +
    theme_bw(base_size = 6) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1),
      strip.background = element_blank()
    )
}

plot_estSonly <- function(est, variable, yname) {
  est2 <- est[S2 == 0 & D %in% c(0, 1)]
  est2$variable <- est2[[variable]]

  ggplot(est2, aes(x = M, y = variable, col = factor(D), group = D)) +
    geom_point(shape = 20, size = 0.5) +
    geom_line() +
    facet_grid(~S) +
    scale_color_manual(values = c("black", "red")) +
    theme_bw(base_size = 6) +
    theme(
      legend.position   = "none",
      axis.text.x       = element_blank(),
      axis.ticks.x      = element_blank(),
      axis.title.x      = element_blank(),
      strip.background  = element_blank(),
      strip.text        = element_blank(),
      plot.background   = element_blank(),
      panel.background  = element_blank()
    )
}

#' Estimate genic mutation rate relative to a reference and correct for NS selection
#'
#' Compares an observed genic mutation rate (mutations / genic bp) to a reference
#' mutation rate (mutations / reference bp), then estimates how much of the
#' difference could be explained by a deficit of nonsynonymous mutations.
#'
#' All class-specific deltas are defined as observed / expected (i.e. multipliers).
#'
#' @param gene_muts Integer. Total mutations in the focal genic region.
#' @param ref_muts Integer. Total mutations in the reference region
#'   (e.g. intergenic or non-focal genes).
#' @param gene_len Numeric. Total bp of the focal genic region.
#' @param ref_len Numeric. Total bp of the reference region.
#' @param nonsyn_muts Integer. Observed nonsynonymous mutations in the focal genes.
#' @param syn_muts Integer. Observed synonymous mutations in the focal genes.
#' @param p_cds Numeric in (0,1]. Fraction of focal gene sequence that is CDS.
#' @param neutral_ns_s Numeric > 0. Neutral NS/S ratio to use for expectations.
#' @param return_reduction Logical. If TRUE, also return 1 - (multiplier) columns.
#'
#' @return data.table
du_estimate <- function(
    gene_muts,
    ref_muts,
    gene_len,
    ref_len,
    nonsyn_muts,
    syn_muts,
    p_cds,
    neutral_ns_s,
    return_reduction = TRUE
) {
  # 1) observed rates
  mu_g_obs <- gene_muts / gene_len       # genic observed rate
  mu_ref   <- ref_muts  / ref_len        # reference observed rate

  # overall genic/reference multiplier
  duO <- mu_g_obs / mu_ref               # <1 = hypomutation

  # 2) split genic mutations
  noncoding_muts <- gene_muts - syn_muts - nonsyn_muts

  # 3) expected class proportions under neutrality
  p_ns_exp <- p_cds * neutral_ns_s / (neutral_ns_s + 1)
  p_s_exp  <- p_cds - p_ns_exp
  p_nc_exp <- 1 - p_cds

  # 4) expected counts if genes mutated like reference
  gene_exp  <- mu_ref * gene_len
  ns_exp_mu <- gene_exp * p_ns_exp
  s_exp_mu  <- gene_exp * p_s_exp
  nc_exp_mu <- gene_exp * p_nc_exp

  # 5) independent NS expectation from observed S
  ns_exp_s <- syn_muts * neutral_ns_s

  # 6) correction term: how much of duO is just missing NS
  RNS <- (p_ns_exp * (ns_exp_s - nonsyn_muts) / ns_exp_mu)

  # 7) selection-corrected overall multiplier
  duEst <- duO + RNS

  # 8) class-specific multipliers (observed / expected)
  dmu_ns <- nonsyn_muts    / ns_exp_mu
  dmu_s  <- syn_muts       / s_exp_mu
  dmu_nc <- noncoding_muts / nc_exp_mu

  # build output first (simpler to paste into console)
  out <- data.table::data.table(
    duO       = duO,
    duEst     = duEst,
    RNS       = RNS,
    mu_g_obs  = mu_g_obs,
    mu_ref    = mu_ref,
    gene_exp  = gene_exp,
    ns_exp_mu = ns_exp_mu,
    s_exp_mu  = s_exp_mu,
    nc_exp_mu = nc_exp_mu,
    dmu_ns    = dmu_ns,
    dmu_s     = dmu_s,
    dmu_nc    = dmu_nc
  )

  if (isTRUE(return_reduction)) {
    out[, duO_reduction       := 1 - duO]
    out[, duEst_reduction     := 1 - duEst]
    out[, dmu_ns_reduction    := 1 - dmu_ns]
    out[, dmu_s_reduction     := 1 - dmu_s]
    out[, dmu_nc_reduction    := 1 - dmu_nc]
  }

  out
}

# Estimating expected neutral Ns/S ----------------------------------------

#' Convert CDS FASTA into a codon-level mutation table
#'
#' For each coding sequence (CDS) in a FASTA (read with **seqinr**), enumerate
#' all possible single-nucleotide mutations at each codon position, and annotate
#' whether the resulting codon change is synonymous (S) or nonsynonymous (N).
#'
#' @param CDS A CDS FASTA object as returned by seqinr::read.fasta().
#'
#' @return A data.table with columns REF, ALT, effect (S/N), and gene.
#' @export
CDS_codon_mutation_table <- function(CDS) {
  stopifnot(requireNamespace("data.table", quietly = TRUE))
  DT <- data.table::data.table

  # all single-nuc substitutions
  muts <- DT(table(REF = c("A", "T", "C", "G"),
                   ALT = c("A", "T", "C", "G")))[REF != ALT]

  # codon → AA table from seqinr, uppercased
  codons <- DT(seqinr::SEQINR.UTIL$CODON.AA)
  codons[, CODON := toupper(CODON)]

  # precompute: for every codon, what happens if we mutate pos 1/2/3 to any other base?
  CODON_mutations_table <- data.table::rbindlist(
    lapply(seq_len(nrow(codons)), function(i) {
      codon_i <- codons$CODON[i]

      # split current codon into 3 bases
      codon_bases <- strsplit(codon_i, "")[[1]]
      codon_dt <- DT(REF = codon_bases,
                     POS = 1:3,
                     REFCODON = codon_i)

      # join with all possible nucleotide changes for that REF base
      codon_dt <- merge(codon_dt, muts, by = "REF", allow.cartesian = TRUE)

      # build the mutated codon string
      codon_dt[, ALTCODON := vapply(seq_len(.N), function(j) {
        new_cod <- codon_bases
        new_cod[codon_dt$POS[j]] <- codon_dt$ALT[j]
        paste0(new_cod, collapse = "")
      }, character(1L))]

      # look up AAs
      codon_dt[, REFAA := codons$AA[match(REFCODON, codons$CODON)]]
      codon_dt[, ALTAA := codons$AA[match(ALTCODON, codons$CODON)]]

      # classify effect
      codon_dt[, effect := ifelse(REFAA == ALTAA, "S", "N")]

      codon_dt
    })
  )

  # iterate over each CDS entry
  pb <- utils::txtProgressBar(min = 0, max = length(CDS), style = 3)
  res <- data.table::rbindlist(
    lapply(seq_along(CDS), function(i) {
      seq_i <- toupper(CDS[[i]])
      gene_i <- names(CDS)[i]

      # must be multiple of 3
      if ((length(seq_i) %% 3) != 0L) {
        message("CDS '", gene_i, "' is not a multiple of 3; skipping.")
        utils::setTxtProgressBar(pb, i)
        return(NULL)
      }

      # translate once
      aa_seq <- seqinr::translate(seq_i)

      # build base-level table for this CDS
      cds_dt <- DT(
        POS     = rep(1:3, times = length(seq_i) / 3),
        DNAPOS  = seq_along(seq_i),
        REF     = seq_i,
        AA      = rep(aa_seq, each = 3)
      )

      # make codon strings for each triplet
      codon_strings <- vapply(seq(1, length(seq_i), by = 3), function(start) {
        paste0(seq_i[start:(start + 2)], collapse = "")
      }, character(1L))

      cds_dt[, REFCODON := rep(codon_strings, each = 3)]

      # expand to all possible nucleotide changes at each position
      cds_mut <- merge(cds_dt, muts, by = "REF", allow.cartesian = TRUE)
      data.table::setorder(cds_mut, DNAPOS)

      # annotate effect by joining to precomputed codon mutation table
      cds_mut <- merge(
        cds_mut,
        CODON_mutations_table,
        by = c("REFCODON", "REF", "POS", "ALT"),
        all.x = TRUE,
        sort = FALSE
      )

      cds_mut[, gene := gene_i]

      utils::setTxtProgressBar(pb, i)

      # return minimal columns as in your original
      cds_mut[, .(REF, ALT, effect, gene)]
    }),
    fill = TRUE
  )
  close(pb)

  res
}


mutation_probs<-function(muts, genome, cds=NULL){
  message("Calculating base pair frequencies from genome...")
  bp_freq<-rbindlist(lapply(1:length(genome), function(i){
    message(paste("\tChr",i))
    bp_freq<-data.table(table(REF=genome[[i]]))[REF!="n"]
    bp_freq$REF<-toupper(bp_freq$REF)
    return(bp_freq)
  }))
  bp_freq<-bp_freq[REF %in% c("A","T","C","G"),.(N=sum(N)), by="REF"]

  mutations<-data.table((table(REF=muts$REF, ALT=muts$ALT)))[REF!=ALT]
  mutations<-merge(mutations, bp_freq, by="REF")
  mutations$rawPROB=(mutations$N.x/(mutations$N.y))

  if(!is.null(cds)){
    message("Calculating base pair frequencies from CDS...")
    bp_freqCDS<-rbindlist(lapply(1:length(cds), function(i){

      bp_freq<-data.table(table(REF=cds[[i]]))[REF!="n"]
      bp_freq$REF<-toupper(bp_freq$REF)
      return(bp_freq)
    }))
    bp_freqCDS<-bp_freqCDS[REF %in% c("A","T","C","G"),.(CDSN=sum(N)), by="REF"]

    mutations<-merge(mutations, bp_freqCDS, by="REF")
    mutations$rawPROB=(mutations$N.x/(mutations$N.y-mutations$CDSN))
  }
  return(mutations)
}

sample_neutrals <- function(neutrals, N, i) {
  pb <- txtProgressBar(min = 0, max = i, style = 3)
  out <- numeric(i)
  for (k in seq_len(i)) {
    samp <- neutrals[sample.int(length(neutrals), N)]
    ns   <- sum(samp == "N")
    s    <- sum(samp == "S")
    out[k] <- ns / (s + 1e-9)  # avoid div/0
    setTxtProgressBar(pb, k)
  }
  close(pb)
  out
}

sample_neutrals_groups<-function(group="all", reps){
  if(group!="all"){

    sub<-mutations[,c("mutations_effect", group), with=F]
    colnames(sub)[2]<-"group"
    sub_sum<-sub[,.(NS=sum(mutations_effect=="Non-Syn"), S=sum(mutations_effect=="Syn")), by=.(inout=ifelse(group,"in","out"))][order(inout, decreasing = T)]
    sub_sum$NsS<-sub_sum$NS/sub_sum$S
    group_neutrals<-sample_neutrals(effects, sub_sum$NS[2]+sub_sum$S[2], reps)
    outgroup_neutrals<-sample_neutrals(effects, sub_sum$NS[1]+sub_sum$S[1], reps)
    dt<-data.table(group=group, s=c(group_neutrals, outgroup_neutrals), inout=rep(c("in","out"), each=reps))
    dt$NsSO<-sub_sum$NsS[match(dt$inout, sub_sum$inout)]
  } else {
    sub_sum<-mutations[,.(NS=sum(mutations_effect=="Non-Syn"), S=sum(mutations_effect=="Syn"))]
    group_neutrals<-sample_neutrals(effects, sub_sum$NS[1]+sub_sum$S[1], reps)
    dt<-data.table(group=group, s=c(group_neutrals), inout=rep(c("in"), each=reps))
    dt$NsSO<-sub_sum$NS/sub_sum$S
  }
  return(dt)
}

# helper: count NS, S, and NS/S from a data.table of mutations
get_ns_s <- function(dt) {
  NS <- sum(dt$mutations_effect == "Non-Syn", na.rm = TRUE)
  S  <- sum(dt$mutations_effect == "Syn",     na.rm = TRUE)
  list(NS = NS, S = S, total = NS + S, ratio = NS / S)
}

# Mutation Accumulation Experiments ---------------------------------------


# simple helper: bootstrap one source
bootstrap_src_mutations <- function(dt_src, B = 1000L) {
  n <- nrow(dt_src)
  # preallocate
  out <- vector("list", B)
  for (b in seq_len(B)) {
    samp <- dt_src[sample.int(n, n, replace = TRUE)]
    NS   <- sum(samp$mutations_effect == "Non-Syn", na.rm = TRUE)
    S    <- sum(samp$mutations_effect == "Syn",     na.rm = TRUE)
    out[[b]] <- data.table(
      nss   = NS / S,
      genic = sum(samp$gene_body, na.rm = TRUE) / nrow(samp)
    )
  }
  rbindlist(out)
}

plot_mutsrc_boot <- function(mutations,
                             B = 1000L,
                             neutral_nss = 2.74,
                             genome_genic = 0.51) {

  # 1) observed per source
  summary_obs <- mutations[
    ,
    .(
      NS    = sum(mutations_effect == "Non-Syn", na.rm = TRUE),
      S     = sum(mutations_effect == "Syn",     na.rm = TRUE),
      genic = sum(gene_body, na.rm = TRUE) / .N,
      N     = .N
    ),
    by = src
  ]
  summary_obs[, nss := NS / S]

  # 2) bootstrap per source
  boot_list <- lapply(summary_obs$src, function(s) {
    dt_src <- mutations[src == s]
    bt     <- bootstrap_src_mutations(dt_src, B = B)
    bt[, src := s]
    bt
  })
  boot_dt <- rbindlist(boot_list, fill = TRUE)

  # 3) CIs
  ci_dt <- boot_dt[
    ,
    .(
      nss_lo   = quantile(nss,   0.025, na.rm = TRUE),
      nss_hi   = quantile(nss,   0.975, na.rm = TRUE),
      genic_lo = quantile(genic, 0.025, na.rm = TRUE),
      genic_hi = quantile(genic, 0.975, na.rm = TRUE)
    ),
    by = src
  ]

  plot_dt <- merge(summary_obs, ci_dt, by = "src", all.x = TRUE)

  # 4) plot
  p <- ggplot() +
    geom_errorbar(
      data = plot_dt,
      aes(x = nss, ymin = genic_lo, ymax = genic_hi),
      width = 0,
      linewidth = 0.25,
      col = "gray90"
    ) +
    geom_errorbar(
      data = plot_dt,
      aes(y = genic, xmin = nss_lo, xmax = nss_hi),
      orientation = "y",
      height = 0,
      linewidth = 0.25,
      col = "gray90"
    ) +
    geom_point(
      data = plot_dt,
      aes(x = nss, y = genic),
      size = 1.5,
      col  = "gray40"
    ) +
    ggrepel::geom_text_repel(
      data = plot_dt,
      aes(x = nss, y = genic, label = src),
      size = 2.2,
      max.overlaps = 40,
      color = "gray15"
    ) +
    # vertical "neutral" line
    geom_vline(
      xintercept = neutral_nss,
      linetype   = "dashed",
      linewidth  = 0.3,
      color      = "darkred"
    ) +
    annotate(
      "text",
      x     = neutral_nss,
      y     = min(plot_dt$genic_lo, na.rm = TRUE) - 0.015,
      label = "Neutral expectation\n(no selection)",
      color = "darkred",
      size  = 2.2,
      vjust = 0,
      hjust = -0.1
    ) +
    # horizontal "genome expectation" line
    geom_hline(
      yintercept = genome_genic,
      linetype   = "dashed",
      linewidth  = 0.3,
      color      = "darkred"
    ) +
    annotate(
      "text",
      x     = max(plot_dt$nss, na.rm = TRUE),
      y     = genome_genic,
      label = "Genome expectation\n(no mutation bias)",
      color = "darkred",
      size  = 2.2,
      hjust = 0,
      vjust = 1.2
    ) +
    scale_y_continuous(name = "Fraction of mutations in genes") +
    scale_x_log10(name = "Nonsynonymous / Synonymous mutations") +
    labs(
      title = "Selection and mutation bias signals in Arabidopsis MA datasets\n(bootstrapped 95% CI)"
    ) +
    theme_classic(base_size = 6)+
    theme(panel.border = element_rect(linewidth = 0.5))

  return(p)
}



bootstrap_src_mutations <- function(dt_src, B = 1000L) {
  # dt_src: data.table of mutations for a single source (one src level)
  n <- nrow(dt_src)
  res <- vector("list", B)

  for (b in seq_len(B)) {
    idx <- sample.int(n, n, replace = TRUE)
    x   <- dt_src[idx]

    NSb    <- sum(x$mutations_effect == "Non-Syn", na.rm = TRUE)
    Sb     <- sum(x$mutations_effect == "Syn",     na.rm = TRUE)
    genicb <- sum(x$gene_body, na.rm = TRUE) / n

    res[[b]] <- list(
      nss   = if (Sb > 0) NSb / Sb else NA_real_,
      genic = genicb
    )
  }

  out <- rbindlist(res)
  out
}

make_src_trimer_plots <- function(mutations, genome_triN) {
  srcs <- sort(unique(mutations$src))

  tri_plots <- lapply(srcs, function(s) {
    x <- mutations[src == s & !is.na(trimer)]

    tri_obj <- plot_tricontexts(
      x$trimer,
      full        = TRUE,
      trimer_freq = genome_triN
    )

    tri_obj$plot + ggtitle(s)
  })

  names(tri_plots) <- srcs
  tri_plots
}

make_genic_intergenic_trimer_plots <- function(mutations,
                                               genes_triN,
                                               intergenic_triN) {
  # genic contexts
  tri_genic <- plot_tricontexts(
    mutations[gene_body == TRUE]$trimer,
    full        = TRUE,
    trimer_freq = genes_triN
  )

  # intergenic contexts
  tri_nongenic <- plot_tricontexts(
    mutations[gene_body == FALSE]$trimer,
    full        = TRUE,
    trimer_freq = intergenic_triN
  )

  # pull context tables
  ct_genic    <- tri_genic$context_table[, .(context_only, mut, N_genic = N)]
  ct_intergen <- tri_nongenic$context_table[, .(context_only, mut, N_intergenic = N)]

  # merge and make ratio
  ct_ratio <- merge(
    ct_genic,
    ct_intergen,
    by = c("context_only", "mut"),
    all = TRUE
  )
  ct_ratio[is.na(N_genic),       N_genic       := 0]
  ct_ratio[is.na(N_intergenic),  N_intergenic  := 0]
  ct_ratio[, ratio := (N_intergenic + 1e-6) / (N_genic + 1e-6)]

  # plots
  tri_genic_plot <- tri_genic$plot + ggtitle("Genic mutations")
  tri_intergenic_plot <- tri_nongenic$plot + ggtitle("Intergenic mutations")

  tri_ratio_plot <- ggplot(
    ct_ratio,
    aes(x = context_only, y = ratio, fill = mut)
  ) +
    geom_bar(stat = "identity", width = 0.5) +
    facet_grid(. ~ mut, scales = "free_x", space = "free_x") +
    scale_x_discrete(name = "Context") +
    scale_y_continuous(name = "Intergenic / genic") +
    theme_classic(base_size = 6) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_blank()
    ) +
    ggtitle("Relative enrichment (intergenic / genic)") +
    scale_fill_manual(
      values = c("cyan3", "black", "red4", "gray", "green3", "pink3"),
      guide = "none"
    )

  list(
    genic      = tri_genic_plot,
    intergenic = tri_intergenic_plot,
    ratio      = tri_ratio_plot
  )
}

# annotate: is each mutation inside a given gene/feature set?
annot_mut_geneset <- function(mutations, gene_set) {
  hits <- features_in_sites(gene_set, mutations)
  mutations$inset <- hits$overlaps[match(mutations$ID, hits$ID)]
  mutations$inset
}

# (optional / legacy) how many CDS models does this mutation overlap?
# kept for compatibility; can be dropped if unused
annot_mod_number <- function(i, mutations, CDS) {
  genic <- mutations$CDS[i]
  if (!genic) return(0L)

  row <- mutations[i]
  overlap <- sites_in_features(CDS[CHROM == row$CHROM], row, mode = "any")[any == TRUE]
  mod <- unique(CDS[ID %in% overlap$ID]$model)
  length(mod)
}


#' Annotate a single mutation with its coding effect
#'
#' Given one mutation ID and the full mutation table, this function finds the
#' overlapping CDS model on the same chromosome, rebuilds the transcript
#' sequence from the CDS exons, applies the nucleotide change, and reports
#' whether the change is synonymous or nonsynonymous. Non-coding mutations and
#' cases where we can't map the CDS are returned with informative `effect`
#' labels.
#'
#' @param mut_id Integer (scalar). The mutation ID (i.e. the value in
#'   `mutations$ID`) to annotate.
#' @param mutations `data.table`. Table of mutations containing at least the
#'   columns `ID`, `CHROM`, `POS`, `REF`, `ALT`, and a logical `CDS` column
#'   indicating whether the site is inside a CDS feature.
#' @param CDS `data.table`. CDS features from the GFF, already filtered to the
#'   representative models. Must contain at least `CHROM`, `START`, `STOP`,
#'   `DIRECTION` (strand), `model`, and `ID`.
#' @param cds List of CDS sequences as read by **seqinr** (e.g.
#'   `seqinr::read.fasta(...)`), named by transcript/model. Used as the
#'   reference coding sequence for translation.
#' @param genome_fasta List of chromosome sequences (FASTA) indexed by chromosome
#'   name (e.g. `"1"`, `"2"`, …).
#'
#' @return A `data.table` with one row and columns:
#'   - `CHROM`, `POS`, `ID`
#'   - `GENE_model` (the model/transcript hit)
#'   - `REF_AA`, `ALT_AA`
#'   - `effect` ∈ {`"Syn"`, `"Non-Syn"`, `"nonCDS"`, `"CDS_model_not_found"`,
#'     `"CDS_pos_BP_mismatch"`}
#'
#' @details
#' On minus-strand CDS, genomic positions are reversed and the alleles are
#' reverse-complemented using your local `revcomp()` helper.
annot_mut_effect <- function(mut_id, mutations, CDS, cds, genome_fasta) {
  # 1) pull the mutation
  row   <- mutations[ID == mut_id]
  CHROM <- row$CHROM
  POS   <- row$POS
  REF   <- row$REF
  ALT   <- row$ALT
  coding <- row$CDS
  this_chrom <- row$CHROM

  # 2) non-coding: return fast
  if (!isTRUE(coding)) {
    return(data.table(
      CHROM      = CHROM,
      POS        = POS,
      ID         = mut_id,
      GENE_model = NA_character_,
      REF_AA     = NA_character_,
      ALT_AA     = NA_character_,
      effect     = "nonCDS"
    ))
  }

  # 3) find CDS feature(s) on this chromosome
  overlap <- sites_in_features(CDS[CHROM == this_chrom], row, mode = "any")[any == TRUE]
  mod     <- unique(CDS[ID %in% overlap$ID]$model)

  if (length(mod) == 0L) {
    return(data.table(
      CHROM      = CHROM,
      POS        = POS,
      ID         = mut_id,
      GENE_model = NA_character_,
      REF_AA     = NA_character_,
      ALT_AA     = NA_character_,
      effect     = "CDS_model_not_found"
    ))
  }

  # 4) collect exons for that model
  CDS_model <- CDS[model == mod][order(START)]
  DIR       <- CDS_model$DIRECTION[1]

  # 5) flatten exons -> genomic CDS positions
  CDS_POS <- unlist(lapply(seq_len(nrow(CDS_model)), function(j) {
    CDS_model$START[j]:CDS_model$STOP[j]
  }))

  # 6) reference CDS from fasta (expected transcript sequence)
  BPs <- toupper(cds[[mod]])

  # 7) reference genome only at CDS positions
  genome_chr <- genome_fasta[[CHROM]]       # no toupper on whole chr
  genome_BPs <- genome_chr[CDS_POS]
  genome_AA  <- seqinr::translate(toupper(genome_BPs))

  # 8) check REF base
  REF_genome <- toupper(genome_chr[POS])
  if (REF_genome != toupper(REF)) {
    warning(mut_id, " REF_genome_mismatch: genome=", REF_genome, " vs REF=", REF)
  }

  # 9) length must match to map codons
  if (length(CDS_POS) != length(BPs)) {
    return(data.table(
      CHROM      = CHROM,
      POS        = POS,
      ID         = mut_id,
      GENE_model = mod,
      REF_AA     = NA_character_,
      ALT_AA     = NA_character_,
      effect     = "CDS_pos_BP_mismatch"
    ))
  }

  # 10) build per-base CDS table
  cds_dt <- data.table(
    cPOS = CDS_POS,
    BP   = toupper(BPs)
  )

  # 11) minus strand: reverse order + use your revcomp()
  if (DIR == "-") {
    cds_dt$cPOS <- rev(cds_dt$cPOS)
    REF <- revcomp(toupper(REF))
    ALT <- revcomp(toupper(ALT))
  }

  # 12) original AA (per codon) repeated per 3 bp
  cds_dt$AA <- rep(seqinr::translate(cds_dt$BP), each = 3)

  # 13) apply mutation
  cds_dt$BP2 <- cds_dt$BP
  cds_dt$BP2[cds_dt$cPOS == POS] <- ALT

  # 14) translate mutated CDS
  cds_dt$newAA <- rep(seqinr::translate(cds_dt$BP2), each = 3)

  # 15) AA actually changed
  OGAA  <- cds_dt[cPOS == POS]$AA
  NEWAA <- cds_dt[cPOS == POS]$newAA

  effect <- if (all(cds_dt$AA == cds_dt$newAA)) "Syn" else "Non-Syn"

  data.table(
    CHROM      = CHROM,
    POS        = POS,
    ID         = mut_id,
    GENE_model = mod,
    REF_AA     = OGAA,
    ALT_AA     = NEWAA,
    effect     = effect
  )
}


calculate_dUest <- function(loci_muts, other_muts, loci_genes, loci_cds, other_L=NULL, DnDsN=NULL, genes=NULL) {

  #message("Calculating genes_list")
  genes_list <- loci_genes[, {
    positions = unlist(lapply(1:.N, function(i) START[i]:STOP[i]))
    unique_positions = unique(positions)
    list(unique_positions = list(unique_positions))
  }, by = CHROM]$unique_positions

  #message("Calculating total number of basepairs in genes (LG)")
  LG <- sum(unlist(lapply(genes_list, length)))

  #message("Calculating CDS_loci_list")
  CDS_loci_list <- loci_cds[, {
    positions = unlist(lapply(1:.N, function(i) START[i]:STOP[i]))
    unique_positions = unique(positions)
    list(unique_positions = list(unique_positions))
  }, by = CHROM]$unique_positions

  #message("Calculating total number of basepairs in CDS (LG_cdsL)")
  LG_cdsL <- sum(unlist(lapply(CDS_loci_list, length)))

  loci_nonCDSL <- LG - LG_cdsL

  if (is.null(DnDsN)) {
   # message("Calculating DnDsN")
    DnDsN <- sum(other_muts$mutations_effect == "Non-Syn") / sum(other_muts$mutations_effect == "Syn")
  }

  if (is.null(other_L)) {
    #message("Calculating genes_list for other_L")
    genes_list <- genes[!gene %in% loci_genes$gene, {
      positions = unlist(lapply(1:.N, function(i) START[i]:STOP[i]))
      unique_positions = unique(positions)
      list(unique_positions = list(unique_positions))
    }, by = CHROM]$unique_positions

   # message("Calculating total number of basepairs in other_L")
    other_L <- sum(unlist(lapply(genes_list, length)))
  }

  #message("Calculating loci_mutN")
  loci_mutN <- nrow(loci_muts)

  #message("Calculating NCO")
  NCO <- sum(loci_muts$CDS == FALSE)

  #message("Calculating other_mutN")
  other_mutN <- nrow(other_muts)

  #message("Calculating u_O and µI_O")
  u_O <- loci_mutN / LG
  µI_O <- other_mutN / other_L

  #message("Calculating duO")
  duO <- 1 - (u_O / µI_O)

 # message("Performing chi-squared test for duO")
  duO_chi <- chisq.test(matrix(c(loci_mutN, other_mutN, LG, other_L), nrow = 2))

  #message("Calculating SO and NSO")
  SO <- sum(loci_muts$mutations_effect == "Syn")
  NSO <- sum(loci_muts$mutations_effect == "Non-Syn")

  #message("Calculating DnDsO")
  DnDsO <- NSO / SO

  #message("Calculating PCDS and PNS_exp")
  PCDS <- LG_cdsL / LG
  PNS_exp <- (DnDsN) / (DnDsN + 1) * PCDS

  #message("Calculating NSO_exp_µ and NSO_exp_S")
  NSO_exp_µ <- µI_O * LG * PNS_exp
  NSO_exp_S <- SO * (DnDsN)

  #message("Calculating RNS")
  RNS <- -(PNS_exp * (NSO_exp_S - NSO) / NSO_exp_µ)

  #message("Calculating duEst")
  duEst <- duO + RNS

  #message("Calculating µNS")
  dµNS <- (NSO_exp_µ - NSO) / NSO_exp_µ

  #message("Calculating PS_exp and SO_exp_µ")
  PS_exp <- PCDS - PNS_exp
  SO_exp_µ <- µI_O * LG * PS_exp

  #message("Calculating µS")
  dµS <- (SO_exp_µ - SO) / SO_exp_µ

 # message("Calculating PNC_exp and NCO_exp_µ")
  PNC_exp <- 1 - PCDS
  NCO_exp_µ <- µI_O * LG * PNC_exp

  #message("Calculating µNC")
  dµNC <- (NCO_exp_µ - NCO) / NCO_exp_µ

  #message("Preparing output data.table")
  out <- data.table(DnDsO, DnDsN, loci_mutN, other_mutN, LG, other_L, u_O, µI_O, duO, duO_chi_p = duO_chi$p.value, duEst,
                    RNS, NSO, NSO_exp_µ, dµNS,
                    SO, SO_exp_µ, dµS,
                    NCO, NCO_exp_µ, dµNC, PCDS)

  #message("Returning output")
  return(out)
}


