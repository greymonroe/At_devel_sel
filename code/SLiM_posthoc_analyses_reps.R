source("code/functions.R")   # must contain du_estimate(), plot_estSonly(), etc.
library(data.table)
library(ggplot2)
library(cowplot)

## constants taken from your SLiM setup
neutral_ns_s_global <- 2
p_cds_global        <- 3/4

## read SLiM mutations (must have S, D, M, S2, S3, D2, D3, rep, genic, Dn, Ds)
muts <- fread("data/slim_out/SLiM_mutsM.csv")

## helper: mean + se
summ_se <- function(x) {
  m  <- mean(x, na.rm = TRUE)
  se <- sd(x,  na.rm = TRUE) / sqrt(sum(is.finite(x)))
  list(mean = m, se = se)
}

## 1) compute per-replicate summaries (one row per param combo per rep)
##    we keep the exact same 4000/4000 genome layout you used
rep_level <- muts[
  ,
  .(
    genespace    = 4000,
    nongenspace  = 4000,
    gene         = sum(genic),
    nongenic     = sum(!genic),
    Dn           = sum(Dn),
    Ds           = sum(Ds)
  ),
  by = .(S, D, M, S2, S3, D2, D3, rep)
]

## 2) run du_estimate on every replicate
rep_est_list <- lapply(seq_len(nrow(rep_level)), function(i) {
  r <- rep_level[i]

  est_i <- du_estimate(
    gene_muts       = r$gene,
    ref_muts        = r$nongenic,
    gene_len        = r$genespace,
    ref_len         = r$nongenspace,
    nonsyn_muts     = r$Dn,
    syn_muts        = r$Ds,
    p_cds           = p_cds_global,
    neutral_ns_s    = neutral_ns_s_global,
    return_reduction = TRUE
  )

  cbind(r, est_i)
})

rep_est <- rbindlist(rep_est_list, fill = TRUE)

## 3) aggregate over replicates to get mean and SE per parameter combo
##    here we keep only the columns we care to plot
est_summ <- rep_est[
  ,
  .(
    gene_rate_mean      = mean(mu_g_obs) * 1e9,
    gene_rate_se        = sd(mu_g_obs) / sqrt(.N) * 1e9,
    ref_rate_mean       = mean(mu_ref)   * 1e9,
    ref_rate_se         = sd(mu_ref) / sqrt(.N) * 1e9,
    duO_mean            = mean(duO),
    duO_se              = sd(duO) / sqrt(.N),
    duEst_mean          = mean(duEst),
    duEst_se            = sd(duEst) / sqrt(.N),
    dmu_ns_mean         = mean(dmu_ns),
    dmu_ns_se           = sd(dmu_ns) / sqrt(.N),
    dmu_s_mean          = mean(dmu_s),
    dmu_s_se            = sd(dmu_s) / sqrt(.N),
    dmu_nc_mean         = mean(dmu_nc),
    dmu_nc_se           = sd(dmu_nc) / sqrt(.N),
    DnDs_mean           = mean(Dn / pmax(Ds, 1L)),  # crude, but matches your diag
    reps                = .N
  ),
  by = .(S, D, M, S2, S3, D2, D3)
]

## 4) plotting helpers ---------------------------------------------------

plot_estSonly_mean_se <- function(df, ymean, yse, ylab) {
  df2 <- df[S2 == 0 & D %in% c(0, 1)]
  ggplot(df2, aes(x = M, y = .data[[ymean]], color = factor(D), group = D)) +
    geom_line() +
    geom_point(size = 0.6) +
    geom_errorbar(aes(ymin = .data[[ymean]] - .data[[yse]],
                      ymax = .data[[ymean]] + .data[[yse]]),
                  width = 0.02, linewidth = 0.25, alpha = 0.7) +
    facet_grid(~ S) +
    scale_y_continuous(name = ylab) +
    scale_color_manual(values = c("black", "red")) +
    theme_bw(base_size = 6) +
    theme(
      legend.position   = "none",
      strip.background  = element_blank(),
      axis.text.x       = element_text(angle = 90, hjust = 1)
    )
}

## 5) make figures -------------------------------------------------------

pdf("Figures/SLiM_meristem_analyses_reps.pdf", width = 5, height = 4)

lims <- range(c(est_summ$gene_rate_mean, est_summ$ref_rate_mean), na.rm = TRUE)

p1 <- plot_estSonly_mean_se(
  est_summ,
  ymean = "gene_rate_mean",
  yse   = "gene_rate_se",
  ylab  = expression(mu[Genic] %*% 10^-9)
) + scale_y_continuous(limits = lims)

p2 <- plot_estSonly_mean_se(
  est_summ,
  ymean = "ref_rate_mean",
  yse   = "ref_rate_se",
  ylab  = expression(mu[Ref] %*% 10^-9)
) + scale_y_continuous(limits = lims)

p3 <- plot_estSonly_mean_se(
  est_summ,
  ymean = "duO_mean",
  yse   = "duO_se",
  ylab  = expression(Delta * mu[Obs])
)

p4 <- plot_estSonly_mean_se(
  est_summ,
  ymean = "duEst_mean",
  yse   = "duEst_se",
  ylab  = expression(Delta * mu[Est])
)

cowplot::plot_grid(
  p1, p2, p3, p4,
  ncol = 1,
  rel_heights = c(1.1, 1, 1, 1.2)
)

dev.off()

pdf("Figures/SLiM_meristem_analyses_du_components_reps.pdf", width = 3.5, height = 2)

pc1 <- plot_estSonly_mean_se(
  est_summ,
  ymean = "dmu_ns_mean",
  yse   = "dmu_ns_se",
  ylab  = expression(Delta * mu[NS])
)

pc2 <- plot_estSonly_mean_se(
  est_summ,
  ymean = "dmu_nc_mean",
  yse   = "dmu_nc_se",
  ylab  = expression(Delta * mu[NC])
)

pc3 <- plot_estSonly_mean_se(
  est_summ,
  ymean = "dmu_s_mean",
  yse   = "dmu_s_se",
  ylab  = expression(Delta * mu[S])
)

cowplot::plot_grid(pc1, pc2, pc3, ncol = 1, rel_heights = c(1.1, 1, 1.2))

dev.off()
