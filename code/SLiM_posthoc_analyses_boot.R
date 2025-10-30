library(data.table)
library(ggplot2)
library(cowplot)

source("code/functions.R")  # must contain du_estimate()

## globals
neutral_ns_s_global <- 2
p_cds_global        <- 3/4
B <- 500  # bootstraps per param combo

muts <- fread("data/slim_out/SLiM_mutsM.csv")

## one bootstrap on one param combo --------------------------------------
boot_once_param <- function(dt_param,
                            genespace   = 4000,
                            nongenspace = 4000,
                            neutral_ns_s = neutral_ns_s_global,
                            p_cds        = p_cds_global) {

  n <- nrow(dt_param)
  samp <- dt_param[sample.int(n, n, replace = TRUE)]

  gene_muts <- sum(samp$genic)
  ref_muts  <- sum(!samp$genic)
  Dn        <- sum(samp$Dn)
  Ds        <- sum(samp$Ds)

  du_estimate(
    gene_muts     = gene_muts,
    ref_muts      = ref_muts,
    gene_len      = genespace,
    ref_len       = nongenspace,
    nonsyn_muts   = Dn,
    syn_muts      = Ds,
    p_cds         = p_cds,
    neutral_ns_s  = neutral_ns_s,
    return_reduction = TRUE
  )
}

## outer loop with progress bar ------------------------------------------
param_keys <- unique(muts[, .(S, D, M, S2, S3, D2, D3)])
pb <- txtProgressBar(min = 0, max = nrow(param_keys), style = 3)

boot_results <- lapply(seq_len(nrow(param_keys)), function(i) {
  key <- param_keys[i]

  dt_param <- muts[
    S  == key$S  &
      D  == key$D  &
      M  == key$M  &
      S2 == key$S2 &
      S3 == key$S3 &
      D2 == key$D2 &
      D3 == key$D3
  ]

  if (!nrow(dt_param)) {
    setTxtProgressBar(pb, i)
    return(NULL)
  }

  boots <- rbindlist(
    lapply(seq_len(B), function(b) {
      boot_once_param(dt_param)
    }),
    fill = TRUE
  )

  boots[, `:=`(
    S  = key$S,
    D  = key$D,
    M  = key$M,
    S2 = key$S2,
    S3 = key$S3,
    D2 = key$D2,
    D3 = key$D3
  )]

  setTxtProgressBar(pb, i)
  boots
})

close(pb)

boot_results <- rbindlist(boot_results, fill = TRUE)

## summarise --------------------------------------------------------------
boot_summ <- boot_results[
  ,
  .(
    duO_mean    = mean(duO, na.rm = TRUE),
    duO_lo      = quantile(duO, 0.025, na.rm = TRUE),
    duO_hi      = quantile(duO, 0.975, na.rm = TRUE),

    duEst_mean  = mean(duEst, na.rm = TRUE),
    duEst_lo    = quantile(duEst, 0.025, na.rm = TRUE),
    duEst_hi    = quantile(duEst, 0.975, na.rm = TRUE),

    mu_g_mean   = mean(mu_g_obs, na.rm = TRUE),
    mu_g_lo     = quantile(mu_g_obs, 0.025, na.rm = TRUE),
    mu_g_hi     = quantile(mu_g_obs, 0.975, na.rm = TRUE),

    mu_ref_mean = mean(mu_ref, na.rm = TRUE),
    mu_ref_lo   = quantile(mu_ref, 0.025, na.rm = TRUE),
    mu_ref_hi   = quantile(mu_ref, 0.975, na.rm = TRUE),

    dmu_ns_mean = mean(dmu_ns, na.rm = TRUE),
    dmu_ns_lo   = quantile(dmu_ns, 0.025, na.rm = TRUE),
    dmu_ns_hi   = quantile(dmu_ns, 0.975, na.rm = TRUE),

    dmu_s_mean  = mean(dmu_s, na.rm = TRUE),
    dmu_s_lo    = quantile(dmu_s, 0.025, na.rm = TRUE),
    dmu_s_hi    = quantile(dmu_s, 0.975, na.rm = TRUE),

    dmu_nc_mean = mean(dmu_nc, na.rm = TRUE),
    dmu_nc_lo   = quantile(dmu_nc, 0.025, na.rm = TRUE),
    dmu_nc_hi   = quantile(dmu_nc, 0.975, na.rm = TRUE)
  ),
  by = .(S, D, M, S2, S3, D2, D3)
]

## plotting helper (same as before) --------------------------------------
plot_boot <- function(df, ymean, ylo, yhi, ylab) {
  df2 <- df[S2 == 0 & D %in% c(0, 1)]
  ggplot(df2, aes(x = M, y = .data[[ymean]], color = factor(D), group = D)) +
    geom_line() +
    geom_point(size = 0.6) +
    geom_ribbon(aes(ymin = .data[[ylo]], ymax = .data[[yhi]], fill = factor(D)),
                alpha = 0.15, colour = NA) +
    facet_grid(~ S) +
    scale_y_continuous(name = ylab) +
    scale_color_manual(values = c("black", "red")) +
    scale_fill_manual(values = c("black", "red")) +
    theme_bw(base_size = 6) +
    theme(
      legend.position  = "none",
      strip.background = element_blank(),
      axis.text.x      = element_text(angle = 90, hjust = 1)
    )
}


## make the plots ---------------------------------------------------

pdf("Figures/SLiM_meristem_analyses_boot.pdf", width = 5, height = 4)

p1 <- plot_boot(
  boot_summ,
  ymean = "mu_g_mean",
  ylo   = "mu_g_lo",
  yhi   = "mu_g_hi",
  ylab  = expression(mu[Genic])
)

p2 <- plot_boot(
  boot_summ,
  ymean = "mu_ref_mean",
  ylo   = "mu_ref_lo",
  yhi   = "mu_ref_hi",
  ylab  = expression(mu[Ref])
)

p3 <- plot_boot(
  boot_summ,
  ymean = "duO_mean",
  ylo   = "duO_lo",
  yhi   = "duO_hi",
  ylab  = expression(Delta * mu[Obs])
)

p4 <- plot_boot(
  boot_summ,
  ymean = "duEst_mean",
  ylo   = "duEst_lo",
  yhi   = "duEst_hi",
  ylab  = expression(Delta * mu[Est])
)

cowplot::plot_grid(p1, p2, p3, p4, ncol = 1, rel_heights = c(1.1, 1, 1, 1.2))

dev.off()

## components plot
pdf("Figures/SLiM_meristem_analyses_du_components_boot.pdf", width = 3.5, height = 2)

pc1 <- plot_boot(
  boot_summ,
  ymean = "dmu_ns_mean",
  ylo   = "dmu_ns_lo",
  yhi   = "dmu_ns_hi",
  ylab  = expression(Delta * mu[NS])
)

pc2 <- plot_boot(
  boot_summ,
  ymean = "dmu_nc_mean",
  ylo   = "dmu_nc_lo",
  yhi   = "dmu_nc_hi",
  ylab  = expression(Delta * mu[NC])
)

pc3 <- plot_boot(
  boot_summ,
  ymean = "dmu_s_mean",
  ylo   = "dmu_s_lo",
  yhi   = "dmu_s_hi",
  ylab  = expression(Delta * mu[S])
)

cowplot::plot_grid(pc1, pc2, pc3, ncol = 1, rel_heights = c(1.1, 1, 1.2))

dev.off()
