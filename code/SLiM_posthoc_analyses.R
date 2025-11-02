
source("code/functions.R")
## 1. Constants / global params
## ------------------------------------------------------------
neutral_ns_s_global <- 2     # true simulated neutral NS/S
p_cds_global        <- 3/4   # true SLiM genome layout
pns     <- p_cds_global * neutral_ns_s_global / (neutral_ns_s_global + 1)  # expected nonsyn fraction
ps      <- p_cds_global - pns                     # expected syn fraction

##  Read SLiM output and aggregate
## ------------------------------------------------------------
muts <- fread("data/slim_out/SLiM_mutsM.csv.gz")

# collapse to parameter-level means

all_means <- muts[
  ,
  .(
    genespace    = 4000,
    nongenspace  = 4000,
    gene         = sum(genic),
    nongenic     = sum(!genic),
    DnDs         = sum(Dn)/sum(Ds),
    Dn           = sum(Dn),
    Ds           = sum(Ds)
  ),
  by = .(S, D, M, S2, S3, D2, D3)
]

## Apply mutation rate estimator to every param combo
## ------------------------------------------------------------

est_list <- lapply(seq_len(nrow(all_means)), function(i) {
  r <- all_means[i]

  est_i <- du_estimate(
    gene_muts     = r$gene,
    ref_muts      = r$nongenic,
    gene_len      = r$genespace,
    ref_len       = r$nongenspace,
    nonsyn_muts   = r$Dn,
    syn_muts      = r$Ds,
    p_cds         = p_cds_global,
    neutral_ns_s  = neutral_ns_s_global,
    return_reduction = TRUE
  )

  cbind(r, est_i)
})

est <- rbindlist(est_list, fill = TRUE)

## (optional) a quick chi on Dn/Ds against 2:1 like you had
est[, chi_DnDs := sapply(seq_len(.N), function(i) {
  chisq.test(c(est$Dn[i], est$Ds[i]), p = prop.table(c(2, 1)))$p.value
})]
est[, DnDs_p := -log10(chi_DnDs)]


pdf("Figures/SLiM_meristem_analyses.pdf", width = 5, height = 4)

## 1) use the new rate columns from du_estimate ------------------------
## est$mu_g_obs = genic muts / genic bp
## est$mu_ref   = ref muts   / ref bp
est$gene_rate  <- est$mu_g_obs /4000/3/15000/10^-9     # mutations/bp/generation (~rescaled to At genome)
est$ref_rate   <- est$mu_ref  /4000/3/15000/10^-9

## subset to target parameters
est2 <- est[S2 == 0 & D %in% c(0, 1)]

fwrite(est2,"tables/SLiM_out.csv")

## y-range shared by first two plots
lims <- range(c(est2$ref_rate, est2$gene_rate), na.rm = TRUE)

## 2) genic mutation rate vs M, faceted by S (like before)
p1 <- plot_estSonly(est, "gene_rate", "Genic mutations") +
  scale_y_continuous(
    limits = lims,
    name   = expression(mu[Genic] %*% 10^-9)
  ) +
  theme(strip.text = element_text())

## 3) reference (non-genic) mutation rate vs M
p2 <- plot_estSonly(est, "ref_rate", "Intergenic / reference mutations") +
  scale_y_continuous(
    limits = lims,
    name   = expression(mu[Ref] %*% 10^-9)
  )

## 4) observed genic reduction (this is 1 - duO)
p3 <- plot_estSonly(est, "duO", "Δµ (observed)") +
  scale_y_continuous(name = expression(Delta * mu[Obs]))

## 5) Ns/S diagnostic (same as before)
p4 <- plot_estSonly(est, "DnDs", "Ns/S") +
  geom_hline(yintercept = 2, linetype = "dashed", alpha = 0.5, size = 0.5) +
  scale_y_log10(name = "Ns/S")

## 6) -log10 p for Ns/S (if you kept that calc)
p5 <- plot_estSonly(est, "DnDs_p", "-log10(p) Ns/S") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5) +
  scale_y_continuous(name = expression(log[10](p) * " Ns/S"))

p6 <- plot_estSonly(est, "duEst", "Δµ (Est)") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", alpha = 0.5, size = 0.5) +
  scale_x_continuous(name = expression(Delta * mu[Genic])) +  # purely cosmetic; x is still M inside plot_estSonly
  scale_y_continuous(name = expression(Delta * mu[Est])) +
  theme(axis.text.x = element_text(), axis.title.x = element_text())

plot_grid2(
  plotlist     = list(p1, p2, p3, p4, p5, p6),
  type         = "rows",
  relsize  = c(1.1, 1, 1, 1, 1, 1.1)
)
dev.off()

pdf("Figures/SLiM_meristem_analyses_du_components.pdf", width = 3.5, height = 2)

est_sim <- est[S > 0.03]

pc1 <- plot_estSonly(est_sim, "dmu_ns", expression(Delta * mu[NS])) +
  theme(strip.text = element_text())

pc2 <- plot_estSonly(est_sim, "dmu_nc", expression(Delta * mu[NC]))

pc3 <- plot_estSonly(est_sim, "dmu_s", expression(Delta * mu[S])) +
  scale_x_continuous(name = expression(Delta * mu[Obs])) +  # x is M in plot_estSonly, but label is just cosmetic
  theme(axis.text.x = element_text(), axis.title.x = element_text())

plot_grid2(
  plotlist     = list(pc1, pc2, pc3),
  type         = "rows",
  relsize  = c(1.1, 1, 1)
)


dev.off()
