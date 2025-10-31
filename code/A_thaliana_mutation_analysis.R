source("code/functions.R")
source("code/load_MA_data.R")

# selection/bias across sources --------------------------

p_mutsrc <- plot_mutsrc_boot(mutations, B = 1000)

pdf("Figures/mutation_source_At.pdf", width=3.5, height=3.5)
p_mutsrc
dev.off()

# trimer contexts ---------------------------------------------------------


## usage
tri_plots <- make_genic_intergenic_trimer_plots(mutations, genes_triN, intergenic_triN)

pdf("Figures/Trimer_contet_genic_intergenic.pdf", width=7, height=5)
plot_grid2(tri_plots, type="rows")
dev.off()

pdf("Figures/Trimer_contet_sources.pdf", width=7, height=10)
tri_plots <- make_src_trimer_plots(mutations, genome_triN)
plot_grid2(tri_plots, type = "rows")
dev.off()

# Mono contexts -----------------------------------------------------------


mono_genic<-plot_tricontexts(contexts = mutations[gene_body==T]$trimer, full = F, trimer_freq = genes_triN)
mono_nongenic<-plot_tricontexts(mutations[gene_body==F]$trimer, full = F, trimer_freq = intergenic_triN)
mono_nongenic$context_table$genic<-mono_genic$context_table$N

plot_mono_nongenic<-mono_nongenic$plot+scale_y_continuous(limits=c(0, 2e-4))+ggtitle("Intergenic")
plot_mono_genic<-mono_genic$plot+scale_y_continuous(limits=c(0, 2e-4))+ggtitle("Genic")

plot_grid2(list(plot_mono_nongenic, plot_mono_genic))



# Region comparisons and pie charts ---------------------------------------


# gene body bp per chromosome
gene_pos_by_chr <- genes[, {
  pos <- unlist(lapply(seq_len(.N), function(i) START[i]:STOP[i]))
  list(unique_positions = list(unique(pos)))
}, by = CHROM]$unique_positions
names(gene_pos_by_chr) <- genes$CHROM[!duplicated(genes$CHROM)]
total_gene_length_genome <- sum(vapply(gene_pos_by_chr, length, integer(1L)))

# CDS bp per chromosome
CDS_pos_by_chr <- CDS[, {
  pos <- unlist(lapply(seq_len(.N), function(i) START[i]:STOP[i]))
  list(unique_positions = list(unique(pos)))
}, by = CHROM]$unique_positions
names(CDS_pos_by_chr) <- CDS$CHROM[!duplicated(CDS$CHROM)]
total_CDS_length_genome <- sum(vapply(CDS_pos_by_chr, length, integer(1L)))

# region sizes
non_coding_gene_size <- total_gene_length_genome - total_CDS_length_genome
genome_size          <- sum(vapply(TAIR10, length, integer(1L)))
intergenic_size      <- genome_size - total_gene_length_genome

# mutation counts per region
CDS_muts             <- sum(mutations$CDS == TRUE, na.rm = TRUE)
noncoding_genic_muts <- sum(mutations$gene_body == TRUE & mutations$CDS == FALSE, na.rm = TRUE)
nongenic_muts        <- sum(mutations$gene_body == FALSE, na.rm = TRUE)

pie_dt <- data.table(
  Region = factor(c("CDS", "NonCoding", "Intergenic"),
                  levels = c("CDS", "NonCoding", "Intergenic")),
  L = c(total_CDS_length_genome,
        non_coding_gene_size,
        intergenic_size),
  M = c(CDS_muts,
        noncoding_genic_muts,
        nongenic_muts)
)

# enrichment tests
chisq.test(pie_dt[c(1, 3), .(L, M)])  # CDS vs intergenic
chisq.test(pie_dt[c(2, 3), .(L, M)])  # non-CDS vs intergenic
chisq.test(pie_dt[c(1, 2), .(L, M)])  # CDS vs non-CDS

# pie: genome space
p_genome_space <- ggplot(pie_dt, aes(x = 1, y = L, fill = Region)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.25) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = c(
    CDS        = "green3",
    NonCoding  = "green2",
    Intergenic = "gray70"
  )) +
  theme_void(base_size = 6) +
  theme(
    legend.key.size  = unit(0.25, "cm"),
    panel.background = element_blank(),
    plot.background  = element_blank()
  )

# pie: mutation counts
p_genome_muts <- ggplot(pie_dt, aes(x = 1, y = M, fill = Region)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.25) +
  coord_polar(theta = "y") +
  scale_fill_manual(values = c(
    CDS        = "green3",
    NonCoding  = "green2",
    Intergenic = "gray70"
  )) +
  theme_void(base_size = 6) +
  theme(
    legend.position  = "none",
    panel.background = element_blank(),
    plot.background  = element_blank()
  )

pdf("Figures/Arab_piecharts.pdf", height = 1.5, width = 1.75)
print(p_genome_space)
dev.off()

pdf("Figures/Arab_piecharts_M.pdf", height = 1, width = 1)
print(p_genome_muts)
dev.off()


# Mutation rates in gene types ---------------------------------------------
genic_mut<-mutations[gene_body==T]

gene_annt$mutations<-sites_in_features(gene_annt, mutations, mode = "counts")$counts

essential<-gene_annt[,.(length=sum(STOP-START), mutations=sum(mutations)), by=.(ESN=essentiality=="ESN")]
essential$rate<-essential$mutations/essential$length
chisq.test(essential[,2:3])


lethal<-gene_annt[,.(length=sum(STOP-START), mutations=sum(mutations)), by=.(ESN=lethal=="Lethal")]
lethal$rate<-lethal$mutations/lethal$length
chisq.test(lethal[,2:3])

enrich_H3K4me1<-gene_annt[,.(length=sum(STOP-START), mutations=sum(mutations)), by=.(ESN=enrich_H3K4me1>1)]
enrich_H3K4me1$rate<-enrich_H3K4me1$mutations/enrich_H3K4me1$length
chisq.test(enrich_H3K4me1[,2:3])

# plot Ns/S components points ---------------------------------------------

pdf("Figures/Arab_mutationeffects.pdf", width=2.5, height=1.75)
ggplot(CDS_dt_all_prob_means, aes(x=REF, y=ALT, fill=log(NSS), size=log10(prob)))+
  geom_point(shape=21)+
  scale_fill_gradient2(midpoint = 1, name="log(Ns/S)", low=c("blue"), high=c("orange"),
                       guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"))+
  theme_bw(base_size = 6)+
  scale_size_continuous(range = c(1,10), name="log10(prob)", breaks = c(-5,-4.5, -4))+
  theme(legend.key.size = unit(0.25,"cm"),
        panel.background = element_blank(),
        plot.background = element_blank(),
        legend.box = "vertical")+
  coord_fixed()
dev.off()


set.seed(1)


# Neutral vs observed Ns/S ------------------------------------------------


# simulate neutral effects once, using genome-wide codon/mutation spectrum
effects <- sample(
  CDS_dt_all_prob_means_split$effect,
  size    = 1e6,
  replace = TRUE,
  prob    = prop.table(CDS_dt_all_prob_means_split$totalprob)
)
# prop.table(table(effects))  # optional check


# observed sets
all_coding  <- mutations[mutations$mutations_effect %in% c("Non-Syn", "Syn")]
genic_mut   <- mutations[gene_body == TRUE]
esn_sub     <- genic_mut[essential == TRUE]
lethal_sub  <- genic_mut[lethal == TRUE]
h3k4_sub    <- genic_mut[enrich_H3K4me1 == TRUE]
hk_sub      <- genic_mut[tissue_genes == TRUE]

# observed counts
all_counts    <- get_ns_s(all_coding)
esn_counts    <- get_ns_s(esn_sub)
lethal_counts <- get_ns_s(lethal_sub)
h3k4_counts   <- get_ns_s(h3k4_sub)
hk_counts     <- get_ns_s(hk_sub)

# sample from neutral model; each returns a vector of NS/S ratios
n_iter <- 10000
all_neutrals     <- sample_neutrals(effects, all_counts$total,    n_iter)
ESN_neutrals     <- sample_neutrals(effects, esn_counts$total,    n_iter)
lethal_neutrals  <- sample_neutrals(effects, lethal_counts$total, n_iter)
h3k4_neutrals    <- sample_neutrals(effects, h3k4_counts$total,   n_iter)
hk_neutrals      <- sample_neutrals(effects, hk_counts$total,     n_iter)

# stack neutral draws
merge_neutrals <- data.table(
  s = c(
    all_neutrals,
    h3k4_neutrals,
    lethal_neutrals,
    ESN_neutrals,
    hk_neutrals
  ),
  Genes = rep(
    c("All", "H3K4me1", "'Lethal'", "'Essential'", "'Housekeeping'"),
    each = n_iter
  )
)

merge_neutrals$Genes <- factor(
  merge_neutrals$Genes,
  levels = c("All", "H3K4me1", "'Lethal'", "'Essential'", "'Housekeeping'")
)

# observed Dn/Ds per set
obs_DnDs <- rbindlist(list(
  data.table(Genes = "All",
             NS = all_counts$NS,    S = all_counts$S),
  data.table(Genes = "'Essential'",
             NS = esn_counts$NS,    S = esn_counts$S),
  data.table(Genes = "'Lethal'",
             NS = lethal_counts$NS, S = lethal_counts$S),
  data.table(Genes = "H3K4me1",
             NS = h3k4_counts$NS,   S = h3k4_counts$S),
  data.table(Genes = "'Housekeeping'",
             NS = hk_counts$NS,     S = hk_counts$S)
))

obs_DnDs[, s := NS / S]
obs_DnDs[, N := NS + S]

# empirical percentile: how extreme is observed vs neutral draws
obs_DnDs[, pctile := sapply(seq_len(.N), function(i) {
  g   <- Genes[i]
  s0  <- s[i]
  sub <- merge_neutrals[Genes == g]
  1 - sum(sub$s > s0) / nrow(sub)
})]

# melt for table figure
obs_DnDs_melt <- melt(obs_DnDs[, .(Genes, NS, S, s)], id.vars = "Genes" )

obs_DnDs_melt$Genes <- factor(
  obs_DnDs_melt$Genes,
  levels = levels(merge_neutrals$Genes)
)
obs_DnDs_melt$variable <- factor(
  obs_DnDs_melt$variable,
  levels = rev(c("NS", "S", "s"))
)

pdf("Figures/DnDs_Arabidopsis_table.pdf", width = 2, height = 0.5)
ggplot(obs_DnDs_melt, aes(x = Genes, y = variable, label = round(value, 2))) +
  geom_tile(fill = "white", col = "black") +
  geom_text(size = 2) +
  scale_y_discrete(labels = c("NS/S", "S", "NS")) +
  theme_classic(base_size = 6) +
  theme(
    axis.line   = element_blank(),
    axis.ticks  = element_blank(),
    axis.text.x = element_blank(),
    axis.title  = element_blank()
  )
dev.off()

pdf("Figures/DnDs_Arabidopsis.pdf", width = 2, height = 1.5)
ggplot(merge_neutrals, aes(x = Genes, y = s)) +
  geom_jitter(height = 0, size = 0.1, alpha = 0.02, col = "gray", width = 0.25) +
  scale_y_log10(name = "Ns/S") +
  geom_point(
    data  = obs_DnDs,
    aes(x = Genes, y = s),
    fill  = "green3",
    size  = 1,
    col   = "black",
    shape = 21
  ) +
  theme_classic(base_size = 6) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()



# mutation delta  -------------------------------------

ESN_Du<-calculate_dUest(loci_muts=mutations[essential==T],
                        other_muts = mutations[gene_body==T & essential==F],
                        loci_genes = gene_annt[essentiality=="ESN"],
                        loci_cds = CDS[gene %in%  gene_annt[essentiality=="ESN"]$gene],
                        genes=gene_annt
)

lethal_Du<-calculate_dUest(loci_muts=mutations[lethal==T],
                           other_muts = mutations[gene_body==T & lethal==F],
                           loci_genes = gene_annt[lethal=="Lethal"],
                           loci_cds = CDS[gene %in%  gene_annt[lethal=="Lethal"]$gene],
                           genes=gene_annt
)

H3K4me1_Du<-calculate_dUest(loci_muts=mutations[enrich_H3K4me1==T & gene_body==T],
                            other_muts = mutations[gene_body==T & enrich_H3K4me1==F],
                            loci_genes = gene_annt[enrich_H3K4me1>0],
                            loci_cds = CDS[gene %in%  gene_annt[enrich_H3K4me1>0]$gene],
                            genes=gene_annt
)

Tissue_Du<-calculate_dUest(loci_muts=mutations[tissue_genes==T & gene_body==T],
                           other_muts = mutations[gene_body==T & tissue_genes==F],
                           loci_genes = gene_annt[tissue_breadth==54],
                           loci_cds = CDS[gene %in%  gene_annt[tissue_breadth==54]$gene],
                           genes=gene_annt
)



genome_size<-sum(unlist(lapply(TAIR10, length)))
all_Du<-calculate_dUest(loci_muts=mutations[gene_body==T],
                        other_muts = mutations[gene_body==F],
                        loci_genes = gene_annt,
                        loci_cds = CDS[gene %in%  gene_annt$gene],
                        other_L = genome_size-total_gene_length_genome,
                        DnDsN=mean(all_neutrals)
)

H3K4me1_IG<-calculate_dUest(loci_muts=mutations[enrich_H3K4me1==T & gene_body==T],
                            other_muts = mutations[gene_body==F ],
                            loci_genes = gene_annt[enrich_H3K4me1>0],
                            loci_cds = CDS[gene %in%  gene_annt[enrich_H3K4me1>0]$gene],
                            other_L = genome_size-total_gene_length_genome,
                            DnDsN=mean(all_neutrals)
)

lethal_IG<-calculate_dUest(loci_muts=mutations[lethal==T],
                           other_muts = mutations[gene_body==F],
                           loci_genes = gene_annt[lethal=="Lethal"],
                           loci_cds = CDS[gene %in%  gene_annt[lethal=="Lethal"]$gene],
                           other_L = genome_size-total_gene_length_genome,
                           DnDsN=mean(all_neutrals)
)

ESN_IG<-calculate_dUest(loci_muts=mutations[essential==T],
                        other_muts = mutations[gene_body==F],
                        loci_genes = gene_annt[essentiality=="ESN"],
                        loci_cds = CDS[gene %in%  gene_annt[essentiality=="ESN"]$gene],
                        other_L = genome_size-total_gene_length_genome,
                        DnDsN=mean(all_neutrals)
)

Tissue_IG<-calculate_dUest(loci_muts=mutations[tissue_genes==T],
                           other_muts = mutations[gene_body==F],
                           loci_genes = gene_annt[tissue_breadth==54],
                           loci_cds = CDS[gene %in%  gene_annt[tissue_breadth==54]$gene],
                           other_L = genome_size-total_gene_length_genome,
                           DnDsN=mean(all_neutrals)
)




merge_Du<-rbindlist(list(all_Du, H3K4me1_Du, lethal_Du, Tissue_Du, ESN_Du, H3K4me1_IG, lethal_IG, ESN_IG, Tissue_IG))
merge_Du$Genes<-factor(c("All_IG","H3K4me1","Lethal","Housekeeping","Essential","H3K4me1_IG","lethal_IG", "ESN_IG","Housekeeping_IG"))
merge_Du$GvsIG<-factor(ifelse(grepl("IG", merge_Du$Genes),"Genic vs Intergenic","Genic vs Genic"), levels=c("Genic vs Intergenic","Genic vs Genic"))
merge_Du$chi<-sapply(1:nrow(merge_Du), function(i){
  loci_mutN=merge_Du$loci_mutN[i]
  other_mutN=merge_Du$other_mutN[i]
  LG=merge_Du$LG[i]
  other_L=merge_Du$other_L[i]
  chi<-chisq.test(matrix(c(loci_mutN, other_mutN, LG, other_L), nrow=2))
  chi$p.value
})

labels<-c(All="All genes", IG="Intergenic", All_IG="All genes",H3K4me1="H3K4me1 enriched genes",Lethal="'Lethal genes'",Essential="'Essential' genes",H3K4me1_IG="H3K4me1 enriched genes",lethal_IG="'Lethal genes'",ESN_IG="'Essential' genes", Housekeeping_IG="'Housekeeping' genes", Housekeeping="'Housekeeping' genes")

pdf("Figures/du_arabidopsis.pdf", width=7, height=2.5)

merge_Du_melt<-melt(merge_Du[,.(GvsIG, Genes, duO, dµNS, dµS, dµNC, duEst)], id.vars = c("Genes","GvsIG"))

dodge <- position_dodge(width = 0.75)

p1 <- ggplot(merge_Du_melt, aes(x = Genes, fill = variable)) +
  # draw from 1 down to (1 - value)
  geom_linerange(
    aes(ymin = 1 - value, ymax = 1, group = variable, color = variable),
    position = dodge,
    linewidth = 1.5
  ) +
  geom_linerange(
    aes(ymin = 1 - value, ymax = 1, group = variable, color = variable),
    position = dodge,
    linewidth = 1.5
  ) +
  scale_y_continuous(name = expression(Delta * mu), limits = c(0, 1), expand=c(0,0)) +
  facet_grid(~GvsIG, scales="free", space="free")+
  scale_color_viridis_d(
    labels = c(
      expression(Delta * mu[Obs]),
      expression(Delta * mu[NS]),
      expression(Delta * mu[S]),
      expression(Delta * mu[NC]),
      expression(Delta * mu[Est])
    ),
    name = ""
  ) +
  scale_fill_viridis_d(guide = "none") +
  scale_x_discrete(labels = labels) +
  theme_classic(base_size = 6) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.key.size = unit(0.25, "cm"),
    strip.background = element_blank(),
    plot.background  = element_rect(),
    axis.line.x.top = element_line(),
    axis.line.x.bottom = element_blank(),
    panel.border =   element_rect(linewidth = 0.5)
  )+
  geom_hline(yintercept = 1, linewidth=0.5)

merge_Du_melt2<-melt(merge_Du[,.(GvsIG, Genes, u_O, µI_O)], id.vars = c("Genes","GvsIG"))
IG_rate<-data.table(variable="IG", Genes="Intergenic", value=unique(merge_Du_melt2[GvsIG=="Genic vs Intergenic" & variable=="µI_O"]$value))
all_rate<-data.table(variable="All", Genes="All genes", value=unique(merge_Du_melt2[Genes=="All_IG" & GvsIG=="Genic vs Intergenic" & variable=="u_O"]$value))
merge_Du_melt2<-rbindlist(list( IG_rate, all_rate, merge_Du_melt2[GvsIG!="Genic vs Intergenic"]), fill=T)

merge_Du_melt2$variable<-factor(merge_Du_melt2$variable, levels=c("All", "u_O","µI_O","IG"))

p2 <- ggplot(
  merge_Du_melt2[is.na(GvsIG) | variable == "u_O"],
  aes(x = Genes,
      fill = ifelse(variable == "IG", "IG", "G"),
      y = value)
) +
  geom_bar(stat = "identity", position = "dodge",
           col = "black", linewidth = 0.25, width = 0.5) +
  scale_y_continuous(name = "µ (mutations/b.p.)", breaks=c(0, 2.5e-5, 5e-5, 7.5e-5), labels=real_sci_format(c(0, 2.5e-5, 5e-5, 7.5e-5))) +
  theme_classic(base_size = 6) +
  scale_fill_manual(
    values = c(G = "green4", IG = "gray70"),
    labels = c("Genes", "Intergenic"),
    name   = "",
    guide  = guide_legend(nrow = 2, byrow = TRUE)
  ) +
  scale_x_discrete(labels = labels) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1),
    legend.position  = "top",
    legend.key.size  = unit(0.25, "cm"),
    strip.background = element_blank(),
    panel.background = element_blank()
  )


p2.2 <- ggplot(
  merge_Du_melt2[!is.na(GvsIG)],
  aes(x = Genes, fill = variable, y = value)
) +
  geom_bar(stat = "identity", position = "dodge",
           col = "black", linewidth = 0.25, width = 0.5) +
  scale_y_continuous(name = "µ (mutations/b.p.)",breaks=c(0, 2e-5, 4e-5, 6.5e-5), labels=real_sci_format(c(0, 2e-5, 4e-5, 6e-5))) +
  theme_classic(base_size = 6) +
  scale_fill_manual(
    values = c(All = "green4", u_O = "green4", µI_O = "gray70", IG = "gray70"),
    labels = c("Test genes", "Other genes", "Other genes", "Intergenic"),
    name   = "",
    guide  = guide_legend(nrow = 4, byrow = TRUE)
  ) +
  scale_x_discrete(labels = labels) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1),
    legend.position  = "top",
    legend.key.size  = unit(0.25, "cm"),
    strip.background = element_blank(),
    panel.background = element_blank()
  )

merge_Du_melt3<-melt(merge_Du[GvsIG=="Genic vs Genic",.(Genes,DnDsO, DnDsN)], id.vars = c("Genes"))

cowplot::plot_grid(plotlist = list(p2,p2.2, p1), ncol=3, rel_widths = c(.5,.6,2))
dev.off()


# H3K4me1 -----------------------------------------------------------------

pdf("Figures/Arab_H3K4me1_genes.pdf", width=2, height=2)

gene_annt_mean<-gene_annt[,.(mut=sum(germline_mutations_meta), LENGTH=sum(STOP-START), rate=sum(germline_mutations_meta)/sum(STOP-START)), by=.(grp=as.numeric(cut2(enrich_H3K4me1, g=10)))]

ggplot(gene_annt_mean, aes(x=factor(grp*20), y=rate, fill=grp*20))+
  geom_bar(stat="identity", col="black", width=0.55)+
  scale_fill_gradient(low="gray90", high="green3", name="H3K4me1 %ile",
                      guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"))+
  theme_classic(base_size = 6)+
  scale_y_continuous(name="Mutations/bp")+
  scale_x_discrete(name="H3K4me1 %ile")+
  theme(legend.key.size = unit(0.25,"cm"),
        panel.background = element_blank(),
        plot.background = element_blank(),
        legend.box = "vertical")
dev.off()

# archive -----------------------------------------------------------------

pdf("Figures/Arab_NsS_genes_comparison.pdf", width=2, height=2)
p3<-ggplot(merge_Du_melt3, aes(x=Genes, fill=variable, y=(value)))+
  geom_point(shape=21, col="black", size=1)+
  scale_y_continuous(name="Ns/S")+
  theme_classic(base_size = 6)+
  scale_fill_manual(values=c(DnDsO="green4",DnDsN="gray70"), labels=c("Test genes","Other genes"), name="") +
  scale_x_discrete(labels=labels)+
  theme(axis.text.x = element_text(angle=45, hjust=1),legend.position = "top",
        legend.key.size = unit(0.25,"cm"),strip.background = element_blank(), panel.background = element_blank()
  )
p3
dev.off()



neutral_groups<-rbindlist(lapply(c("essential","lethal","enrich_H3K4me1","all","tissue_genes"),function(g){
  message(g)
  sample_neutrals_groups(g, reps=10000)
}))

pdf("Figures/DnDs_Arabidopsis.pdf", width=1.5, height=2)
dodge <- position_dodge(width = 0.75)
ggplot(neutral_groups, aes(x=group, y=(s), col=inout))+
  #geom_jitter(height=0, size=0.1, alpha=0.01, col="orange")+
  geom_violin(fill="gray90", alpha=0.1, size=0.25, position = dodge)+
  scale_y_log10(name="Ns/S")+
  geom_point(aes(y=NsSO), size=0.5,position = dodge)+
  theme_classic(base_size = 6)+
  theme(axis.text.x = element_text(angle=45, hjust=1))
dev.off()



