source("code/functions.R")

## ------------------------------------------------------------------------
## data loading
## ------------------------------------------------------------------------

# 1) coding sequences and genome FASTAs
cds    <- read.fasta("data/TAIR10/TAIR10_cds_20110103_representative_gene_model_updated.fa.gz")
TAIR10 <- read.fasta("data/TAIR10/TAIR10_chr_all.fa.gz")
names(TAIR10) <- gsub("Chr", "", names(TAIR10))

# 2) gene annotation
gene_annt <- fread("data/TAIR10/arabidopsis_genes.csv")

# 3) GFF3 (genes + CDS)
gff <- read.GFF("data/TAIR10/TAIR10_GFF3_genes.gff")
gff$CHROM <- gsub("Chr", "", gff$CHROM)

# CDS features only
CDS <- gff[TYPE == "CDS"]
CDS$gene <- gsub("\\..+", "", gsub("Parent=(.+),.+", "\\1", CDS$INFO))

# normalize model names to CDS fasta names
CDS$model <- gsub("Parent=(.+),.+", "\\1", CDS$INFO)
CDS$model <- gsub("-Protein", "", CDS$model)
CDS <- CDS[model %in% names(cds)]
CDS$ID <- seq_len(nrow(CDS))

# gene features
genes <- gff[TYPE == "gene"]
genes$gene <- gsub(".+Name=", "", genes$INFO)
genes <- genes[gene %in% CDS$gene]
genes <- genes[CHROM %in% 1:5]
genes$ID <- seq_len(nrow(genes))

## ------------------------------------------------------------------------
## define intergenic regions (between genes on same chromosome)
## ------------------------------------------------------------------------
genes <- genes[order(CHROM, START)]
intergenic <- genes[, .(
  SOURCE = "TAIR10",
  TYPE   = "intergenic",
  START  = STOP + 1L,
  STOP   = data.table::shift(START, type = "lead") - 1L
), by = CHROM]

intergenic <- intergenic[!is.na(STOP)]
intergenic <- intergenic[START <= STOP]
intergenic$ID <- seq_len(nrow(intergenic))

## ------------------------------------------------------------------------
## mutations table
## ------------------------------------------------------------------------
mutations <- fread("data/germline_mutations_meta.csv")
mutations$ID <- seq_len(nrow(mutations))

## ------------------------------------------------------------------------
## annotate mutations by feature sets
## ------------------------------------------------------------------------
mutations$CDS           <- annot_mut_geneset(mutations, CDS)
mutations$intergenic    <- annot_mut_geneset(mutations, intergenic)
mutations$gene_body     <- annot_mut_geneset(mutations, genes)
mutations$essential     <- annot_mut_geneset(mutations, gene_annt[essentiality == "ESN"]) #Lloyd et al. 2011
mutations$enrich_H3K4me1 <- annot_mut_geneset(mutations, gene_annt[enrich_H3K4me1 > 0]) #H3k4me1 enrichment Niu et al. 2021
mutations$tissue_genes  <- annot_mut_geneset(mutations, gene_annt[tissue_breadth == 54]) #Mergner et al. 2016
mutations$lethal        <- annot_mut_geneset(mutations, gene_annt[lethal == "Lethal"]) #Lloyed et al. 2015

# tri-nucleotide context
mutations$trimer <- tricontexts(mutations, TAIR10)

# reassign CDS IDs (already done above; keep consistent)
CDS$ID <- seq_len(nrow(CDS))

mutations <- mutations[order(ID)]

## ------------------------------------------------------------------------
## per-mutation coding effect
## ------------------------------------------------------------------------
mutation_effects <- rbindlist(
  pblapply(mutations$ID, function(i) {
    annot_mut_effect(
      mut_id       = i,
      mutations    = mutations,
      CDS          = CDS,
      cds          = cds,
      genome_fasta = TAIR10
    )
  })
)

# quick summary
print(table(mutation_effects$effect))


## ------------------------------------------------------------------------
## merge mutation metadata with effects and save
## ------------------------------------------------------------------------
mutation_annotated <- merge(mutations, mutation_effects, by = c("ID", "CHROM","POS"))

#QC check
table(mutation_annotated$CHROM, substr(mutation_annotated$GENE_model,1,3))

print(table(mutation_annotated$src))
# keep original effect in a dedicated column
mutation_annotated$mutations_effect <- mutation_annotated$effect
mutation_annotated[, effect := NULL]

fwrite(mutation_annotated, "Tables/mutations_Athal_annotated.csv")




# other data prep ---------------------------------------------------------
## genic trimers
genes_tri  <- Nmerfrequency(genes, TAIR10, Nmer = 3, mode = "counts")
genes_triN <- lapply(genes_tri, sum)
genes_triN <- data.table(TRI = names(genes_triN), N = unlist(genes_triN))[TRI != "ID"]
# keep only A/T/C/G trimers
genes_triN <- genes_triN[grepl("^[ACGT]{3}$", TRI)]
fwrite(genes_triN, "data/TAIR10/genes_trimers.csv")

## intergenic trimers
intergenic_tri  <- Nmerfrequency(intergenic, TAIR10, Nmer = 3, mode = "counts")
intergenic_triN <- lapply(intergenic_tri, sum)
intergenic_triN <- data.table(TRI = names(intergenic_triN), N = unlist(intergenic_triN))[TRI != "ID"]
intergenic_triN <- intergenic_triN[grepl("^[ACGT]{3}$", TRI)]
fwrite(intergenic_triN, "data/TAIR10/intergenic_trimers.csv")

## whole-genome trimers (all TAIR10 chromosomes loaded)
# combine genic + intergenic features
# full join on TRI, missing → 0
genome_triN <- merge(
  genes_triN, intergenic_triN,
  by = "TRI", all = TRUE,
  suffixes = c("_genes", "_intergenic")
)

genome_triN[is.na(N_genes), N_genes := 0L]
genome_triN[is.na(N_intergenic), N_intergenic := 0L]

# sum
genome_triN[, N := N_genes + N_intergenic]

# keep only A/C/G/T trimers (in case)
genome_triN <- genome_triN[grepl("^[ACGT]{3}$", TRI)]

# final minimalist table
genome_triN <- genome_triN[, .(TRI, N)]

fwrite(genome_triN, "data/TAIR10/genome_trimers.csv")




# build mutation probability+effect tables --------------------------------

## section commented due to size and time required to generate.
# t<-Sys.time()
# CDS_dt_all<-CDS_codon_mutation_table(cds)
# Sys.time()-t
# fwrite(CDS_dt_all, "Tables/TAIR_10_codons.csv")


# #CDS_dt_all<-fread("dir/tables/TAIR_10_codons.csv") # very slow to produce, and huge file
# CDS_dt_all$REF<-as.character(CDS_dt_all$REF)
# mut_probs<-mutation_probs(mutations[CDS==F], TAIR10, cds=cds)
# CDS_dt_all_prob<-merge(CDS_dt_all[,.(REF,ALT,effect)], mut_probs[,.(REF, ALT, rawPROB)], by=c("REF","ALT"))
#
# CDS_dt_all_prob_means<-CDS_dt_all_prob[,.(NSS=sum(effect=="N")/sum(effect=="S"), prob=mean(rawPROB)), by=.(REF=as.character(REF), ALT)]
#
# CDS_dt_all_prob_means_split<-CDS_dt_all_prob[,.(N=.N, prob=mean(rawPROB)), by=.(REF=as.character(REF), ALT, effect)]
# CDS_dt_all_prob_means_split$totalprob<-CDS_dt_all_prob_means_split$N*CDS_dt_all_prob_means_split$prob
#
# fwrite(CDS_dt_all_prob_means, "tables/mutation_probabilities_effects.csv")
# fwrite(CDS_dt_all_prob_means_split, "tables/mutation_probabilities_effects_split.csv")
