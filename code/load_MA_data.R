## ------------------------------------------------------------
## load_MA_data.R
## Load TAIR10 references, features, intergenic regions,
## precomputed 3-mer counts, and annotated mutation table
## ------------------------------------------------------------

message("[load_MA_data] Loading TAIR10 CDS and genome FASTA...")

## 1) TAIR10 sequences ----------------------------------------------------

# representative CDS sequences (TAIR10)
cds <- read.fasta("data/TAIR10/TAIR10_cds_20110103_representative_gene_model_updated.fa.gz")

# whole-genome chromosomes (TAIR10)
TAIR10 <- read.fasta("data/TAIR10/TAIR10_chr_all.fa.gz")
names(TAIR10) <- gsub("Chr", "", names(TAIR10))  # "Chr1" -> "1"

message("[load_MA_data] Loading gene annotations and GFF3...")

## 2) gene-level annotation tables (literature sets, essential, etc.) -----

gene_annt <- fread("data/TAIR10/arabidopsis_genes.csv")

## 3) TAIR10 GFF3: genes + CDS -------------------------------------------

gff <- read.GFF("data/TAIR10/TAIR10_GFF3_genes.gff")
gff$CHROM <- gsub("Chr", "", gff$CHROM)

## 3a) CDS features, harmonized to CDS fasta names
CDS <- gff[TYPE == "CDS"]
# extract gene / model names from INFO
CDS$gene  <- gsub("\\..+", "", gsub("Parent=(.+),.+", "\\1", CDS$INFO))
CDS$model <- gsub("Parent=(.+),.+", "\\1", CDS$INFO)
CDS$model <- gsub("-Protein", "", CDS$model)

# keep only CDS entries for which we actually have a sequence
CDS <- CDS[model %in% names(cds)]
CDS$ID <- seq_len(nrow(CDS))

message(sprintf("[load_MA_data] Loaded %d CDS features.", nrow(CDS)))

## 3b) gene features (restricted to chr 1–5 and to genes with CDS)
genes <- gff[TYPE == "gene"]
genes$gene <- gsub(".+Name=", "", genes$INFO)
genes <- genes[gene %in% CDS$gene]   # only genes we have CDS for
genes <- genes[CHROM %in% 1:5]       # keep nuclear/chrom 1–5
genes$ID <- seq_len(nrow(genes))

message(sprintf("[load_MA_data] Loaded %d gene features (chr1–5).", nrow(genes)))

## 4) intergenic regions --------------------------------------------------
## defined as gaps between ordered genes on each chromosome

genes <- genes[order(CHROM, START)]

intergenic <- genes[, .(
  SOURCE = "TAIR10",
  TYPE   = "intergenic",
  START  = STOP + 1L,
  STOP   = data.table::shift(START, type = "lead") - 1L
), by = CHROM]

# drop trailing NA and inverted intervals
intergenic <- intergenic[!is.na(STOP)]
intergenic <- intergenic[START <= STOP]
intergenic$ID <- seq_len(nrow(intergenic))

message(sprintf("[load_MA_data] Defined %d intergenic intervals.", nrow(intergenic)))

## 5) precomputed 3-mer frequencies --------------------------------------

message("[load_MA_data] Loading precomputed 3-mer counts...")
genes_triN      <- fread("data/TAIR10/genes_trimers.csv")
intergenic_triN <- fread("data/TAIR10/intergenic_trimers.csv")
genome_triN     <- fread("data/TAIR10/genome_trimers.csv")

## 6) combined genic + intergenic table (for quick lookups) --------------

genic_intergenic <- rbind(genes, intergenic, fill = TRUE)
genic_intergenic$ID <- seq_len(nrow(genic_intergenic))

message(sprintf("[load_MA_data] Combined genic+intergenic: %d rows.", nrow(genic_intergenic)))

## 7) annotated mutation table -------------------------------------------

message("[load_MA_data] Loading annotated mutation table...")
mutations <- fread("Tables/mutations_Athal_annotated.csv")
message(sprintf("[load_MA_data] Loaded %d mutations.", nrow(mutations)))

message("[load_MA_data] Done.")


# mutation probabilities and effects --------------------------------------
codon_mut_counts<-fread("tables/codon_mut_counts.csv")

