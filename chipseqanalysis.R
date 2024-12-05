### chiqseq analysis
setwd("/Users/cristalvillalba/Downloads/jimmychipseq/rep1H3K27me3")
# Load libraries
library(ChIPseeker)
library(clusterProfiler)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(GenomicRanges)
library(ggplot2)
library(rtracklayer)

# Read the narrowPeak file
narrowPeak <- read.table("chipseq_analysis_peaks.narrowPeak", header = FALSE)

# Assign column names
colnames(narrowPeak) <- c("chrom", "start", "end", "name", "score", "strand", 
                          "signal", "pval", "qval", "summit")

# Convert peaks to GRanges object
peaks <- GRanges(
  seqnames = narrowPeak$chrom,
  ranges = IRanges(start = narrowPeak$start, end = narrowPeak$end),
  strand = "*"
)


# Load the genome annotation database
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Annotate peaks to genes
peak_annotation <- annotatePeak(
  peaks, 
  tssRegion = c(-3000, 3000), # Define upstream/downstream TSS region
  TxDb = txdb,
  annoDb = "org.Hs.eg.db"
)

# View summary of annotation
print(peak_annotation)
# Pie chart of genomic annotations
plotAnnoPie(peak_annotation)

# Barplot of genomic annotations
plotAnnoBar(peak_annotation)

# Distribution of peaks relative to TSS
plotDistToTSS(peak_annotation, title = "Peak Distribution Around TSS")
# Convert annotated peaks to a data frame
annotated_peaks_df <- as.data.frame(peak_annotation)

# Save to a CSV file
write.csv(annotated_peaks_df, "annotated_peaks.csv", row.names = FALSE)
# Extract Entrez IDs of genes associated with peaks
gene_ids <- unique(annotated_peaks_df$geneId)

# Perform GO enrichment analysis
go_enrichment <- enrichGO(
  gene = gene_ids,
  OrgDb = org.Hs.eg.db,
  ont = "BP", # Biological Process
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = TRUE
)

# View enrichment results
head(go_enrichment)

# Plot top 10 enriched GO terms
barplot(go_enrichment, showCategory = 10)
# Perform KEGG pathway enrichment analysis
kegg_enrichment <- enrichKEGG(
  gene = gene_ids,
  organism = 'hsa', # 'hsa' for human; use 'mmu' for mouse
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05
)

# View KEGG enrichment results
head(kegg_enrichment)

# Plot top 10 enriched pathways
barplot(kegg_enrichment, showCategory = 10)
# Export peaks to BED format
export(peaks, "chipseq_peaks.bed", format = "BED")
# Visualize peak distribution across chromosomes
covplot(peaks, weightCol = "score")
# Example: Compare two peak sets
peaks2 <- GRanges(
  seqnames = c("chr1", "chr2"),
  ranges = IRanges(start = c(15000, 25000), end = c(20000, 30000)),
  strand = "*"
)

# Find overlaps
overlaps <- findOverlaps(peaks, peaks2)

# Extract overlapping regions
overlap_peaks <- peaks[queryHits(overlaps)]



########### using the DEGs provided by Oscar
# Load libraries
library(ChIPseeker)
library(clusterProfiler)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(JASPAR2022)
library(TFBSTools)
library(rtracklayer)
library(GenomicRanges)

## To test I will use just the top 100 degs from the list
# Define gene list
gene_list <- c("TSPAN6", "DPM1", "CFH", "FUCA2", "GCLC", "NFYA", "NIPAL3", "LAS1L", 
               "SEMA3F", "RAD52", "BAD", "LAP3", "CD99", "HS3ST1", "LASP1", "SNX11", 
               "KLHL13", "CYP26B1", "ICA1", "ALS2", "CASP10", "CFLAR", "RBM5", "SLC7A2", 
               "ARF5", "SARM1", "PLXND1", "FKBP4", "KDM1A", "RBM6", "CAMKK1", "ARHGAP33", 
               "PDK4", "SLC25A13", "CDC27", "SKAP2", "LIG3", "RPAP3", "REXO5", "CIAPIN1", 
               "SPPL2B", "PRKAR2B", "MSL3", "WDR54", "CROT", "FBXL3", "PDK2", "ITGA3", 
               "LAMP2", "GDE1", "TMEM98", "MAP3K14", "TMEM132A", "AP2B1", "CX3CL1", 
               "SPATA20", "TNFRSF12A", "RALA", "BAIAP2L1", "KDM7A", "ETV1", "PHTF2", 
               "FARP2", "GGCT", "TBXA2R", "IFRD1", "VPS41", "ELAC2", "ADIPOR2", "CCDC124", 
               "BLTP2", "TSR3", "PIGQ", "TEAD3", "SELE", "DNAJC11", "MYLIP", "NADK", 
               "CYTH3", "SYPL1", "SPAG9", "CELSR3", "MGST1", "CRY1", "NFIX", "ST3GAL1", 
               "IL32", "PKD1", "HEATR5B", "REV3L", "POMT2", "VTA1", "ETV7", "METTL13", 
               "STARD3NL", "NCAPD2", "SEMA3G", "NISCH", "STAB1", "ZNF200")

# Retrieve gene annotations (hg38 genome)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
gene_annot <- AnnotationDbi::select(org.Hs.eg.db, keys = gene_list, 
                                    keytype = "SYMBOL", columns = c("SYMBOL", "ENTREZID"))
head(gene_annot)

# Filter genes with available genomic regions
genes_gr <- genes(txdb, columns = c("gene_id"))
selected_genes <- genes_gr[genes_gr$gene_id %in% gene_annot$ENTREZID]
# Import peaks (example narrowPeak file from MACS2)
peak_file <- "chipseq_analysis_peaks.narrowPeak"
peaks <- rtracklayer::import(peak_file)

# Annotate peaks near genes
peak_annotation <- annotatePeak(
  peaks,
  tssRegion = c(-3000, 3000), 
  TxDb = txdb, 
  annoDb = "org.Hs.eg.db"
)

# Filter annotated peaks for genes in the list
filtered_peaks <- peak_annotation@anno[peak_annotation@anno$SYMBOL %in% gene_list, ]
head(filtered_peaks)
# Load JASPAR motif database
jaspar_db <- getMatrixSet(
  JASPAR2022, 
  opts = list(collection = "CORE", tax_group = "vertebrates", all_versions = FALSE)
)

# Match motifs to peaks
seq_regions <- getAllPeakSeq(filtered_peaks, genome = "hg38")
motif_enrichment <- findMotifsGenome(
  input = seq_regions,
  genome = "hg38",
  motif_database = jaspar_db,
  pValueCutoff = 0.05
)

# View enriched motifs
head(motif_enrichment)
# Perform transcription factor enrichment analysis
tf_enrichment <- enrichGO(
  gene = gene_annot$ENTREZID,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = TRUE
)

# Plot enriched transcription factors
barplot(tf_enrichment, showCategory = 10)
# Export annotated peaks
write.csv(as.data.frame(filtered_peaks), "filtered_peaks.csv", row.names = FALSE)

# Export enriched motifs
write.csv(as.data.frame(motif_enrichment), "motif_enrichment.csv", row.names = FALSE)

# Export transcription factor enrichment
write.csv(as.data.frame(tf_enrichment), "tf_enrichment.csv", row.names = FALSE)
# Peak annotation pie chart
plotAnnoPie(peak_annotation)

# Barplot of enriched motifs
barplot(tf_enrichment, showCategory = 10)

# Genomic distribution of peaks
plotDistToTSS(peak_annotation, title = "Distribution of Peaks Relative to TSS")
