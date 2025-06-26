# Differential Peak Analysis with DiffBind in R

This guide walks you through the process of performing differential peak analysis using the **DiffBind** package in R, followed by downstream annotation and gene ontology enrichment.

---

## 🧾 Inputs Required

- **Sample sheet CSV**: Metadata for your samples including paths to BAM and BED/peak files.
- **Peak files** (e.g., BED or narrowPeak): One per sample.
- **Optional**: BAM files for read counting.

---

## 1. Install and Load Required Packages

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("DiffBind", "ChIPseeker", "TxDb.Hsapiens.UCSC.hg19.knownGene", 
                       "org.Hs.eg.db", "clusterProfiler", "rtracklayer", "ggplot2", "VennDiagram"))

library(DiffBind)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(clusterProfiler)
library(rtracklayer)
library(ggplot2)
```

---

## 2. Prepare the Sample Sheet

Create a `sample_sheet.csv` file with this structure:

```
SampleID,Tissue,Factor,Condition,Replicate,bamReads,ControlID,Peaks,PeakCaller
ASL1,DRG,H3K4me3,ASL,1,ASL1.bam,NA,ASL1.bed,bed
CTRL1,DRG,H3K4me3,CTRL,1,CTRL1.bam,NA,CTRL1.bed,bed
```

---

## 3. Load the Sample Sheet and Create DBA Object

```r
samples <- read.csv("sample_sheet.csv")
dba_obj <- dba(sampleSheet = samples)
```

**Input**: sample_sheet.csv  
**Output**: DiffBind object (`dba_obj`)

---

## 4. Count Reads (Optional if BAM files are unavailable)

```r
dba_obj <- dba.count(dba_obj)
```

**Input**: BAM files  
**Output**: Read counts per peak

---

## 5. Set Contrasts

```r
dba_obj <- dba.contrast(dba_obj, categories = DBA_CONDITION, minMembers = 2)
```

**Input**: DBA object with read counts  
**Output**: Contrast information added

---

## 6. Perform Differential Analysis

```r
dba_obj <- dba.analyze(dba_obj)
dba_report <- dba.report(dba_obj, th = 0.05)
```

**Input**: DBA object with contrast  
**Output**: Differential peak GRanges object

---

## 7. Export Differential Peaks

```r
write.csv(as.data.frame(dba_report), "DiffBind_results.csv")
```

**Output**: `DiffBind_results.csv`

---

## 8. Annotate Peaks

```r
peakAnno <- annotatePeak(dba_report, TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene, 
                         tssRegion = c(-3000, 3000), annoDb = "org.Hs.eg.db")
plotAnnoPie(peakAnno)
```

**Output**: Annotated peaks, pie chart

---

## 9. Classify Peaks: Gained / Lost / Unchanged

```r
peak_df <- as.data.frame(dba_report)
peak_df$classification <- ifelse(peak_df$Fold > 0 & peak_df$FDR < 0.05, "Gained",
                           ifelse(peak_df$Fold < 0 & peak_df$FDR < 0.05, "Lost", "Unchanged"))

table(peak_df$classification)
write.csv(peak_df, "classified_peaks.csv")
```

**Output**: `classified_peaks.csv`

---

## 10. PCA and Heatmap

```r
plot(dba_obj)                  
dba.plotHeatmap(dba_obj)
```

**Output**: PCA and correlation heatmap plots

---

## 11. GO Enrichment for Gained Peaks

```r
gained_genes <- as.data.frame(peakAnno)$geneId[peak_df$classification == "Gained"]
gained_genes <- na.omit(unique(gained_genes))

ego <- enrichGO(gene = gained_genes, OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH")

dotplot(ego, showCategory = 10)
```

**Output**: Dotplot of enriched GO terms

---

## Optional: Export GRanges Object

```r
saveRDS(dba_report, file = "DiffBind_GRanges.rds")
```

**Output**: `DiffBind_GRanges.rds`

---

## 🔚 Outputs Summary

- `DiffBind_results.csv`: Differential peak results
- `classified_peaks.csv`: Peaks classified as gained/lost/unchanged
- `DiffBind_GRanges.rds`: GRanges object
- Plots: PCA, heatmap, GO enrichment

