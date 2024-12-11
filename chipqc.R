### chip quality control retrieved from https://hbctraining.github.io/Intro-to-ChIPseq/lessons/06_combine_chipQC_and_metrics.html
library(ChIPQC)
samples <- read.csv('sample.csv')
## Create ChIPQC object
chipObj <- ChIPQC(sample, annotation="hg38") 
ChIPQCreport(chipObj, reportName="ChIP QC report: H3K4Me3 and H3K27me3 rep 1 and 2", reportFolder="ChIPQCreport")
