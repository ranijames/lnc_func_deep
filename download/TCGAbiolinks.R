if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("BioinformaticsFMRP/TCGAbiolinks")
library(TCGAbiolinks)
library("limma")
library("edgeR")
library(SummarizedExperiment)


###
source("http://www.bioconductor.org/biocLite.R")
library(TCGAbiolinks)
cancer <- "BRCA"
PlatformCancer <- "RNA-seq"
dataType <- "HTSeq - FPKM-UQ"
pathCancer <- "/home/alva/Documents/TCGS/braca"

datQuery <- TCGAquery(tumor = cancer, platform = PlatformCancer, level = "3")  
lsSample <- TCGAquery_samplesfilter(query = datQuery)

# get subtype information
dataSubt <- TCGAquery_subtype(tumor = cancer)
lumA <- dataSubt[which(dataSubt$PAM50.mRNA == "Luminal A"),1]
allSamples <- lsSample$IlluminaHiSeq_RNASeqV2 #1218 total samples
lumASamples <- allSamples[grep(x = allSamples, pattern = paste(lumA, collapse = "|"))] # 263 luminal samples found

# get the subtypes 13.07.20
## ********************* https://bioc.ism.ac.jp/packages/3.2/bioc/vignettes/TCGAbiolinks/inst/doc/tcgaBiolinks.html#tcgaquery_subtype-working-with-molecular-subtypes-data **********
BRCA_path_subtypes <- TCGAquery_subtype(tumor = "brca")

### https://costalab.ukaachen.de/open_data/Bioinformatics_Analysis_in_R_2019/BIAR_D3/handout.html

query_TCGA = GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  experimental.strategy = "RNA-Seq",
  workflow.type = "HTSeq - FPKM-UQ",
  sample.type = "Primary Tumor")

# download the TCGA BRCA data for gene expression
GDCdownload(query = query_TCGA)
tcga_data = GDCprepare(query_TCGA)
# gene expression matrix
dim(assay(tcga_data))
data_matrix =assay(tcga_data)
write.table(data_matrix, "/home/alva/Documents/TCGS/braca/BRCA_download/GDCdata/TCGA-BRCA/harmonized/Transcriptome_Profiling/Gene_Expression_Quantification/BRCA_all.txt",row.names=T,quote=F,sep="\t")

