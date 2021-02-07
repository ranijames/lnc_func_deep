library("TCGAbiolinks")
library("limma")
library("edgeR")
library("glmnet")
library("factoextra")
library("FactoMineR")
library("caret")
library("SummarizedExperiment")
library("gplots")
library("survival")
library("survminer")
library("RColorBrewer")
library("gProfileR")
library("genefilter")
GDCprojects = getGDCprojects()
args         = commandArgs(TRUE)
# loading R object

tcga_data = readRDS(file = "tcga_data.RDS")

# the limma pipeline for BRCA  DE genes

limma_pipeline = function(
  tcga_data,
  condition_variable,
  reference_group=NULL){

  design_factor = colData(tcga_data)[, condition_variable, drop=T]
 
# alternative
design=colData(tcga_data)[,c("definition","project_id"),drop=T]

  group = factor(design_factor)
  if(!is.null(reference_group)){group = relevel(group, ref=reference_group)}

  design = model.matrix(~ group)

  dge = DGEList(counts=assay(tcga_data),
                 samples=colData(tcga_data),
                 genes=as.data.frame(rowData(tcga_data)))

  # filtering
  keep = filterByExpr(dge,design)
  dge = dge[keep,,keep.lib.sizes=FALSE]
  rm(keep)

  # Normalization (TMM followed by voom)
  dge = calcNormFactors(dge)
  v = voom(dge, design, plot=TRUE)

  # Fit model to data given design
  fit = lmFit(v, design)
  fit = eBayes(fit)

  # Show top genes
  topGenes = topTable(fit, coef=ncol(design), number = nrow(tcga_data), sort.by="none")
#result_2<- topTable(fit, number = nrow(data), sort.by ="none")

  write.table(topGenes,file="/home/alva/Documents/lncRNAs_project_2020/Analysis_2021/Tcga_R/GDCdata/TCGA-BRCA/harmonized/Transcriptome_Profiling/Gene_Expression_Quantification/r_data/PRAD_DE",quote=FALSE,sep="\t",row.names = F)

  return(
    list(
      voomObj=v, # normalized data
      fit=fit, # fitted model and statistics
      topGenes=topGenes # All differentially expressed genes
    )
  )
}

# Run the above pipeline


limma_res = limma_pipeline(
  tcga_data=tcga_data,
  condition_variable="definition",
  reference_group="Solid Tissue Normal"
)




