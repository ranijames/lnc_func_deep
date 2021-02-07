library("TCGAbiolinks")
library("limma")
library("edgeR")

library("SummarizedExperiment")
#library("gplots")
library("survival")
library("survminer")
library("RColorBrewer")
#library("gProfileR")
library("genefilter")
GDCprojects  = getGDCprojects()
args         = commandArgs(TRUE)
# runing example cat cancertypes | xargs -I% echo Rscript test.R % /home/alva/Documents/lncRNAs_project_2020/Analysis_2021/Tcga_R/GDCdata/TCGA-BRCA/harmonized/Transcriptome_Profiling/Gene_Expression_Quantification/r_data_new/ /home/alva/Documents/lncRNAs_project_2020/Analysis_2021/Tcga_R/GDCdata/TCGA-BRCA/harmonized/Transcriptome_Profiling/Gene_Expression_Quantification/r_data_new/%.RDS
# the limma pipeline for BRCA  DE genes

cancertype      = as.character(args[1])
outdir          = args[2]
input           = readRDS(file=args[3])
#
outdir     = '/home/alva/Documents/lncRNAs_project_2020/Analysis_2021/Tcga_R/GDCdata/TCGA-BRCA/harmonized/Transcriptome_Profiling/Gene_Expression_Quantification/r_data_new/'
#cancertype = c("TCGA-LUSC")
#input      = readRDS(file = "/home/alva/Documents/lncRNAs_project_2020/Analysis_2021/Tcga_R/GDCdata/TCGA-BRCA/harmonized/Transcriptome_Profiling/Gene_Expression_Quantification/r_data_new/LIHC_LUSC.RDS")

coldata      = as.data.frame(colData(input))
metadata     = coldata[,c("barcode","project_id","sample_type","tissue_or_organ_of_origin")]
met          = metadata[metadata$project_id==cancertype,]
write.table(met,file=paste0(outdir,"metadata/",cancertype,"_metadata"),quote=FALSE,sep="\t")

#################
limma_pipeline = function(
  tcga_data,
  condition_variable,
  reference_group=NULL,cancertype,outdir){
   g            = tcga_data[,tcga_data$project_id == cancertype]

  design_factor =  colData(g)[, condition_variable, drop=T] 
  group         = factor(design_factor)

  if(!is.null(reference_group)){group = relevel(group, ref=reference_group)}

  design = model.matrix(~ group)

  dge    = DGEList(counts=assay(g),
                 samples=colData(g),
                 genes=as.data.frame(rowData(g)))

  # filtering
  keep = filterByExpr(dge,design)
  dge  = dge[keep,,keep.lib.sizes=FALSE]
  rm(keep)

  # Normalization (TMM followed by voom)
  dge = calcNormFactors(dge)
  v   = voom(dge, design, plot=TRUE)

  # Fit model to data given design
  fit = lmFit(v, design)
  fit = eBayes(fit)

  # Show top genes
  topGenes = topTable(fit, coef=ncol(design), number = nrow(g), sort.by="none")
  normalized_matrix=v$E
  resdata_single = merge(as.data.frame(topGenes), normalized_matrix, by = 'row.names', sort = FALSE)

  #result_2<- topTable(fit, number = nrow(data), sort.by ="none")

  write.table(resdata_single,file=paste0(outdir,"de/",cancertype,"_DE"),quote=FALSE,sep="\t",row.names = F)

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
  tcga_data=input,
  condition_variable="definition",
  reference_group="Solid Tissue Normal",
  cancertype = cancertype,
  outdir=outdir
)


