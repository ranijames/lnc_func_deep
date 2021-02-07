library("TCGAbiolinks")
library("limma")
library("edgeR")
library("TCGAbiolinks")
library("SummarizedExperiment")
library("survival")
library("survminer")
library("RColorBrewer")

library("genefilter")
GDCprojects  = getGDCprojects()
args         = commandArgs(TRUE)
# TCGA-BRCA, TCGA-KIRC, TCGA-THCA,TCGA-LUAD,TCGA-PRAD"
# query for BRCA helathy and normal
#tcga_cancers =list("TCGA-LIHC","TCGA-LUSC","TCGA-HNSC","TCGA-STAD","TCGA-KIRP","TCGA-COAD","TCGA-KICH","TCGA-BLCA","TCGA-ESCA","TCGA-PAAD")

# cd /home/alva/Documents/lncRNAs_project_2020/Analysis_2021/Tcga_R/GDCdata/TCGA-BRCA/harmonized/Transcriptome_Profiling/Gene_Expression_Quantification/
# Rscript /home/alva/Documents/lncRNAs_project_2020/scripts_2021/de_analysis/download_tcga_data.R TCGA-LIHC /home/alva/Documents/lncRNAs_project_2020/Analysis_2021/Tcga_R/GDCdata/TCGA-BRCA/harmonized/Transcriptome_Profiling/Gene_Expression_Quantification/r_data_new/   
                       
#"TCGA-STAD","TCGA-KIRP","TCGA-COAD","TCGA-KICH","TCGA-BLCA","TCGA-ESCA","TCGA-PAAD
# Function to download the files
cancertype= c("TCGA-KIRP")
cancertype           = as.character(args[1])
output_path          = args[2]
tcga_cancers         = list("TCGA-LIHC","TCGA-LUSC","TCGA-HNSC","TCGA-STAD","TCGA-KIRP","TCGA-COAD","TCGA-KICH","TCGA-BLCA","TCGA-ESCA","TCGA-PAAD")

query_TCGA = GDCquery(project  = cancertype, data.category = "Transcriptome Profiling", experimental.strategy = "RNA-Seq",workflow.type = "HTSeq - Counts", sample.type = c("Primary Tumor", "Solid Tissue Normal"))
project   = query_TCGA$project
#GDCdownload(query = query_TCGA)
tcga_data = GDCprepare(query_TCGA)

# sample ids or barcode
bar          = tcga_data@colData$barcode
Solid_tumor  = TCGAquery_SampleTypes(bar,"TP")
Solid_normal = TCGAquery_SampleTypes(bar,"NT")
# create a vector of Tumor and normal samples names
x                 = c(Solid_tumor,Solid_normal)
matrix            = assay(tcga_data)
norm_tum_matrix   = matrix[,x]

write.table(met,file=paste0(output_path,"matrix/",cancertype,".tsv"),quote=FALSE,sep="\t")

# Get the project ids
#project      = tcga_data@colData$project_id
#unique(project) 


# get the row read count for each cancertype

coldata      = as.data.frame(colData(tcga_data))
metadata     = coldata[,c("barcode","project_id","sample_type","tissue_or_organ_of_origin")]

met          = metadata[metadata$project_id==cancertype,]
write.table(met,file=paste0(output_path,"metadata/",cancertype,"_metadata"),quote=FALSE,sep="\t")

# saving rds
tcga = tcga_data[,tcga_data$project_id == cancertype]
saveRDS(object = tcga,
        file   = paste0(output_path,"metadata/",cancertype,".RDS"),
        compress = FALSE)

# for normal and tumor




