#-------------------------------------- 
# DNA methylation data
#-------------------------------------- 
# DNA methylation aligned to hg38
library("TCGAbiolinks")

setwd("/home/alva/Documents/lncRNAs_project_2020/Analysis_2021/DNA_methy_tcga")
cancertype =list("TCGA-LIHC","TCGA-LUSC","TCGA-HNSC","TCGA-STAD","TCGA-KIRP","TCGA-COAD","TCGA-KICH","TCGA-BLCA","TCGA-ESCA","TCGA-PAAD")
can =list("TCGA-BLCA","TCGA-ESCA","TCGA-PAAD")
query_met.hg38 <- GDCquery(project= "TCGA-BRCA", 
                           data.category = "DNA Methylation", 
                           platform = "Illumina Human Methylation 450")
GDCdownload(query_met.hg38)
data.hg38 <- GDCprepare(query_met.hg38)

query <- GDCquery(project = "TCGA-GBM",
                  data.category = "DNA methylation",
                  platform = "Illumina Human Methylation 27",
                  legacy = TRUE, 
                  barcode = c("TCGA-02-0058-01A-01D-0186-05", "TCGA-12-1597-01B-01D-0915-05",
                              "TCGA-12-0829-01A-01D-0392-05", "TCGA-06-0155-01B-01D-0521-05",
                              "TCGA-02-0099-01A-01D-0199-05", "TCGA-19-4068-01A-01D-1228-05",
                              "TCGA-19-1788-01A-01D-0595-05", "TCGA-16-0848-01A-01D-0392-05"))
GDCdownload(query, method = "api")
data <- GDCprepare(query)

#----------- 8.3 Identification of Regulatory Enhancers   -------
library(TCGAbiolinks)
# Samples: primary solid tumor w/ DNA methylation and gene expression
project =c("TCGA-BRCA")

matched_met_exp <- function(project){
    # get primary solid tumor samples: DNA methylation
    message("Download DNA methylation information")
    met450k <- GDCquery(project = project,
                        data.category = "DNA Methylation",
                        platform = "Illumina Human Methylation 450",
                        sample.type = c("Primary Tumor", "Solid Tissue Normal"))
    met450k.tp <-  met450k$results[[1]]$cases
    
    # get primary solid tumor samples: RNAseq
    message("Download gene expression information")
    exp <- GDCquery(project = project,
                    data.category = "Transcriptome Profiling",
                    experimental.strategy = "RNA-Seq",
                    workflow.type = "normalized_results",,
                    sample.type = c("Primary Tumor", "Solid Tissue Normal"))
                    
    exp.tp <- exp$results[[1]]$cases
    # Get patients with samples in both platforms
    patients <- unique(substr(exp.tp,1,15)[substr(exp.tp,1,12) %in% substr(met450k.tp,1,12)] )
    #if(!is.null(n)) 
     #{patients <- patients[1:n]} # get only n samples
    return(patients)
}

samples <- matched_met_exp("TCGA-BRCA")

#-----------------------------------
# 1 - Methylation
# ----------------------------------
query.met <- GDCquery(project = "TCGA-BRCA",
                      data.category = "DNA Methylation",
                      platform = "Illumina Human Methylation 450",
                      barcode = samples)
GDCdownload(query.met)
met <- GDCprepare(query.met, save = FALSE)
#met <- subset(met)

#-----------------------------------
# 2 - Expression
# ----------------------------------
query.exp <- GDCquery(project = "TCGA-BRCA",
                     data.category = "Transcriptome Profiling",
                     platform = "Illumina HiSeq", experimental.strategy = "RNA-Seq",
                     file.type  = "normalized_results", 
                     barcode =  samples)

GDCdownload(query.exp)
exp <- GDCprepare(query.exp, save = FALSE)
save(exp, met, gbm.samples, lgg.samples, file = "brcamethexp.rda", compress = "xz")
"TCGA-STAD"
"TCGA-KICH"

#samples <- c(lgg.samples,gbm.samples)
