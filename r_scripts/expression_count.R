tcga_data = readRDS(file = "/home/alva/Documents/lncRNAs_project_2020/Analysis_2021/Tcga_R/GDCdata/TCGA-BRCA/harmonized/Transcriptome_Profiling/Gene_Expression_Quantification/r_data/tcga_data.RDS")

tcga <- getResults(tcga_data,cols=c("cases"))

dataSmTP.miR <- TCGAquery_SampleTypes(barcode = tcga,
                                  typesample = "TP")

dataSmNT.miR <- TCGAquery_SampleTypes(barcode = tcga,
                                  typesample = "NT")

read_countData <-  colnames(dataAssy.miR)[grep("count", colnames(dataAssy.miR))]
dataAssy.miR <- dataAssy.miR[,read_countData]
colnames(dataAssy.miR) <- gsub("read_count_","", colnames(dataAssy.miR))


