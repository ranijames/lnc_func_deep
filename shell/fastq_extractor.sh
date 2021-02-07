#!/bin/bash
while read gene; do
     echo "faidx  --regex $gene /home/alva/Documents/lncRNAs_project_2020/scripts_2021/TriplexFPP_data/data/gencode.v36.lncRNA_transcripts.fa >> /home/alva/Documents/lncRNAs_project_2020/Analysis_2021/Tcga_R/GDCdata/TCGA-BRCA/harmonized/Transcriptome_Profiling/Gene_Expression_Quantification/r_data_new/fasta_lncRNAs/accross_nine_cancer_fastq/$gene.fasta"
# "faidx  --regex $gene /home/alva/rgtdata/hg38/gencode.v36.lncRNA_transcripts.fa >>  /home/alva/Documents/lncRNAs_project_2020/Analysis_2021/Tcga_R/GDCdata/TCGA-BRCA/harmonized/Transcriptome_Profiling/Gene_Expression_Quantification/r_data/fasta_lncRNAs/$gene.fasta" > cmd.sh
    #faidx  --regex $gene /home/alva/rgtdata/hg38/gencode.v36.lncRNA_transcripts.fa >>  /home/alva/Documents/lncRNAs_project_2020/Analysis_2021/Tcga_R/GDCdata/TCGA-BRCA/harmonized/Transcriptome_Profiling/Gene_Expression_Quantification/r_data/fasta_lncRNAs/$gene.fasta;

done < /home/alva/Documents/lncRNAs_project_2020/Analysis_2021/Tcga_R/GDCdata/TCGA-BRCA/harmonized/Transcriptome_Profiling/Gene_Expression_Quantification/r_data_new/de/lncRNA_gene_sym.txt >cmd.sh

 # for each genes save seperately

