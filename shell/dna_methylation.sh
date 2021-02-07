cat /home/alva/Documents/lncRNAs_project_2020/Analysis_2021/DNA_methy_tcga/Differential-CpGs-2021-02-06_BRCA_hyper.csv| sed 's/,/\t/g'|sed 's/"//g'| cut -f9-11 >BRCA_hyper.bed
cat Differential-CpGs-2021-02-06_BRCA_hypo.csv |sed 's/,/\t/g'|sed 's/"//g'| cut -f9-11 >BRCA_hypo.bed
cat BRCA_hypo.bed BRCA_hyper.bed > BRCA.bed
AnnotatePeaks.pl BRCA_hyper.bed /home/alva/Documents/lncRNAs_project_2020/tools/data/genomes/hg38 -gene >BRCA_annaotted 
