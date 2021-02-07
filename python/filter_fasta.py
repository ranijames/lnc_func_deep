#!/usr/bin/env python
# python filter_fasta.py /home/alva/Documents/lncRNAs_project_2020/scripts_2021/TriplexFPP_data/data/header_remove /home/alva/Documents/lncRNAs_project_2020/scripts_2021/TriplexFPP_data/data/gencode.v36.lncRNA_transcripts.fa /home/alva/Documents/lncRNAs_project_2020/scripts_2021/TriplexFPP_data/data/gencode.v36.lncRNA_transcripts_filtered.fasta
import sys
#my_set = set([])
from Bio import SeqIO
# https://www.biostars.org/p/157811/
def list_ids():
    """
    Return a set containing the identifiers presented in a file,
    line by line, starting with ">"
    """

    # read the first file given and generate a set (faster iteration respect lists

    identifiers = list()

    with open(sys.argv[1], 'r') as fi:
        for line in fi:
            line = line.strip()
            identifiers.append(str(line).replace(">", ""))

    return identifiers

def filter():
    """
    Writes a file containing only the sequences with identifiers NOT
    present in a set of identifiers
    """

    identifiers = list_ids()

    with open(sys.argv[2]) as original_fasta, open(sys.argv[3], 'w') as corrected_fasta:
        records = SeqIO.parse(original_fasta, 'fasta')
        for record in records:
            #print (record.id)
            if record.id in identifiers:
                SeqIO.write(record, corrected_fasta, 'fasta')

if __name__ == '__main__':
    filter()