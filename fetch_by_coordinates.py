import copy
from Bio import Entrez, SeqIO
from Scripts.fetch_genomic_seq import fetch_dna_coordinates

def fetch_by_coordinates(genome,chrom,mut_pos,temp_dir):
    sequence=fetch_dna_coordinates(genome, chrom, mut_pos - 25, mut_pos + 25, temp_dir)
    return sequence


# aa=fetch_by_coordinates("eboVir3", 16, 77367589 ,'D:/Documents/ASchool/year 3/sem2/blog/user_files')
# print (aa)
# print (len(aa))