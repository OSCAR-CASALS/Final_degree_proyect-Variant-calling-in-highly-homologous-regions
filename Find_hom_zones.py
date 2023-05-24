import os
import sys
import argparse
from Divide_fasta_3 import Flagg_high_homology_regions

parser = argparse.ArgumentParser()
parser.add_argument('VCF', help = 'VCF file directory')
parser.add_argument('Al', help = 'Align_to bed file directory.')
parser.add_argument('Res', help = 'Prefix of the generated vcf files.')
args = parser.parse_args()

hom_regions = [] # List of the regions extracted from the align_to file that are considered to be of high homology.
directory_align_to = os.listdir(args.Al) # All Align_to files
al_to_dir = args.Al # Directory were hte align to files are kept.

for bed in directory_align_to:
    #In the columns 6 and 7 of each align_to file we can find the exact region that is aligned to the pseudogene when creating the align_to and realign files.
    #In this loop we save the regions in a list of high homology regions.
    b_file = open(al_to_dir + '/' + bed, 'r')
    line = b_file.readline()
    while line != '':
        l = line.split('\t')
        hom_regions += [[l[6], l[7]]]
        line = b_file.readline()
    b_file.close()

#Finally I iterate through the vcfs adding a flag to those variants that are in a high homology region.
directory_vcf = args.VCF + '/'
vcfs = os.listdir(directory_vcf)

for v in vcfs:
    if (v.split('.')[-1] != 'idx'):
        pl = v.split('_')[-1] 
        Flagg_high_homology_regions(directory_vcf + v, hom_regions, args.Res + pl)