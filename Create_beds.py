import os
import sys
import argparse
from Functions import *

parser = argparse.ArgumentParser()
parser.add_argument('-fp', '--Fasta_pseudo', help = 'Pseudogene bed file.')
parser.add_argument('-bl', '--blast_database', help = 'Blast database of only genes.')
parser.add_argument('-OM', '--output_masked', help = 'Where the bed files with camo regions will be.')
parser.add_argument('-OR', '--output_realigned', help = 'Where the realign bed files will be created.')
parser.add_argument('-Id_threshold', '--Identity_threshold', help = 'Maximum sequence identity the homologous regions taken into account can have.', type = float, default = 95.0)
parser.add_argument('-pl', '--ploidy', help = 'Maximum number of ploidy (10 by default)', type = int, default = 10)
args = parser.parse_args()

Temporary_file = 'TEMP_DIR_B_AL'
os.system('mkdir ' + Temporary_file)
Pseudo_fasta = DivideFasta(args.Fasta_pseudo)

pseudogenes = Pseudo_fasta.get_seq_names()

for pseudo in pseudogenes:
    print('pseudo gene = ', pseudo)
    #Find those genes that have as sequence identity higher or equal to 90.
    Make_beds(pseudo, Pseudo_fasta, Temporary_file, args.blast_database, args.Identity_threshold, args.output_realigned, args.output_masked, args.ploidy)

os.system('rm -r ' + Temporary_file)





