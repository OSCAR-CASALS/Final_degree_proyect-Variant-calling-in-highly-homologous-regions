import os
import sys
import argparse
from Divide_fasta_3 import *

parser = argparse.ArgumentParser()
parser.add_argument('-fp', '--Fasta_pseudo', help = 'Pseudogene bed file.')
parser.add_argument('-bl', '--blast_database', help = 'Blast database of only genes.')
parser.add_argument('-OM', '--output_masked', help = 'Where the bed file with camo regions will be.')
parser.add_argument('-OR', '--output_realigned', help = 'Where the realign bed will be created.')
parser.add_argument('-Reg', '--Region', help = 'Region belonging to the gene the reads come from(format -> chromosome:Start-End, for example NC_000001.11:155234452-155244627)(Optional but recommended for more exact results).', default = '')
parser.add_argument('-max_pl', '--max_ploidy', help = 'All regions with a ploidy higher to this one will nnot be taken into account(with the exeption of the specified region.)', type = int, default = 4)
parser.add_argument('-Id_thres', '--Percentage_identity_threshold', help = 'Only the regions with equal or below percentage of identity to a pseudogene will be take into account', type = float, default=95.0)
args = parser.parse_args()

#Checking that all required files are in the input of the script.

if args.Fasta_pseudo == None:
    print('ERROR: fasta file with the sequences belonging to pseudogenes required. You can use Init.sh to create it.')
    sys.exit()

if args.blast_database == None:
    print('ERROR: Blast database of genes required, you can use Init.sh to create it.')
    sys.exit()

if args.output_masked == None:
    print('ERROR: File where the Align_to regions will be kept required.')
    sys.exit()

if args.output_realigned == None:
    print('ERROR: File were the realign regions will be kept required.')
    sys.exit()

if args.Region == '':
    print('WARNING: No region specified, all mutations in the camouflaged groups that resemble your reads will be returned (high number of false positives).')
    sys.exit()

#Creating temp0orary directory
Temporary_file = 'TEMP_DIR_B_AL'
os.system('mkdir ' + Temporary_file)
#Creating a DivideFasta object that will allow us to extract sequenes from the Pseudogene fasta file, this will allow us to obtain the query sequence of the allignments.
Pseudo_fasta = DivideFasta(args.Fasta_pseudo)
#We align each pseudogene to the blast database, those alignments with a percentage of identity higher than 90 and lower or equal than the threshold will be considered as homologous regions and kept in the align_to and realign regions.
pseudogenes = Pseudo_fasta.get_seq_names()

for pseudo in pseudogenes:
    print('pseudo gene = ', pseudo)
    #Create Align_to and realign files.
    Get_best_alignment_blastn(pseudo, Pseudo_fasta, Temporary_file, args.blast_database, args.Percentage_identity_threshold, args.output_realigned, args.output_masked, args.max_ploidy, args.Reg)

os.system('rm -r ' + Temporary_file)





