import sys
import os
import argparse
from Functions import Flagg_high_homology_regions

parser = argparse.ArgumentParser()
parser.add_argument('-I', '--input', help = 'Samples where we want to locate variants in pseudogenes (Only one fasta)')
parser.add_argument('-ID', '--ID_RG_group', help = 'ID of the sample to be used')
parser.add_argument('-SM', '--SM_RG_group', help = 'SM of the sample to be used')
parser.add_argument('-t', '--threads', help = 'Number of threads that will be used.')
parser.add_argument('-O', '--output', help = 'Path were the gvcf files will be created.')
parser.add_argument('-real', '--realign_bed', help = 'Directory where the realign files are.')
parser.add_argument('-Al_to', '--Align_to', help = 'Directory where the align to beds are.')
parser.add_argument('-R', '--REFERENCE', help = 'Reference genome.')
parser.add_argument('-Qual', '--Quality', help = 'Quality threshold', default='900.0')
args = parser.parse_args()

#Aligning samples to reference genome so I can extract the regions that need to be realigned.
rG_group = '@RG' + r'\t' + 'ID:' + args.ID_RG_group + r'\t' + 'SM:' + args.SM_RG_group + r'\t' + 'PL:ILLUMINA'
RG_group = "'" + rG_group + "'"

print('RG group = ', RG_group)

TENPORARY_DICT = 'tmp_Variant_calling_pseudogenes'

os.system('mkdir ' + TENPORARY_DICT)

print('Aligning to reference genome...')

os.system('bwa mem -M ' + '-R ' + RG_group + ' -t ' + args.threads + ' ' + args.REFERENCE + ' ' + args.input + ' | samtools view -hb | samtools sort -@ ' + args.threads + ' > ' + TENPORARY_DICT + '/basic_al.bam')

print('Checking if the bam, file is correct...')
os.system('samtools quickcheck ' + TENPORARY_DICT + '/basic_al.bam')

print('Finished checking', 'Indexing...', sep = '\n')
os.system('samtools index ' + TENPORARY_DICT + '/basic_al.bam')


def extract_ploidy(name):
    '''
    This function extracts the ploidy of the align_to and realign_to files.
    '''
    return int(name.split('_')[-1][0:-4])

#Making sure the align_to and realign files have the same positions in both lists
directory_realign = args.realign_bed
realign_files = sorted(os.listdir(directory_realign), key=extract_ploidy)

directory_align_to = args.Align_to
align_to_files = sorted(os.listdir(directory_align_to), key=extract_ploidy)

for rea in range(0, len(realign_files)):
    al_to = directory_align_to + '/' + align_to_files[rea]
    real_to = directory_realign + '/' + realign_files[rea]
    ploidy = str(extract_ploidy(real_to))
    print('Ploidy', ploidy)
    os.system('samtools view -h -L ' + real_to + ' ' + TENPORARY_DICT + '/basic_al.bam' + ' | ' + 'samtools view -hb | samtools sort -n -@ ' + args.threads + ' > ' + TENPORARY_DICT + '/sorted_bam_' + ploidy + '.bam')
    TEMP_FASTA = TENPORARY_DICT + '/' + 'ploidy_' + ploidy + '.fasta'
    os.system('bedtools bamtofastq -i ' + TENPORARY_DICT + '/sorted_bam_' + ploidy + '.bam ' + '-fq ' + TEMP_FASTA)
    #sys.exit()
    os.system('rm ' + TENPORARY_DICT + '/sorted_bam_' + ploidy + '.bam')
    #Masking
    out_c = TENPORARY_DICT + '/' + 'Masked_' + ploidy
    os.system('mkdir ' + out_c)
    TEMPORARY_FILE_MASKED = TENPORARY_DICT + '/' + 'Temporary_masked_ploidy_' + ploidy
    os.system('bash Mask_2.sh ' + al_to + ' ' + args.REFERENCE + ' ' + out_c + ' ' + TEMPORARY_FILE_MASKED)
    MASKED_AL = out_c + '/Masked_alignment.bam'
    os.system('bwa mem -M -R ' + RG_group + ' ' + '-t ' + args.threads + ' ' + out_c + '/Masked_gen.fasta' + ' ' +  TEMP_FASTA + ' | samtools view -hb | samtools sort -@ ' + args.threads + ' > ' + MASKED_AL)
    os.system('samtools quickcheck ' + MASKED_AL)
    os.system('samtools index ' + MASKED_AL)
    #Haplotype caller
    HAPL_creation = out_c + '/Hapl_' + ploidy + '.g.vcf'
    os.system('gatk HaplotypeCaller -R ' + out_c + '/Masked_gen.fasta ' + '-I ' + MASKED_AL + ' -ERC GVCF --dont-use-soft-clipped-bases -O ' + HAPL_creation + ' -ploidy ' + ploidy)
    #Genotyping
    GEN_out = TENPORARY_DICT + '/Genotyping_' + ploidy + '.vcf'
    os.system('gatk GenotypeGVCFs -R ' + out_c + '/Masked_gen.fasta ' + '-O ' + GEN_out + ' --variant ' + HAPL_creation + ' -ploidy ' + ploidy + ' -stand-call-conf ' + args.Quality)
    #Removing files
    os.system('rm -r ' + out_c)
print('Fusing vcfs of all ploidys...')
os.system('vcf-concat ' + TENPORARY_DICT + '/Genotyping_*.vcf > ' + args.output)
print('Removing temporal file')
os.system('rm -r ' + TENPORARY_DICT)

