# Variant calling in highly homologous regions
This pipeline is designed to perform variant calling on regions with high homology to a pseudogene by masking all of them with the exeption of one that will represent the rest.

## Required Input
This pipeline requires:
  - A reference genome.
  - A gff3 annotation of the reference genome.
  - A fastq file with the reads we want to do variant calling on.
  
## How to use it?
First it is required to execute the script **Init.sh** to create a blast database of the genes present in the reference genome and a fasta file with the sequences of the pseudogene, this script requires a gff3 annotation file, a reference genome and the output directory as input(in that order).

Once Init.sh has finished creating the necessary files for the pipeline the script **Create_beds.py** must be executed. This script will create for each ploidy a bed file with the regions that need to be realligned(realign_to_ploidy.bed) and another one with their representative regions(align_to_ploidy.bed). The script requires the pseudogene fasta file and the blast database created earlier as well as the output directories were the bed files will be created, optionally the maximum percentage of sequence identity the homologous regions taken into account can have (95% by default) and the maximum number of ploidys considered (by default 10) can be specified. 

With the bed files created **Find_variants.py** must be executed to mask the genome and look for variants in it. This script requires as input: the reads we want to look for variants in, the ID and SM of the samples, the number of threads that will be used in the allignments and variant calling inside the script, the name of the vcf output file, the path to the directory of the align_to and realign_to files, and the reference genome; optionally the phred-scaled confidence threshold can be specified (900 by default).

Finally the script **Find_hom_zones.py** must be executed to flag the variants in the vcf output who present high homology to a pseudogene. This script requires the vcf file created earlier, the directory containing the Align_to files and the name of the generated vcfs.

Instead of executing each script by order this github contains a python script called **Execute_all.py** that will automatically execute all of them returning only the final vcf file. This script requires the reference genome, the annotation file of the genome used,reads were to look for variants, the ID and SM of the sample, the name of the vcf output file, the number of threads (1 by default). Optionally the phred-scaled confidence threshold of the variants, the maximum sequence identity a high homology region can have to a pseudogene and the maximum of ploidies considered can be specified (10 by default).
