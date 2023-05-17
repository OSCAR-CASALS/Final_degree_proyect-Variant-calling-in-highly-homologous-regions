# Variant Calling on genes with high homology to pseudogenes
This pipeline is designed to perform variant calling on regions with high homology to a pseudogene.

# Required Input
This pipeline requires:
  - A reference genome.
  - A gff3 annotation of the reference genome.
  - A fastq file with the reads we want to do variant calling on.

# How does it work?
First the pipeline will create a blast database 
