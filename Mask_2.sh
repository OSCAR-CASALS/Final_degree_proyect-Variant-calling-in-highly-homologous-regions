#!/bin/bash
align_to_bed=$1
ref=$2
out=$3
tmp=$4
if [ -z $align_to_bed ]
then
    echo "Bedfile with regions were pseudogenes should align to required."
fi
mkdir $tmp
RETURN_CODE=$?
if [[ $RETURN_CODE != 0 ]]
then
    exit $RETURN_CODE
fi
#Sorting bed
bedtools sort -i $align_to_bed > $tmp/Sorted_bed.bed
#Locating thos zones that are not in the bed file and masking them.
echo "bedtools complement and maskfasta."
bedtools complement -i $tmp/Sorted_bed.bed -g $ref.fai | bedtools maskfasta -fi $ref -bed - -fo $out/Masked_gen.fasta
#Indexing the newly masked genome.
echo "Indexing"
bwa index -a bwtsw $out/Masked_gen.fasta
samtools faidx $out/Masked_gen.fasta
echo "Dictionary file: "
picard CreateSequenceDictionary REFERENCE=$out/Masked_gen.fasta OUTPUT=$out/Masked_gen.dict
echo "Deleting Temporal file"
rm -r $tmp
