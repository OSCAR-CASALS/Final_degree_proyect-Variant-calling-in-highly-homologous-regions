#!/bin/bash

gff_file=$1
ref=$2
Res_dir=$3

awk -F '\t' '{if($3=="pseudogene") print $1"\t"$4"\t"$5}' $gff_file > $Res_dir/Pseudogenes.bed
awk -F '\t' '{if($3=="gene") print $1"\t"$4"\t"$5}' $gff_file > $Res_dir/genes.bed

bedtools getfasta -fi $ref -bed $Res_dir/Pseudogenes.bed > $Res_dir/Pseudogenes.fasta
bedtools getfasta -fi $ref -bed $Res_dir/genes.bed > $Res_dir/genes.fasta

rm $Res_dir/Pseudogenes.bed
rm $Res_dir/genes.bed

mkdir $Res_dir/Blast_Database

makeblastdb -in $Res_dir/genes.fasta -out $Res_dir/Blast_Database/GENE_DATABASE -dbtype nucl

rm $Res_dir/genes.fasta

