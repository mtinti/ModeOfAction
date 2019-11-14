#!/bin/bash -l
#$ -cwd
#$ -V
#$ -j y
#$ -pe smp 8
#$ -o /dev/null
#$ -e /dev/null


# expecting a path/to/ file.fasta
# ie genomes/ld39
# fasta file names same as path basename
path_to=$1

grid_out='/homes/mtinti/RNAseq/viper-test/rit-seq/'$path_to

exec >$grid_out'/make_index.out.txt' 2>$grid_out'/make_index.err.txt'

conda activate ritSeq

in_file=$path_to$(basename $path_to)'.fasta'
out_file=$path_to$(basename $path_to)
echo $in_file 
echo $out_file
bowtie2-build $in_file $out_file


in_file=$path_to$(basename $path_to)'.gff'
out_file=$path_to$(basename $path_to)'.gtf'
echo $in_file 
echo $out_file
gffread -E $in_file  -T -o $out_file

convert2bed -i gff < $in_file > $in_file'.bed'

grep 'artemis	gene\|EuPathDB	gene' $in_file'.bed' > $in_file'.gene_bed'

