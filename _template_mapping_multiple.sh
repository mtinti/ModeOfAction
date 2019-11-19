#!/bin/bash -l
#$ -cwd
#$ -V
#$ -j y
#$ -pe smp 8
#$ -o /dev/null
#$ -e /dev/null


grid_out='/homes/mtinti/RNAseq/viper-test/rit-seq/{experiment}'

exec >$grid_out'/{base_fastq}.out.txt' 2>$grid_out'/{base_fastq}.err.txt'

echo 'redirect grid outputs to:' $grid_out

conda activate ritSeq

echo 'conda env loaded'

#for the future, if always one fastqfile * experiment
#g_version=$1
#experiment=$2
#base_fastq=$3
#echo $g_version $experiment $base_fastq


echo 'load variable'
##variables
path_genome_index='genomes/{g_version}/{g_version}'
path_fastq_files='{experiment}/data/'
base_fastq='{base_fastq}'
path_out=$path_fastq_files$base_fastq'/'
gtf_file='genomes/{g_version}/{g_version}.gtf'
gff_file='genomes/{g_version}/{g_version}.gff'

echo 'make out dirs'
mkdir $path_out
mkdir $path_out'fastqc'
##qc of fastq 
echo 'run 0.1' 


{remove_1}cat {reads_1} >  '{experiment}/data/{base_fastq}1.fq.gz'
{remove_1}cat {reads_2} >  '{experiment}/data/{base_fastq}2.fq.gz'

inFileLength=$(echo $(zcat $path_fastq_files$base_fastq'1.fq.gz'|wc -l)/4|bc)
echo 'inFileLength: '$inFileLength

echo 'run 0.2 parse fastq for barcodes' 
python parse_barcode2.py '{library}' $path_fastq_files$base_fastq $inFileLength

echo 'run 1.2 align ff barcode reads' 
fastq_1=$path_fastq_files$base_fastq'1_ff_barcode.fq.gz'
fastq_2=$path_fastq_files$base_fastq'2_ff_barcode.fq.gz'
bowtie2 --very-sensitive-local -p 8 -x $path_genome_index -1 $fastq_1 -2 $fastq_2 | samtools view -bS - | samtools sort -o $path_out$base_fastq'ff_barcode_sorted.bam'
bedtools genomecov -ibam $path_out$base_fastq'ff_barcode_sorted.bam' -bg > $path_out$base_fastq'ff_barcode_coverage_bg.bed'
bedtools genomecov -ibam $path_out$base_fastq'ff_barcode_sorted.bam'  -d > $path_out$base_fastq'ff_barcode_coverage_d.bed'

echo 'run 1.3 align fr barcode reads' 
fastq_1=$path_fastq_files$base_fastq'1_fr_barcode.fq.gz'
fastq_2=$path_fastq_files$base_fastq'2_fr_barcode.fq.gz'
bowtie2 --very-sensitive-local -p 8 -x $path_genome_index -1 $fastq_1 -2 $fastq_2 | samtools view -bS - | samtools sort -o $path_out$base_fastq'fr_barcode_sorted.bam'
bedtools genomecov -ibam $path_out$base_fastq'fr_barcode_sorted.bam' -bg > $path_out$base_fastq'fr_barcode_coverage_bg.bed'
bedtools genomecov -ibam $path_out$base_fastq'fr_barcode_sorted.bam'  -d > $path_out$base_fastq'fr_barcode_coverage_d.bed'

echo 'run 1.4 align rr barcode reads' 
fastq_1=$path_fastq_files$base_fastq'1_rr_barcode.fq.gz'
fastq_2=$path_fastq_files$base_fastq'2_rr_barcode.fq.gz'
bowtie2 --very-sensitive-local -p 8 -x $path_genome_index -1 $fastq_1 -2 $fastq_2 | samtools view -bS - | samtools sort -o $path_out$base_fastq'rr_barcode_sorted.bam'
bedtools genomecov -ibam $path_out$base_fastq'rr_barcode_sorted.bam' -bg > $path_out$base_fastq'rr_barcode_coverage_bg.bed'
bedtools genomecov -ibam $path_out$base_fastq'rr_barcode_sorted.bam'  -d > $path_out$base_fastq'rr_barcode_coverage_d.bed'

echo 'run 1.5 align rf barcode reads' 
fastq_1=$path_fastq_files$base_fastq'1_rf_barcode.fq.gz'
fastq_2=$path_fastq_files$base_fastq'2_rf_barcode.fq.gz'
bowtie2 --very-sensitive-local -p 8 -x $path_genome_index -1 $fastq_1 -2 $fastq_2 | samtools view -bS - | samtools sort -o $path_out$base_fastq'rf_barcode_sorted.bam'
bedtools genomecov -ibam $path_out$base_fastq'rf_barcode_sorted.bam' -bg > $path_out$base_fastq'rf_barcode_coverage_bg.bed'
bedtools genomecov -ibam $path_out$base_fastq'rf_barcode_sorted.bam'  -d > $path_out$base_fastq'rf_barcode_coverage_d.bed'

##align the fastq files to the transcript and convert to bam
 
fastq_1=$path_fastq_files$base_fastq'1.fq.gz'
fastq_2=$path_fastq_files$base_fastq'2.fq.gz'

echo 'run fastqc' 
fastqc  $fastq_1 -o $path_out'fastqc'
fastqc  $fastq_2 -o $path_out'fastqc'

echo 'run 2.0 align all reads'
bowtie2 --very-sensitive-local -p 8 -x $path_genome_index -1 $fastq_1 -2 $fastq_2 -S $path_out$base_fastq'.sam'

echo 'run 2.1'
samtools view -bS $path_out$base_fastq'.sam' > $path_out$base_fastq'.bam'
rm $path_out$base_fastq'.sam' ##saving space

echo 'run 4'
samtools sort $path_out$base_fastq'.bam' -o $path_out$base_fastq'.sorted.bam'

echo 'run 4'
samtools sort -n $path_out$base_fastq'.bam' -o $path_out$base_fastq'.sorted_name.bam'

echo 'run 5'
samtools index $path_out$base_fastq'.sorted.bam'

echo 'run 6'
rm $path_out$base_fastq'.bam'

echo 'run 7'
bedtools genomecov -ibam $path_out$base_fastq'.sorted.bam' -bg > $path_out$base_fastq'_coverage_bg.bed'
echo 'run 7.0'
bedtools genomecov -ibam $path_out$base_fastq'.sorted.bam' -d > $path_out$base_fastq'_coverage_d.bed'
echo 'run 7.1'
python plot_chr_coverage.py $path_out$base_fastq'_coverage_d.bed' '{experiment}/'$base_fastq'coverage_d.png'

#echo 'run 7.1'
#python mod_bed_graph.py $path_out$base_fastq'_coverage_bg.bed'
#echo 'run 7.2'
#bedmap --mean $gff_file'.bed' $path_out$base_fastq'_coverage_bg.bed' > $path_out$base_fastq'_coverage_reference_mean.bed'

unset DISPLAY

echo 'run 8'
mkdir $path_out'qualimap_results_rnaseq'
qualimap rnaseq -bam $path_out$base_fastq'.sorted.bam'  -gtf $gtf_file -outdir $path_out'qualimap_results_rnaseq' -pe

echo 'run 8'
mkdir $path_out'qualimap_results_bamqc'
qualimap bamqc  -bam $path_out$base_fastq'.sorted.bam' -outdir $path_out'qualimap_results_bamqc'

echo 'run 9'
mkdir $path_out'unmap/'
samtools view -u -f 12 -F 256 $path_out$base_fastq'.sorted.bam' > $path_out'unmap/'$base_fastq'.sorted_unmap.bam'

echo 'run 10'
samtools sort -n $path_out'unmap/'$base_fastq'.sorted_unmap.bam' -o $path_out'unmap/'$base_fastq'.sorted_unmap.bam'

echo 'run 11'
bamToFastq -i $path_out'unmap/'$base_fastq'.sorted_unmap.bam' -fq $path_out'unmap/'$base_fastq'.sorted_unmap_1.fastq' -fq2 $path_out'unmap/'$base_fastq'.sorted_unmap_2.fastq'


echo 'run 12'

script_name='count_'$base_fastq'.R'
echo 'running' $script_name
Rscript --vanilla $script_name

echo 'run 13'
mv $script_name {experiment}

echo 'run 14'
mkdir '{experiment}/data/{base_fastq}/fastq_screen'
fastq_screen --nohits $path_out'unmap/'$base_fastq'.sorted_unmap_1.fastq' $path_out'unmap/'$base_fastq'.sorted_unmap_2.fastq' \
--conf 'FastQ_Screen_Genomes/fastq_screen.conf' --outdir '{experiment}/data/{base_fastq}/fastq_screen'

rm $path_out'unmap/'$base_fastq'.sorted_unmap_1.fastq'
rm $path_out'unmap/'$base_fastq'.sorted_unmap_2.fastq'

echo 'run 16'
mkdir '{experiment}/macs2'
macs2 callpeak -t $path_out$base_fastq'.sorted.bam' -f BAMPE --min-length 10000 --max-gap 5000 -q 0.01 --broad --outdir '{experiment}/macs2_'$base_fastq

echo 'run 15'
multiqc '{experiment}/data/{base_fastq}' '{experiment}/macs2' -o '{experiment}/multi_qc_'$base_fastq

python plot_reagion_coverage.py $gff_file '{experiment}/macs2_'$base_fastq'/NA_peaks.xls' \
$path_out$base_fastq'_coverage_d.bed' $path_out$base_fastq'ff_barcode_coverage_d.bed' \
$path_out$base_fastq'fr_barcode_coverage_d.bed' $path_out$base_fastq'rf_barcode_coverage_d.bed' \
$path_out$base_fastq'rr_barcode_coverage_d.bed' '{experiment}/'$base_fastq '{experiment}/'$base_fastq'count.txt'




