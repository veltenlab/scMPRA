#!/bin/bash
#$ -N map
#$ -cwd 
#$ -e mapiseq_BIG.err
#$ -o mapiseq_BIG.out
#$ -l virtual_free=20G,h_rt=24:00:00
#$ -q long-centos79

#For Library H the memory footprint was 58Gig and it took 13 hours. The memory is larger, if there are many barcodes.

module use /software/as/el7.2/EasyBuild/CRG/modules/all
module load BWA/0.7.17
module load SAMtools/1.10-GCC-9.3.0

FASTQDIR=/PATHTOFASTQ/

GUIDE="$FASTQDIR/JRBETA_lib_03750AAF_TGGTCTACGT_R1_001.fastq.gz";
CRS="$FASTQDIR/JRBETA_lib_03750AAF_TGGTCTACGT_R3_001.fastq.gz";
BC="$FASTQDIR/JRBETA_lib_03750AAF_TGGTCTACGT_R2_001.fastq.gz";

OUT="/OUTPATH/"

REF=/PATHTOREF/example_fasta_CRE_reference.fa

THREADS=1

#__________________STEP 1______________________#
#creates index for your CRE reference 

bwa index $REF

#__________________STEP 2______________________#
#you are creating a .bam file with aligned reads against the CRE reference 

bwa mem -t $THREADS $REF $CRS | samtools view -b > $OUT/aligned_BIG.bam 2> $OUT/alignment_BIG.log 

#__________________STEP 3______________________#
#here you are running the custom .pl script, that counts reads and matches the barcodes to create the look-up table  

perl ./002_map_crs_bc_BAM_2_singlecell.pl  $BC $GUIDE $OUT/aligned_BIG.bam $OUT/mapped_BIG.csv.gz 2> $OUT/mapping_BIG.log

#__________________STEP 4______________________#
# run 003_barcode_assoc_report.Rmd to perform QC of reads 



#this is a useful command for finding all the alignments of a given barcode
#zcat $BC | grep -B1 GTAACCACAGAGTTG | grep "M03766" | perl -pe 's/^@(\S+)\s.+/$1/' > test.txt
#samtools view $OUT/aligned.bam | grep -f test.txt 
