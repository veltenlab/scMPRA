#!/bin/bash
#$ -N scMPRA
#$ -cwd
#$ -q long-centos79
#$ -l virtual_free=80G,h_vmem=80G,h_rt=48:00:00
#$ -pe smp 8
#$ -o ./scMPRA_filtered_bcs.out
#$ -e ./scMPRA_filtered_bcs.err
#$ -m be

dir=/scMPRA_data_processing/002_screen_data_preprocessing/

export PATH=/path_to_cellranger/software/cellranger-6.0.1:$PATH

cd "$dir"

#Run CellRanger
cellranger count --id=scMPRA_filtered_bcs --libraries=library_scMPRA.csv --transcriptome=/path_to_reference_genome/refdata-gex-mm10-2020-A --feature-ref=example_feature_ref_filtered.csv

