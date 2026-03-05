#!/bin/bash
#$ -N scMPRA7
#$ -cwd
#$ -q long-centos79
#$ -l virtual_free=80G,h_vmem=80G,h_rt=48:00:00
#$ -pe smp 8
#$ -o /users/lvelten/project/SCG4SYN/Experiments/240315_scMPRA_LibAlpha_Wilkinson/screen_HSC/001_cellranger/scMPRA9_filtered_bcs.out
#$ -e /users/lvelten/project/SCG4SYN/Experiments/240315_scMPRA_LibAlpha_Wilkinson/screen_HSC/001_cellranger/scMPRA9_filtered_bcs.err
#$ -m be

#Add EasyBuild configuration
#module use /software/as/el7.2/EasyBuild/CRG/modules/all

dir=/scMPRA_data_processing/002_screen_data_preprocessing/

export PATH=/path_to_cellranger/software/cellranger-6.0.1:$PATH

cd "$dir"

#Run CellRanger
cellranger count --id=scMPRA9_filtered_bcs --libraries=library_scMPRA9_TAP.csv --transcriptome=/path_to_reference_genome/refdata-gex-mm10-2020-A --feature-ref=example_feature_ref_filtered.csv

