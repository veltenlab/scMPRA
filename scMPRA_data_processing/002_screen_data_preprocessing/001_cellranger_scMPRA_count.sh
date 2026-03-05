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

dir=/users/lvelten/project/SCG4SYN/Experiments/240315_scMPRA_LibAlpha_Wilkinson/screen_HSC/001_cellranger

export PATH=/users/lvelten/sbeneyto/software/cellranger-6.0.1:$PATH

cd "$dir"

#Run CellRanger
cellranger count --id=scMPRA9_filtered_bcs --libraries=library_scMPRA9_TAP.csv --transcriptome=/nfs/users/lvelten/project/SCG4SYN/Experiments/231127_scMPRA_LibAlpha_HSC/screen_HSC/001_cellranger/reference/refdata-gex-mm10-2020-A --feature-ref=feature_ref_filtered.csv

