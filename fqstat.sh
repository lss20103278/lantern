#!/bin/bash

#SBATCH -J jobname
#SBATCH -p BP10
#SBATCH -N 1
#SBATCH -n 4

sample=list

/home/ana005/anaconda2/bin/iTools Fqtools stat -InFq ../raw/$sample\_R1.fastq.gz -InFq ../raw/$sample\_R2.fastq.gz -OutStat $sample.info
