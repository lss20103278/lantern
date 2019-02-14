#!/bin/sh

#SBATCH -J jobname
#SBATCH -p BP10
#SBATCH -N 1
#SBATCH -n 4

sample=list
ID=list

dir=dir
md5=md5
R1=R1
R2=R2
strategy=strategy

mkdir -p /anaData/anaData004/children_hos_genetic/rawdata/2018/$strategy/$ID
awk '{if(match($2,"'$ID'")){print $0}}' $dir/$md5 >> /anaData/anaData004/children_hos_genetic/rawdata/2018/$strategy/$ID/$ID.md5
cp $dir/*$ID*gz /anaData/anaData004/children_hos_genetic/rawdata/2018/$strategy/$ID
cp $dir/*$ID*$R1 ../raw/$sample\_R1.fastq.gz
cp $dir/*$ID*$R2 ../raw/$sample\_R2.fastq.gz


