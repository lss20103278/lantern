#!/bin/sh

#SBATCH -J jobname
#SBATCH -p BP10
#SBATCH -N 1
#SBATCH -n 4

sample=list

sh run2_sentieon.sh $sample $SLURM_NPROCS

cp 2_mapping/$sample.depth.sample_gene_summary ../annotation
cp 3_variants/$sample.{ann.hg19_multianno.txt,CADD,maf,link} ../annotation
cd ../annotation
echo $sample >> list
python /DATA/sslyu/trio_BWA-GATK_3.0/src/score_re.py -sn $sample
python /DATA/sslyu/trio_BWA-GATK_3.0/src/annotation_filt_ver3.0.py -sn $sample -c_f_m c
echo "end `date`" >> ../$sample/finished

