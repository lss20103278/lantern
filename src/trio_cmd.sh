#!/bin/sh

#SBATCH -J jobname
#SBATCH -p BP10
#SBATCH -N 1
#SBATCH -n 4

sample=list


#select from 'hard'/'vqsr'
filtmode='vqsr'

path=`pwd`
var=`echo $path | awk -F '/' '{print $NF}'`
#echo $path
#if echo $path | grep peddy > /dev/null
if [ $var = 'peddy' ]
then
sh /DATA/sslyu/trio_BWA-GATK_3.0/src/trio.sh $sample $SLURM_NPROCS
else
sh /DATA/sslyu/trio_BWA-GATK_3.0/src/trio.sh $sample $SLURM_NPROCS $filtmode
[ -e ../../$sample/2_mapping/$sample.depth.sample_gene_summary ] && cp ../../$sample/2_mapping/$sample.depth.sample_gene_summary ../../annotation
cp segtrio/$sample.{ann.hg19_multianno.txt,CADD,maf,link} ../../annotation
cd ../../annotation
echo $sample >> list
python /DATA/sslyu/trio_BWA-GATK_3.0/src/score_re.py -sn $sample
python /DATA/sslyu/trio_BWA-GATK_3.0/src/annotation_filt_ver3.0.py -sn $sample -c_f_m c_f_m
echo "end `date`" >> ../trio/$sample/finished
fi
