#!/bin/sh

#if [ ! -e ped/peddy.ped ]
#then
#echo "ERROR:peddy.ped must exist,please provide.Exiting....."
#exit
#fi
#
#gname=`awk '{print $2}' ped/peddy.ped`
#for i in $gname
#do
#[ ! -e gvcf/$i.raw.g.vcf ] && mv gvcf/$i.g.vcf gvcf/$i.raw.g.vcf
#./gvcf.sh $i
#done

#hardfilt=$1 

for i in `cat trio/list`
do
[ ! -e trio/$i ] && mkdir -p trio/$i
#cp /DATA/sslyu/trio_BWA-GATK_3.1/src/{trio_cmd.sh,vfilt.sh,vep-mendelscan.sh,trio.sh,sort_sample.py,mendel_to_annovar.py} trio/$i
cd trio/$i 
sed -i "s/^sample.*$/sample=$i/;s/^#SBATCH -J .*$/#SBATCH -J $i/" trio_cmd.sh
#sed -i "s/^filtmode.*$/filtmode='$hardfilt'/" trio_cmd.sh
sbatch trio_cmd.sh
cd ../..
done
