#########################################################################
# File Name: sub_fqstat.sh
# Author: sslyu
# mail: 1573077420@qq.com
# Created Time:Thu 11 Oct 2018 03:08:02 PM CST
#########################################################################
#!/bin/bash

for i in `cat list`; do cp fqstat.sh $i; cd $i; sed -i "s/^sample=.*$/sample=$sample/;s/^#SBATCH -J .*$/#SBATCH -J $sample/" fqstat.sh; sbatch fqstat.sh; done

########################################################################################################
#sn=$1
#a=`cat $sn`
#sample=($a)
#
##rawdata=$2
##b=`cat $rawdata`
##rawdata_path=($b)
#
#num=${#sample[@]}
#
#for i in `seq 0 $[$num-1]`
#do
#	#cp fqstat.sh ${sample[$i]}\_fqstat.sh
#	#c=`cat ${rawdata_path[$i]}`
#	#sample_rawdata_path=($c)
#	#R1=${sample_rawdata_path[0]}
#	#R2=${sample_rawdata_path[1]}
#	#sed -i "s/^sample.*$/sample=${sample[$i]}/" ${sample[$i]}\_fqstat.sh
#	#sed -i "s/^R1.*$/R1=$R1/" ${sample[$i]}\_fqstat.sh # there is / in R1, so the sed interpreter can't interprete $R1
#	#sed -i "s/^R2.*$/R2=$R2/" ${sample[$i]}\_fqstat.sh
#	sbatch ${sample[$i]}\_fqstat.sh
#done
