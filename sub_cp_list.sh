#########################################################################
# File Name: sub_cp.sh
# Author: sslyu
# mail: 1573077420@qq.com
# Created Time:Wed 19 Sep 2018 11:00:45 PM CST
#########################################################################
#!/bin/bash

sn=$1
dir=$2
md5=$3
R1=$4
R2=$5
strategy=$6

for i in `cat $sn |tr "\t" "-"`
do
	ID=${i%-*}; sample=${i#*-}
	cp /DATA/sslyu/trio_BWA-GATK_3.0/src/cp_rawdata.sh $sample
	cd $sample
	sed -i "s/^sample.*$/sample=$sample/;s/^#SBATCH -J .*$/#SBATCH -J $sample/;s/^ID.*$/ID=$ID/" cp_rawdata.sh
	#sed -i "s/^dir=.*$/dir=$dir/;s/^md5=.*$/md5=$md5/;s/^R1=.*$/R1=$R1/;s/^R2=.*$/R2=$R2/;s/^strategy=.*$/strategy=$strategy/" cp_rawdata.sh
	sed -i "s/^dir.*$/dir=$dir/;s/^md5.*$/md5=$md5/;s/^R1.*$/R1=$R1/;s/^R2.*$/R2=$R2/;s/^strategy.*$/strategy=$strategy/" cp_rawdata.sh
	sbatch cp_rawdata.sh 
	cd ..
done	

####################################################################################################
#a=`cut -f 4 $sn |sed '1d'`
#a=`cut -f 2 $sn`
#sample=($a)
#if [ "${sample[0]}" == "nan" ]; then
#	a=`cut -f 3 $sn |sed '1d'`
#	sample=($a)
#fi
#b=`cut -f 3 $sn |sed '1d'`
#b=`cut -f 1 $sn`
#ID=($b)
#num=${#sample[@]}
#
#for i in `seq 0 $[$num-1]`
#do
#	[ ! -e ${sample[$i]} ] && mkdir -p ${sample[$i]}
#	cp cp_rawdata.sh ${sample[$i]}
#	cd ${sample[$i]}
#	sed -i "s/^sample.*$/sample=${sample[$i]}/;s/^#SBATCH -J .*$/#SBATCH -J ${sample[$i]}/" cp_rawdata.sh
#	sed -i "s/^ID.*$/ID=${ID[$i]}/" cp_rawdata.sh
#	sbatch cp_rawdata.sh $dir $md5 $R1 $R2 $strategy
#	cd ..
#done

