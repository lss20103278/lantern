#########################################################################
# File Name: prepare.sh
# Author: sslyu
# mail: 1573077420@qq.com
# Created Time:Wed 19 Sep 2018 11:00:45 PM CST
#########################################################################
#!/bin/bash

#sample=$1
#serial=$2
#strategy=$3

serial=$1
strategy=$2

mkdir -p $serial/$serial\_$strategy

cd $serial/$serial\_$strategy
mkdir -p raw ped gvcf peddy annotation $serial\_$strategy
if [ "$strategy" != "WES" ]; then
	mkdir CNV
fi
cp /DATA/sslyu/trio_BWA-GATK_3.0/src/cp_rawdata.sh .
cp /DATA/sslyu/trio_BWA-GATK_3.0/src/trio_cmd.sh peddy
sed -i "s/^sample.*$/sample=peddy/;s/^#SBATCH -J .*$/#SBATCH -J peddy/" peddy/trio_cmd.sh

###############################################################################################
#### modify the excel filename
#na=$(echo $sample |tr ")" "_")
#cp "$sample" $na

#python /DATA/sslyu/trio_BWA-GATK_3.1/step1_excel2tmp.py -sn $na -serial $serial -strategy $strategy
#mv $serial\_$strategy\_excel.tmp $serial/$serial\_$strategy

#cp -r /DATA/sslyu/trio_BWA-GATK_3.1/src/ .
