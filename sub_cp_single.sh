#########################################################################
# File Name: sub_cp_single.sh
# Author: sslyu
# mail: 1573077420@qq.com
# Created Time:Wed 19 Sep 2018 11:03:04 PM CST
#########################################################################
#!/bin/bash

sample=$1
ID=$2

cp /DATA/sslyu/trio_BWA-GATK_3.1/src/cp_rawdata.sh $sample
cd $sample
sed -i "s/^sample.*$/sample=$sample/;s/^#SBATCH -J .*$/#SBATCH -J $sample/" cp_rawdata.sh
sed -i "s/^ID.*$/ID=$sample/" cp_rawdata.sh
#sbatch cp_rawdata.sh
cd ..
