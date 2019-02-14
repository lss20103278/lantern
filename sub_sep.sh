#!/bin/sh

sample=$1
[ ! -e $sample ] && mkdir -p $sample
cd $sample
sed -i "s/^sample.*$/sample=$sample/;s/^#SBATCH -J .*$/#SBATCH -J $sample/" cmd.sh
sbatch cmd.sh
cd ..


#for i in `cat list`
#do
#[ ! -e $i ] && mkdir -p $i
#cd $i 
##cp ../*.sh ./
##cp ../dbevn.sh ./
#sed -i "s/^sample.*$/sample=$i/;s/^#SBATCH -J .*$/#SBATCH -J $i/" cmd.sh
#sbatch cmd.sh
#cd ..
#done
