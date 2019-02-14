#!/bin/sh

list=$1
for i in `cat $list`
do
[ ! -e $i ] && mkdir -p $i
cd $i 
#cp ../*.sh ./
#cp ../dbevn.sh ./
sed -i "s/^sample.*$/sample=$i/;s/^#SBATCH -J .*$/#SBATCH -J $i/" cmd_step2.sh
sbatch cmd_step2.sh
cd ..
done
