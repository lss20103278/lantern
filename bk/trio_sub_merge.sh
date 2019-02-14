#!/bin/sh

if [ ! -e ped/peddy.ped ]
then
echo "ERROR:peddy.ped must exist,please provide.Exiting....."
exit
fi

#merge:
gname=`awk '{print $2}' ped/peddy.ped`
for i in $gname
do
if [ ! -e gvcf/$i.sorted.vcf.gz ]	
then	
(grep ^# ../$i.vcf; grep -v ^# ../$i.vcf |sort -k1,1 -k2,2n) |bgzip > gvcf/$i.sorted.vcf.gz
tabix -p vcf gvcf/$i.sorted.vcf.gz
fi
done

for i in `cat list`
do
[ ! -e $i ] && mkdir -p $i
cd $i 
#cp ../*.sh ./
if [ -e cmd_merge.sh ]
then	
sed -i "s/^sample.*$/sample=$i/;s/^#SBATCH -J .*$/#SBATCH -J $i/" cmd_merge.sh
sbatch cmd_merge.sh
else
sed -i "s/^sample.*$/sample=$i/;s/^#SBATCH -J .*$/#SBATCH -J $i/" cmd.sh
sbatch cmd.sh
fi
cd ..
done
