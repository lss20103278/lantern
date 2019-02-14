#########################################################################
# File Name: PACS2_statistics.sh
# Author: sslyu
# mail: 1573077420@qq.com
# Created Time:Sat 29 Sep 2018 02:26:57 PM CST
# Usage: This script is for calculating the frequency of one gene in all existing samples
#########################################################################
#!/bin/bash

gene=$1

grep $gene /DATA/sslyu/soft/annovar/humandb/refseq.sorted.txt |awk '{l=split($10,a,",");m=split($11,b,",");for(i=1;i<l;i++){print a[i],b[i]}}' |sort |uniq > $gene\_all_uniq.bed
#split -l 1 -d --additional-suffix .bed $gene\_all_uniq.bed $gene\_

chrom=`grep $gene /DATA/sslyu/soft/annovar/humandb/refseq.sorted.txt |cut -f 3|sort |uniq`
a=`cut -f 1 $gene\_all_uniq.bed`
start=($a)
b=`cut -f 2 $gene\_all_uniq.bed`
end=($b)
num=${#start[@]}
mkdir children_hospital_all_vcf_$gene
for i in `find /DATA/BPshare/opm/遗传病-交大附属儿医/ -type f -name "*.vcf"`; do for j in `seq 0 $[$num-1]`; do awk '{if($2>="'${start[$j]}'" && $2<="'${end[$j]}'" && $1=="'$chrom'"){print $0}}' $i >> children_hospital_all_vcf_$gene/$gene\_${i##*/}; done; done
find children_hospital_all_vcf_$gene/ -name "*" -type f -size 0c | xargs -n 1 rm -f

for i in `ls children_hospital_all_vcf_$gene/`; do a=${i#*_}; awk '{OFS="\t"; print $0,"'${a%.*}'"}' children_hospital_all_vcf_$gene/$i >> $gene\_tmp_sample; done
awk '{l=split($10,a,":");$10=a[1];if($10 != "0/0"){print $0}}' $gene\_tmp > $gene\_for_annovar
/DATA/sslyu/soft/annovar/convert2annovar.pl -format vcf4 $gene\_for_annovar -outfile $gene.avinput
/DATA/sslyu/soft/annovar/table_annovar.pl $gene.avinput /DATA/sslyu/soft/annovar/humandb -buildver hg19 -out PACS2.ann -protocol refGene,gnomad_genome_eas,gnomad_exomes_hom,gnomad_exomes_hom_all,gnomad_exomes_hemi,clinvar_20180603,hgmd_2018_spring -operation g,f,f,f,f,f,f -arg ',,,,,,' -remove -nastring "NA"


#sum=`find /DATA/BPshare/opm/遗传病-交大附属儿医/ -type f -name "*.vcf" |wc -l`
#PACS2=`cut -f 11 PACS2_tmp |sort |uniq |wc -l`

