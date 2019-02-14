#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
version:3.0
by:lss
"""

import numpy as np
import pandas as pd
import sys
reload(sys)
sys.setdefaultencoding("utf8")
import os
import re

sys.path.append("/DATA/sslyu/trio_BWA-GATK_3.0/")
from lib.prepare import *
from lib.GAF import *
from lib.xlsx import *


def generate_cmd_part1(sample,ID,dir,md5,R1,R2,strategy):
    single_cmd = open(sample+'/cmd.sh', 'w')
    single_cmd.write('#!/bin/sh\n\n#SBATCH -J '+sample+'\n#SBATCH -p BP10\n#SBATCH -N 1\n#SBATCH -n 4\n\n')
    single_cmd.write('sample='+sample+'\nID='+ID+'\n\ndir='+dir+'\nmd5='+md5+'\nR1='+R1+'\nR2='+R2+'\nstrategy='+strategy+'\n\n')

    single_cmd.write('mkdir -p /anaData/anaData004/children_hos_genetic/rawdata/2018/$strategy/$ID\n')
    single_cmd.write('awk \'{if(match($2,"\'$ID\'")){print $0}}\' $md5 >> /anaData/anaData004/children_hos_genetic/rawdata/2018/$strategy/$ID/$ID.md5\n')
    single_cmd.write('cp $dir/*$ID*gz /anaData/anaData004/children_hos_genetic/rawdata/2018/$strategy/$ID\n')
    single_cmd.write('cp $dir/*$ID*$R1 ../raw/$sample\_R1.fastq.gz\n')
    single_cmd.write('cp $dir/*$ID*$R2 ../raw/$sample\_R2.fastq.gz\n')

    single_cmd.write('for j in $R1 $R2; do md5sum $dir/*$ID*$j >> ../original_md5; done\n')
    single_cmd.write('for j in R1 R2; do md5sum ../raw/$sample\_$j.fastq.gz >> ../cp_md5; done\n')
    single_cmd.write('for j in $R1 $R2; do md5sum /anaData/anaData004/children_hos_genetic/rawdata/2018/$strategy/$ID/*$ID*$j >> ../backup_md5; done\n')
    single_cmd.write('for j in $R1 $R2; do a=`grep $ID ../original_md5 |grep $j |cut -d" " -f 1`; b=`grep $ID $md5 |grep $j |cut -d" " -f 1`; if [ "$a" != "$b" ]; then echo -e "original "$ID"_"$j" is problematic" >> ../check_md5_log; echo -e "original "$ID"_"$j" is problematic"; fi; done\n')
    single_cmd.write('pattern_original=($R1 $R2); pattern_copied=(R1 R2); for j in 0 1; do original=`grep $ID ../original_md5 |grep ${pattern_original[$j]} |cut -d" " -f 1`; copied=`grep $sample ../cp_md5 |grep ${pattern_copied[$j]} |cut -d" " -f 1`; if [ "$original" != "$copied" ]; then echo -e $ID"_"${pattern_original[$j]}" is not copied correctly" >> ../check_md5_log; echo -e $ID"_"${pattern_original[$j]}" is not copied correctly"; fi; done\n')
    single_cmd.write('for j in $R1 $R2; do original=`grep $ID ../original_md5 |grep $j |cut -d" " -f 1`; backup=`grep $ID ../backup_md5 |grep $j |cut -d" " -f 1`; if [ "$original" != "$backup" ]; then echo -e $ID"_"$j" is not backup correctly" >> ../check_md5_log; echo -e $ID"_"$j" is not backup correctly"; fi; done\n')
    single_cmd.write('sh run2_sentieon.sh $sample $SLURM_NPROCS\n')
    single_cmd.close()

def generate_cmd_part2(sample):
    single_cmd = open(sample+'/cmd.sh', 'a')
    single_cmd.write('cp 2_mapping/$sample.depth.sample_gene_summary ../annotation\n')
    single_cmd.write('cp 3_variants/$sample.{ann.hg19_multianno.txt,CADD,maf,link} ../annotation\n')
    single_cmd.write('cd ../annotation\n')
    single_cmd.write('echo $sample >> list\n')
    single_cmd.write('python /DATA/sslyu/trio_BWA-GATK_3.0/src/score_re.py -sn $sample\n')
    single_cmd.write('python /DATA/sslyu/trio_BWA-GATK_3.0/src/annotation_filt_ver3.0.py -sn $sample -c_f_m c\n')
    single_cmd.write('echo "end `date`" >> ../$sample/finished\n')
    ofile = open('annotation/annotation.sh', 'a')
    ofile.write('python /DATA/sslyu/trio_BWA-GATK_3.0/src/annotation_filt_ver3.0.py -sn '+sample+' -c_f_m c\n')
    ofile.close()
    single_cmd.close()

def generate_trio_cmd(k, filtmode, pedigree_dict):              
    sample = d1.loc[k]['sample']
    path = 'trio/'+sample
    if not os.path.exists(path):
        os.makedirs(path)
    trio_cmd = open(path+'/trio_cmd.sh', 'w')
    trio_cmd.write('#!/bin/sh\n\n#SBATCH -J '+sample+'\n#SBATCH -p BP10\n#SBATCH -N 1\n#SBATCH -n 4\n\nsample='+sample+'\n\n')
    trio_cmd.write('filtmode=\''+filtmode+'\'\npath=`pwd`\nvar=`echo $path | awk -F \'/\' \'{print $NF}\'`\n')
    trio_cmd.write('if [ $var = \'peddy\' ]\nthen\nsh /DATA/sslyu/trio_BWA-GATK_3.0/src/trio.sh $sample $SLURM_NPROCS\n')
    trio_cmd.write('else\nsh /DATA/sslyu/trio_BWA-GATK_3.0/src/trio.sh $sample $SLURM_NPROCS $filtmode\n')
    trio_cmd.write('[ -e ../../$sample/2_mapping/$sample.depth.sample_gene_summary ] && cp ../../$sample/2_mapping/$sample.depth.sample_gene_summary ../../annotation\ncp segtrio/$sample.{ann.hg19_multianno.txt,CADD,maf,link} ../../annotation\ncd ../../annotation\necho $sample >> list\npython /DATA/sslyu/trio_BWA-GATK_3.0/src/score_re.py -sn $sample\npython /DATA/sslyu/trio_BWA-GATK_3.0/src/annotation_filt_ver3.0.py -sn $sample -c_f_m c_f_m\necho "end `date`" >> ../trio/$sample/finished\nfi\n')
    trio_cmd.write('[ -e ../../$sample/2_mapping/$sample.depth.sample_gene_summary ] && cp ../../$sample/2_mapping/$sample.depth.sample_gene_summary ../../annotation\n')
    trio_cmd.write('cp segtrio/$sample.{ann.hg19_multianno.txt,CADD,maf,link} ../../annotation\n')
    trio_cmd.write('cd ../../annotation\n')
    trio_cmd.write('echo $sample >> list\n')
    trio_cmd.write('python /DATA/sslyu/trio_BWA-GATK_3.0/src/score_re.py -sn $sample\n')
    relation = []
    for j in pedigree_dict[k]:
        relation.append(d1.loc[k]['relationship'])
    a = ' '.join(relation)
    annotationfile = open('annotation/annotation.sh', 'a')
    if u'父' in a and u'母' in a:
        trio_cmd.write('python /DATA/sslyu/trio_BWA-GATK_3.0/src/annotation_filt_ver3.0.py -sn $sample -c_f_m c_f_m\n')
        annotationfile.write('python /DATA/sslyu/trio_BWA-GATK_3.0/src/annotation_filt_ver3.0.py -sn '+sample+' -c_f_m c_f_m\n')
    else:
        if u'父' in a:
            trio_cmd.write('python /DATA/sslyu/trio_BWA-GATK_3.0/src/annotation_filt_ver3.0.py -sn $sample -c_f_m c_f\n')
            annotationfile.write('python /DATA/sslyu/trio_BWA-GATK_3.0/src/annotation_filt_ver3.0.py -sn '+sample+' -c_f_m c_f\n')
        else:
            trio_cmd.write('python /DATA/sslyu/trio_BWA-GATK_3.0/src/annotation_filt_ver3.0.py -sn $sample -c_f_m c_m\n')
            annotationfile.write('python /DATA/sslyu/trio_BWA-GATK_3.0/src/annotation_filt_ver3.0.py -sn '+sample+' -c_f_m c_m\n')
    trio_cmd.write('echo "end `date`" >> ../trio/$sample/finished\n')
    trio_cmd.write('fi\n')
    trio_cmd.close()
    annotationfile.close()

print """
Examples: 
prepare config
-sn excelfile 
-strategy IDT-PANEL 
-filtmode hard 
-dir absolute_path_of_rawdata 
-md5 absolute_path_of_rawdata_of_md5_file 
-R1 suffix_of_raw_R1 
-R2 suffix_of_raw_R2
Note:
make sure excel colnames contain 样本编号 reads, bases(Mb), Q30 and they are not null; the rawdata and md5 fileare in the same single dir 
"""
    
kwargs={'-sn':'', '-strategy':'', '-filtmode':'', '-dir':'', '-md5':'', '-R1':'', '-R2':''}

with open('config') as f:
    for l in f:
        l = l.strip('\n').split()
        for j in kwargs.keys():
		    if l[0] == j:
	    		kwargs.update({j:l[1]})
        
sn=kwargs.get('-sn')
strategy=kwargs.get('-strategy')
filtmode=kwargs.get('-filtmode')
dir = kwargs.get('-dir')
md5 = kwargs.get('-md5')
R1 = kwargs.get('-R1')
R2 = kwargs.get('-R2')

generate_dbevn(strategy)
with open('dbevn.sh') as f:
    for l in f:
        if l.startswith('panel'):
            print l.strip('\n')

print '-filtmode: '+filtmode

report_dirname = os.getcwd().split('/')[-1]
for i in ['ped', 'peddy', 'gvcf', 'annotation', 'raw', 'CNV', report_dirname]:
    if not os.path.exists(i):
        os.mkdir(i)

excel = pd.read_excel(sn)
excel = append_sample_txt(excel)
print excel[[u'姓名', u'性别', 'sample']]
d1 = append_gender(excel)
d1 = append_relation(d1)
d1 = append_pedigree(d1)
d1 = append_phenotype(d1)
d1 = append_father_mother(d1)
print d1

##### step3: generate cmd files
pedigree_dict = generate_pedigree(d1)
trio = pedigree_dict.keys()

for i in d1['sample']:
	if not os.path.exists(str(i)):
		os.mkdir(str(i))

run2_sentieon = open('/DATA/sslyu/trio_BWA-GATK_3.0/src/run2_sentieon.sh').readlines()
for i in range(21,31):
    print run2_sentieon[i]

single = generate_single(d1)
if len(single) > 0:
    list_single = open('list_single', 'w')
    if strategy == "WES":
        cmd_single = "clean mapping precalling gvcf index_single gtgvcf vqsr left-normalize depth qcsum qcvcf annotation"
    else:
        cmd_single = "clean mapping precalling gvcf index_single gtgvcf hardfilt left-normalize depth qcsum qcvcf annotation"
    print d1
    for i in single:
        d1.index = d1[u'姓名']
        sample = d1.loc[i]['sample']
        print 'single '+sample
        ID = d1.loc[i][u'样本编号']
        list_single.write(sample+'\n')
        generate_cmd_part1(sample, ID, dir, md5, R1, R2, strategy)
        generate_cmd_part2(sample)
        generate_sinle_ped(sample,d1)
        generate_run2_sentieon(sample,cmd_single)
    list_single.close()

if len(trio) > 0:
    if not os.path.exists('trio'):
        os.mkdir('trio')
    cmd_trio = "clean mapping precalling gvcf index_trio depth qcsum"
    list_trio = open('list_trio', 'w')
    ofile = open('trio/list', 'w')
    for k in pedigree_dict:
        d1.index = d1[u'姓名']
        generate_ped(k,d1)
        d1.index = d1[u'姓名']
        sample = d1.loc[k]['sample']
        ofile.write(sample+'\n')
        generate_trio_cmd(k, filtmode, pedigree_dict)

        for i in pedigree_dict[k]:
            d1.index = d1[u'姓名']
            sample = d1.loc[i]['sample']
            print 'trio '+sample
            ID = d1.loc[k][u'样本编号']
            list_trio.write(sample+'\n')
            generate_cmd_part1(sample, ID, dir, md5, R1, R2, strategy)
            generate_run2_sentieon(sample,cmd_trio)
    list_trio.close()
    ofile.close()

###### step5: generate related files and dirs
d1['sample'].to_csv('list', header=None, index=None)
d1['sample'].to_csv('CNV/id', header=None, index=None)
d1['CNV'] = '/BP12_share/sslyu/bam/'+d1['sample']+'.sort.mkdup.bam'
d1['CNV'].to_csv('CNV/list', header=None, index=None)

if os.path.exists('trio'):
    os.system('sort trio/list |uniq > tmp; mv tmp trio/list')
os.system('cp /DATA/sslyu/trio_BWA-GATK_3.0/src/trio_cmd.sh peddy')
if os.path.exists('list_single'):
    os.system('for i in `cat list_single`; do cat ped/$i.mendel.ped >> peddy_tmp; done')
if os.path.exists('trio/list'):
    os.system('for i in `cat trio/list`; do cat ped/$i.mendel.ped >> peddy_tmp; done')
os.system('sort peddy_tmp |uniq |awk \'{OFS="\t"; print \'0\',$2,$3,$4,$5,$6}\' > ped/peddy.ped; rm peddy_tmp')
with open('ped/peddy.ped') as f:
    for l in f:
        print l.strip('\n')
os.system('sort annotation/annotation.sh |uniq > annotation/tmp; mv annotation/tmp annotation/annotation.sh')        
#####################################################################################
#import shutil
#shutil.copyfile('/DATA/sslyu/trio_BWA-GATK_2.7/sub.sh', 'sub.sh')

