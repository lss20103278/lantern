#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import pandas as pd

ofile = open('step5_backup.sh', 'w')
dirname = os.getcwd().split('/')[-2:]

kwargs_raw = sys.argv[1:]
kwargs={'-sn':'', '-l':''}
for i in range(len(kwargs_raw)):
    for j in kwargs.keys():
        if(kwargs_raw[i]==j):
            kwargs.update({j:kwargs_raw[i+1]})

#read sample list
if(kwargs.get('-sn')!=''):
    list=[]
    list.append(kwargs.get('-sn'))
    print "\nNote:sigle sample mode"
else:
    if kwargs.get('-l')!='':
        list=pd.read_table(kwargs.get('-l'),header=None)[0].tolist()
        print "\nNote:list samples mode"

for sn in list:
    strategy = os.getcwd().split('/')[-1].split('_')[-1]
    print strategy
    print '\n'
    ofile.write('sh /DATA/sslyu/trio_BWA-GATK_3.0/src/unmapped_duplicate.sh list\n')
    ofile.write('mkdir -p '+dirname[-1]+'/vcf\n')
    if os.path.exists('trio'):
        ofile.write('for i in `cat trio/list`; do cp trio/$i/sep/*vcf '+dirname[-1]+'/vcf; done\n')
        ofile.write('sort annotation/list |uniq > annotation/tmp; mv annotation/tmp annotation/list\n')
        ofile.write('for i in `cat annotation/list trio/list |sort |uniq -u`; do cp $i/3_variants/$i.vcf '+dirname[-1]+'/vcf; done\n')
    else:
        ofile.write('for i in `cat list_single`; do cp $i/3_variants/$i.vcf '+dirname[-1]+'/vcf; done\n')
    ofile.write('cp '+sn+' '+dirname[-1]+'\n')
    ofile.write('cp annotation/*_children_hospital_v3.0.xlsx '+dirname[-1]+'\n')
    ofile.write('cd CNV/\n')
    ofile.write('cat /DATA/sslyu/trio_BWA-GATK_3.0/src/CNV_Read.R\n')
    ofile.write('/DATA/ypliu/opt/R-3.4.3/bin/R # note: choose the right bed file\n\n')
    ofile.write('cd ..\n')
    ofile.write('mkdir -p '+dirname[-1]+'/CNV; mv CNV/*csv '+dirname[-1]+'/CNV\n')
    ofile.write('zip -r '+dirname[-1]+'.zip '+dirname[-1]+'\n')
    ofile.write('ls -R '+dirname[-1]+'\n\n')
    ofile.write('if [ ! -e /DATA/BPshare/opm/遗传病-交大附属儿医/2018/'+dirname[-2]+' ]\n')
    ofile.write('then\n')
    ofile.write('mkdir -p /DATA/BPshare/opm/遗传病-交大附属儿医/2018/'+dirname[-2]+'\n')
    ofile.write('fi\n')
    ofile.write('cp -r '+dirname[-1]+' /DATA/BPshare/opm/遗传病-交大附属儿医/2018/'+dirname[-2]+'\n')
    ofile.write('for i in `cat list`; do scp /BP12_share/sslyu/bam/$i.sort.mkdup.bam* klyang@122.112.248.194:/media/bpdata/data/exon_temp2/; done\n')
    ofile.write('scp list klyang@122.112.248.194:/media/bpdata/data/exon_temp2/list_'+dirname[-1]+'\n')
    ofile.write('ssh kly\n\n')
    ofile.close()
    
    with open('step5_backup.sh') as f:
    	for l in f:
    		print l.strip('\n')
    #truth = raw_input('Is the shell script right? If not, please enter no: ')
    #if truth != '':
    #	sys.exit()
    #else:
    #	os.system('bash step5_backup.sh')
