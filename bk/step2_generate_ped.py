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
from lib.GAF import *
from lib.xlsx import *

if len(sys.argv[1:]) != 6:
    print '''
    examples:
    python /DATA/sslyu/trio_BWA-GATK_3.0/src/step2_generate_ped.py -sn filename -strategy (WES/IDT-PANEL...) -filtmode (vqsr/hard) 
    '''
    sys.exit()

kwargs_raw = sys.argv[1:]
kwargs={'-sn':'', '-l':'', '-strategy':'', '-filtmode':'vqsr'}
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

strategy = kwargs.get('-strategy')
filtmode = kwargs.get('-filtmode')
print '-strategy:'+strategy,'-filtmode:'+filtmode

if len(list) != 0:
    if not os.path.exists('dbevn.sh'):
        strategy = raw_input('strategy(WES PANEL Agilent IDT Agilent_wes brain blood_disease, if not exist, please enter none):')
        generate_dbevn(strategy)
    if not os.path.exists('annotation'):
        os.mkdir('annotation')
    for sn in list:
    	print sn
        excel = pd.read_excel(sn, dtype=str)
        excel = append_sample_excel(excel)
        d1 = append_gender(excel)
        print d1[[u'姓名', 'sample', 'gender']]
        d1 = append_relation(d1)
        d1 = append_pedigree(d1)
        d1 = append_phenotype(d1)
        d1.to_csv('tmp', sep="\t")
        d1 = append_father_mother(d1)

        pedigree_dict = generate_pedigree(d1)
        trio = generate_trio(pedigree_dict)
        for i in trio:
            generate_ped(i,d1)
            ID = d1.loc[i]['sample']
            path = 'trio/'+ID
            if not os.path.exists(path):
                os.makedirs(path)
            os.system('cp /DATA/sslyu/trio_BWA-GATK_3.0/src/{mendel_to_annovar.py,sort_sample.py,vep-mendelscan.sh,vfilt.sh} '+path)
            trio_cmd = generate_trio_cmd(path, filtmode) # trio_cmd.sh
            trio_cmd.write('[ -e ../../$sample/2_mapping/$sample.depth.sample_gene_summary ] && cp ../../$sample/2_mapping/$sample.depth.sample_gene_summary ../../annotation\n')
            trio_cmd.write('cp segtrio/$sample.{ann.hg19_multianno.txt,CADD,maf,link} ../../annotation\n')
            if not os.path.exists('annotation'):
                os.makedirs('annotation')
            trio_cmd.write('cd ../../annotation\n')
            trio_cmd.write('echo $sample >> list\n')
            trio_cmd.write('python /DATA/sslyu/trio_BWA-GATK_3.0/src/score_re.py -sn $sample\n')
            relation = []
            for j in pedigree_dict[i]:
                relation.append(d1.loc[j]['relationship'])
            a = ' '.join(relation)
            annotationfile = open('annotation/annotation.sh', 'a')
            if u'父' in a and u'母' in a:
                trio_cmd.write('python /DATA/sslyu/trio_BWA-GATK_3.0/src/annotation_filt_ver3.0.py -sn $sample -c_f_m c_f_m\n')
                annotationfile.write('python /DATA/sslyu/trio_BWA-GATK_3.0/src/annotation_filt_ver3.0.py -sn '+ID+' -c_f_m c_f_m\n')
            else:
                if u'父' in a:
                    trio_cmd.write('python /DATA/sslyu/trio_BWA-GATK_3.0/src/annotation_filt_ver3.0.py -sn $sample -c_f_m c_f\n')
                    annotationfile.write('python /DATA/sslyu/trio_BWA-GATK_3.0/src/annotation_filt_ver3.0.py -sn '+ID+' -c_f_m c_f\n')
                else:
                    trio_cmd.write('python /DATA/sslyu/trio_BWA-GATK_3.0/src/annotation_filt_ver3.0.py -sn $sample -c_f_m c_m\n')
                    annotationfile.write('python /DATA/sslyu/trio_BWA-GATK_3.0/src/annotation_filt_ver3.0.py -sn '+ID+' -c_f_m c_m\n')
            trio_cmd.write('fi\n')
            trio_cmd.close()
            annotationfile.close()
        generate_trio_list(trio, d1)

else:
    print 'please offer an information excel, the name of the excel should include \"分析任务单\", the columns should be 样本编号 原始样本ID 姓名 性别 familyname relationship'
    sys.exit()

os.system('sort trio/list |uniq > tmp; mv tmp trio/list')
if not os.path.exists('peddy'):
    os.mkdir('peddy')
if not os.path.exists('peddy/trio_peddy.sh'):    
    os.system('cp /DATA/sslyu/trio_BWA-GATK_3.0/src/trio_cmd.sh peddy')
os.system('echo peddy > peddy/list')
os.system('for i in `cat trio/list`; do cat ped/$i.mendel.ped >> peddy_tmp; done')
if not os.path.exists('ped'):
    os.mkdir('ped')
os.system('sort peddy_tmp |uniq |awk \'{OFS="\t"; print \'0\',$2,$3,$4,$5,$6}\' > ped/peddy.ped')
with open('ped/peddy.ped') as f:
	for l in f:
		print l.strip('\n')
os.system('sort annotation/annotation.sh |uniq > annotation/tmp; mv annotation/tmp annotation/annotation.sh')        
