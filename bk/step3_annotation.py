#!/usr/bin/env python
# -*- coding: utf-8 -*-
#########################################################################
# File Name: /DATA/sslyu/trio_BWA-GATK_2.7/annotation.py
# Author: sslyu
# mail: 1573077420@qq.com
# Created Time:Tue 07 Aug 2018 03:12:48 PM CST
#########################################################################

import numpy as np
import pandas as pd
import sys
reload(sys)
sys.setdefaultencoding("utf8")
import os
import re
from prepare import *

#if len(sys.argv) < 2:
#	print 'miss the value of thisminedge...'
#	sys.exit()

kwargs_raw = sys.argv[1:]
kwargs={'-sn':'', '-l':''}
for i in range(len(kwargs_raw)):
    for j in kwargs.keys():
        if(kwargs_raw[i]==j):
            kwargs.update({j:kwargs_raw[i+1]})

#filelist = os.listdir('.')
#list = []
#for i in filelist:
#    if 'xlsx.tmp' in i:
#        list.append(i)
#    elif '分析任务单' in i:
#        list.append(i)
#    elif '产品信息表' in i:
#        list.append(i)
##for i in filelist:
##    if '分析任务单' in i and 'tmp' not in i:
##        list.append(i)
##    if '产品信息表' in i and 'tmp' not in i:
##        list.append(i)

#read sample list
if(kwargs.get('-sn')!=''):
    list=[]
    list.append(kwargs.get('-sn'))
    print "\nNote:sigle sample mode"
else:
    if kwargs.get('-l')!='':
        list=pd.read_table(kwargs.get('-l'),header=None)[0].tolist()
        print "\nNote:list samples mode"

if len(list) == 0:
	print 'please offer an information excel, the name of the excel should include\"分析任务单\", the columns should be 样本编号 原始样本ID 姓名 性别 pedigree relationship'
	sys.exit()
	
for sn in list:	
    if 'tmp' in sn:
        d1 = pd.read_table(sn, dtype=str, encoding='utf-8')
    else:
        d1 = append_sample_gender_excel(sn)
        d1 = append_relation(d1)
        d1 = append_pedigree(d1)
        d1 = append_phenotype(d1)
        d1 = append_father_mother(d1)
        infile = sn.split(')')[0]+u'）'+sn.split(')')[1]
        d1.to_csv(infile+'.tmp', index=None, sep="\t", encoding='utf-8') # the encoding of tmp file is problematic
    generate_annotation(d1)

with open('step3_annotation.sh') as f:
	for l in f:
		print l.strip('\n')
#truth = raw_input('Is the shell script correct? If not, please enter no: ')
#if truth == 'no':
#	sys.exit()
#else:
#	os.system('sh step3_annotation.sh')
