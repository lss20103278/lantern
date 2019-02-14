#!/usr/bin/env python
# -*- coding: utf-8 -*-
#########################################################################
# File Name: fqstat_info.py
# Author: sslyu
# mail: 1573077420@qq.com
# Created Time:Thu 27 Sep 2018 09:52:15 AM CST
#########################################################################

import sys
reload(sys)
sys.setdefaultencoding("utf8")
import os
import pandas as pd

kwargs_raw = sys.argv[1:]
if len(kwargs_raw) == 0:
    print """
    Examples:
    python /DATA/sslyu/trio_BWA-GATK_3.0/src/fqstat_info.py -sn sampleID
    python /DATA/sslyu/trio_BWA-GATK_3.0/src/fqstat_info.py -l list
    """
    sys.exit()

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

ReadNum = []
BaseNum = []
Q30 = []

def extract_info(file):
    with open(file) as f:
        Q30_base = 0.0
        for l in f:
            if 'ReadNum' in l:
                readnum = l.split(': ')[1].split('\t')[0]
                if 'BaseNum' in l:
                    basenum = float(l.split(': ')[2].split('\t')[0])*2/1000000.0
            if 'BaseQ' in l and 'Q30' in l:
                Q30_base = Q30_base + float(l.strip('\n').split('Q30: ')[1][:-1])
    info = (readnum,'{:.2f}'.format(basenum),'{:.2f}'.format(Q30_base/2))
    return info

for sn in list:
    info = extract_info(sn+'.info')
    print sn+'\t'+'\t'.join(info)

######################################################################################################
#if(kwargs.get('-sn')!=''):
#    infile = kwargs.get('-sn')
#
#if 'tmp' not in infile:
#    d1 = append_sample_gender_excel(infile)
#else:
#    d1 = pd.read_table(infile, dtype=str, encoding='utf-8')
#
#d1.index = d1['sample']
#sn = d1['sample'].tolist()
#loss = []
#fq_list = d1.columns.tolist()[16:19]
#if 'tmp' in infile:
#    for i in range(len(d1.index)):
#        for j in fq_list:
#            if d1[j].isnull()[i]:
#                loss.append(d1['sample'][i])
#                break
#else:
#    for i in d1[u'sample']:
#        for j in fq_list:
#            if d1.loc[i][j] == 'nan':
#                loss.append(i)
#                break
#
#for i in loss:
#    with open('raw/'+i+'.info') as f:
#        Q30_base = 0.0
#        for l in f:
#            if 'ReadNum' in l:
#                readnum = l.split(': ')[1].split('\t')[0]
#                if 'BaseNum' in l:
#                    basenum = float(l.split(': ')[2].split('\t')[0])*2/1000000000.0
#            if 'BaseQ' in l and 'Q30' in l:
#                Q30_base = Q30_base + float(l.strip('\n').split('Q30: ')[1][:-1])
#        #ReadNum.append('{:.2f}'.format(readnum))
#        ReadNum.append(readnum)
#        BaseNum.append('{:.2f}'.format(basenum))
#        Q30.append('{:.2f}'.format(Q30_base/2)+'%')
#
#df = pd.DataFrame(columns=['Reads', 'bases(Gb)', 'Q30'])
#df['Reads'] = ReadNum
#df['bases(Gb)'] = BaseNum
#df['Q30'] = Q30
#df.index = loss
#
#ex_list = [i for i in d1.index.tolist() if i not in loss]
#if len(ex_list) > 0:
#    d1_tmp = d1.loc[ex_list][d1.columns.tolist()[16:19]]
#    d_all = pd.concat([d1_tmp,df], axis=0)
#else:
#    d_all = df
#d_all.to_csv('all_fqstat', sep="\t")

