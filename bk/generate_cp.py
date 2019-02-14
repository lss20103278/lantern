#!/usr/bin/env python
# -*- coding: utf-8 -*-
#########################################################################
# File Name: generate_cp.py
# Author: sslyu
# mail: 1573077420@qq.com
# Created Time:Fri 28 Sep 2018 05:36:40 PM CST
#########################################################################

import numpy as np
import pandas as pd
import sys
reload(sys)
sys.setdefaultencoding("utf8")
import os
import re
import collections

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
    else:
        print """
        examples:
        python generate_cp.py -sn excelname
        python generate_cp.py -l \'listname of list of excelname\'
        """
        sys.exit()

def append_sample_gender_excel(excel):
    df = pd.read_excel(excel, dtype=str)
    ID = []
    for i in range(len(df[u'原始样本ID'])):
        if df[u'原始样本ID'][i] == 'nan':
            ID.append(df.iloc[i][u'样本编号'])
        else:
            ID.append(df.iloc[i][u'原始样本ID'])
    df['sample'] = ID
    df['gender'] = df[u'性别'].apply(lambda x:'1' if u'男' in x else '2')
    return df

if not os.path.exists('raw'):
    os.system('mkdir raw')    
os.system('cp /DATA/sslyu/trio_BWA-GATK_3.0/sub_cp_{list,single}.sh .')
for sn in list:
    df = append_sample_gender_excel(sn)
    df.index = df['sample']
    redo_judement = raw_input('Is this a redo work? If not, please enter no:  ')
    if redo_judement == 'no':
        single_judgement = raw_input('Is the raw data in a single dir or in a set of subdirs? If single, please enter single, if set, please enter set:')
        if single_judgement == 'single':
            rawdata_dir = raw_input('Please enter the absolute path of the raw data:  ')
            if not rawdata_dir.startswith('/DATA'):
                rawdata_dir = os.getcwd()+'/'+rawdata_dir
            os.system('ls '+rawdata_dir)
            origin = raw_input('please enter the pattern of R1 and R2:  ')
            md5 = raw_input('please enter the absolute path of the md5 file:  ')
            if not md5.startswith('/DATA'):
                md5 = os.getcwd()+'/'+md5
            for i in df['sample']:
                if not os.path.exists(str(i)):
                    os.mkdir(str(i))
                cp_rawdata_file = open(str(i)+'/cp_rawdata.sh','w')
                cp_rawdata_file.write('#!/bin/sh\n\n#SBATCH -J '+i+'\n#SBATCH -p BP10\n#SBATCH -N 1\n#SBATCH -n 4\n\nsample='+i+'\n\nmkdir -p ../raw/bk/'+df.loc[i][u'样本编号']+'\nawk \'{if(match($2,\"\'$sample\'\")){print $0}}\' '+md5+' > ../raw/bk/'+df.loc[i][u'样本编号']+'/'+df.loc[i][u'样本编号']+'.md5\ncp '+rawdata_dir+'/*'+df.loc[i][u'样本编号']+'*'+origin.split(' ')[0]+'* ../raw/bk/'+df.loc[i][u'样本编号']+'\ncp '+rawdata_dir+'/*'+df.loc[i][u'样本编号']+'*'+origin.split(' ')[1]+'* ../raw/bk/'+df.loc[i][u'样本编号']+'\ncp '+rawdata_dir+'/*'+df[u'样本编号'][i]+'*'+origin.split(' ')[0]+'* ../raw/'+i+'_R1.fastq.gz\ncp '+rawdata_dir+'/*'+df[u'样本编号'][i]+'*'+origin.split(' ')[1]+'* ../raw/'+i+'_R2.fastq.gz\n')
                cp_rawdata_file.close()
        else:
            parent_raw_dir = raw_input('please enter the parent directory:  ')
            for i in df[u'样本编号']:
                rawdata_dir = parent_raw_dir+'/*'+i+'*'
                os.system('ls '+rawdata_dir)
            origin = raw_input('please enter the pattern of R1 and R2:  ')
            md5 = raw_input('please enter the pattern of md5:  ')
            for i in df['sample']:
                if not os.path.exists(str(i)):
                    os.mkdir(str(i))
                cp_rawdata_file = open(str(i)+'/cp_rawdata.sh','w')
                cp_rawdata_file.write('#!/bin/sh\n\n#SBATCH -J '+i+'\n#SBATCH -p BP10\n#SBATCH -N 1\n#SBATCH -n 4\n\nsample='+i+'\n\nmkdir -p ../raw/bk/'+df.loc[i][u'样本编号']+'\ncp '+parent_raw_dir+'/*'+df.loc[i][u'样本编号']+'*/*'+md5+'* ../raw/bk/'+df.loc[i][u'样本编号']+'/'+df.loc[i][u'样本编号']+'.md5\ncp '+parent_raw_dir+'/*'+df.loc[i][u'样本编号']+'*/*gz ../raw/bk/'+df.loc[i][u'样本编号']+'/\ncp '+parent_raw_dir+'/*'+df.loc[i][u'样本编号']+'*/*'+df.loc[i][u'样本编号']+'*'+origin.split(' ')[0]+'* ../raw/'+i+'_R1.fastq.gz\ncp '+parent_raw_dir+'/*'+df.loc[i][u'样本编号']+'*/*'+df.loc[i][u'样本编号']+'*'+origin.split(' ')[1]+'* ../raw/'+i+'_R2.fastq.gz\n')
                cp_rawdata_file.close()
    else:
        parent_raw_dir = raw_input('please enter the parent directory:  ')
        for i in df[u'样本编号']:
            rawdata_dir = parent_raw_dir+'/'+i
            os.system('ls '+rawdata_dir)
        origin = raw_input('please enter the pattern of R1 and R2:  ')
        for i in df['sample']:
            if not os.path.exists(str(i)):
                os.mkdir(str(i))
            cp_rawdata_file = open(str(i)+'/cp_rawdata.sh','w')
            cp_rawdata_file.write('#!/bin/sh\n\n#SBATCH -J '+i+'\n#SBATCH -p BP10\n#SBATCH -N 1\n#SBATCH -n 4\n\nsample='+i+'\n\ncp '+parent_raw_dir+'/'+df.loc[i][u'样本编号']+'/*'+df.loc[i][u'样本编号']+'*'+origin.split(' ')[0]+'* ../raw/'+i+'_R1.fastq.gz\ncp '+parent_raw_dir+'/'+df.loc[i][u'样本编号']+'/*'+df.loc[i][u'样本编号']+'*'+origin.split(' ')[1]+'* ../raw/'+i+'_R2.fastq.gz\n')
            cp_rawdata_file.close()
                    

