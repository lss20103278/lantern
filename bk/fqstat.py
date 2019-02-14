#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
reload(sys)
sys.setdefaultencoding("utf8")
import os
import pandas as pd
from prepare import *

def list_file(rootDir):
    all_file = []
    for lists in os.listdir(rootDir):
        path = os.path.join(rootDir, lists)
        if os.path.isdir(path):
            list_file(path)
        else:
            all_file.append(lists)
    return all_file

kwargs_raw = sys.argv[1:]
kwargs={'-sn':'', '-l':''}
for i in range(len(kwargs_raw)):
    for j in kwargs.keys():
        if(kwargs_raw[i]==j):
            kwargs.update({j:kwargs_raw[i+1]})

ReadNum = []
BaseNum = []
Q30 = []
#sn = pd.read_table('list', header=None)
if(kwargs.get('-sn')!=''):
    infile = kwargs.get('-sn')
if 'tmp' not in infile:    
    d1 = append_sample_gender_excel(infile)
else:
    d1 = pd.read_table(infile, dtype=str, encoding='utf-8')
d1.index = d1['sample']
sn = d1['sample'].tolist()
loss = []
if 'tmp' in infile:
    for i in range(len(d1.index)):
        #for j in [u'Reads(M)', 'bases(Gb)', 'Q30']:
        sn_list = d1.columns.tolist()[16:19]
        for j in sn_list:
            if d1[j].isnull()[i]:
                loss.append(d1['sample'][i])
                break
else:                
    for i in d1[u'sample']:
        #for j in [u'Reads(M)', 'bases(Gb)', 'Q30']:
        sn_list = d1.columns.tolist()[16:19]
        for j in sn_list:
            if d1.loc[i][j] == 'nan':
                loss.append(i)
                break
if len(loss) > 0:            
    ofile = open('fqstat.sh', 'w')
    if int(os.popen('ls raw |grep -ci info').readlines()[0].strip()) == len(sn):
        ofile.close()
        loss = sn
        for i in sn:
            with open('raw/'+i+'.info') as f:
                Q30_base = 0.0
                for l in f:
                    if 'ReadNum' in l:
                        #readnum = float(l.split(': ')[1].split('\t')[0])/1000000.0
                        readnum = l.split(': ')[1].split('\t')[0]
                        if 'BaseNum' in l:
                            basenum = float(l.split(': ')[2].split('\t')[0])*2/1000000000.0
                    if 'BaseQ' in l and 'Q30' in l:
                        Q30_base = Q30_base + float(l.strip('\n').split('Q30: ')[1][:-1])
                ReadNum.append('{:.2f}'.format(readnum))
                BaseNum.append('{:.2f}'.format(basenum))
                Q30.append('{:.2f}'.format(Q30_base/2)+'%')
    elif int(os.popen('ls raw |grep -ci info').readlines()[0].strip()) > 0:
        if len(os.popen('ls raw/*gz').readlines())/2 != len(sn):
            rawdata_dir = raw_input('please enter the absolute path of raw data:')
            d1.index = d1[u'样本编号']
            for i in d1[u'样本编号']:
                print os.listdir(rawdata_dir+'/'+i)
            pattern = raw_input('please enter the pattern of R1 and R2:')
            for i in loss:
                ofile.write('/home/ana005/anaconda2/bin/iTools Fqtools stat -InFq '+rawdata_dir+'/'+i+'/*'+pattern.split()[0]+'* -InFq '+rawdata_dir+'/'+i+'/*'+pattern.split()[1]+'* -OutStat raw/'+d1.loc[i]['sample']+'.info\n')
            ofile.close()
        else:    
            for i in loss:
            	ofile.write('/home/ana005/anaconda2/bin/iTools Fqtools stat -InFq raw/'+i+'*R1*gz -InFq raw/'+i+'*R2*gz -OutStat raw/'+i+'.info\n')
            ofile.close()
        print loss
    else:
        loss = sn
        ofile = open('fqstat.sh', 'w')
        if len(os.popen('ls raw/*gz').readlines())/2 != len(sn):
            rawdata_dir = raw_input('please enter the absolute path of raw data:')
            d1.index = d1[u'样本编号']
            for i in d1[u'样本编号']:
                print os.listdir(rawdata_dir+'/'+i)
            pattern = raw_input('please enter the pattern of R1 and R2:')
            for i in loss:
                ofile.write('/home/ana005/anaconda2/bin/iTools Fqtools stat -InFq '+rawdata_dir+'/'+i+'/*'+pattern.split()[0]+'* -InFq '+rawdata_dir+'/'+i+'/*'+pattern.split()[1]+'* -OutStat raw/'+d1.loc[i]['sample']+'.info\n')
            ofile.close()
        else:    
            for i in loss:
            	ofile.write('/home/ana005/anaconda2/bin/iTools Fqtools stat -InFq raw/'+i+'*R1*gz -InFq raw/'+i+'*R2*gz -OutStat raw/'+i+'.info\n')
            ofile.close()
        print loss
if os.path.getsize('fqstat.sh') > 0:        
    with open('fqstat.sh') as f:
        for l in f:
            print l
    judgement = raw_input('Is the shell script correct? If not, please enter no: ')
    if judgement != 'no':
        os.system('bash fqstat.sh')
        for i in loss:
            with open('raw/'+i+'.info') as f:
                Q30_base = 0.0
                for l in f:
                    if 'ReadNum' in l:   
                        #readnum = float(l.split(': ')[1].split('\t')[0])/1000000.0
                        readnum = l.split(': ')[1].split('\t')[0]
                        if 'BaseNum' in l:
                            basenum = float(l.split(': ')[2].split('\t')[0])*2/1000000000.0
                    if 'BaseQ' in l and 'Q30' in l:
                        Q30_base = Q30_base + float(l.strip('\n').split('Q30: ')[1][:-1])
                #ReadNum.append('{:.2f}'.format(readnum))
                ReadNum.append(readnum)
                BaseNum.append('{:.2f}'.format(basenum))	
                Q30.append('{:.2f}'.format(Q30_base/2)+'%')

#df = pd.DataFrame(columns=['Reads(M)', 'bases(Gb)', 'Q30'])
#df['Reads(M)'] = ReadNum
df = pd.DataFrame(columns=['Reads', 'bases(Gb)', 'Q30'])
df['Reads'] = ReadNum
df['bases(Gb)'] = BaseNum
df['Q30'] = Q30
df.index = loss
#loss_ID = [d1.loc[i]['sample'] for i in loss] #
#ex_list = [d1.loc[i]['sample'] for i in d1.index.tolist() if i not in loss]
#d1_tmp = d1.loc[ex_list]
#df.index = loss_ID    #
#d_all = pd.concat([d1_tmp,df], axis=1)
#d_all.to_csv('all_fqstat', sep="\t")             #
ex_list = [i for i in d1.index.tolist() if i not in loss]
if len(ex_list) > 0:
    d1_tmp = d1.loc[ex_list][d1.columns.tolist()[16:19]]
    d_all = pd.concat([d1_tmp,df], axis=0)
else:
    d_all = df
d_all.to_csv('all_fqstat', sep="\t")
