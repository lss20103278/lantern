#!/usr/bin/env python
# -*- coding: utf-8 -*-
#########################################################################
# File Name: /DATA/sslyu/trio_BWA-GATK_2.7/check.py
# Author: sslyu
# mail: 1573077420@qq.com
# Created Time:Wed 08 Aug 2018 01:46:07 PM CST
#########################################################################

import os
import filecmp
import sys
reload(sys)
sys.setdefaultencoding("utf8")
import pandas as pd
import numpy as np

sys.path.append("/DATA/sslyu/trio_BWA-GATK_3.0/")
from lib.prepare import *

def find_file(prefix, name):
    target_file = []
    for relpath, dirs, files in os.walk(prefix):
        if name in files:
            full_path = os.path.join(prefix, relpath, name)
            full_path = os.path.normpath(os.path.abspath(full_path))
            target_file.append(full_path)
    return target_file

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

if len(list) == 0:
    print """
    Examples: python /DATA/sslyu/trio_BWA-GATK_3.0/src/check.py -sn example.xlsx
              python /DATA/sslyu/trio_BWA-GATK_3.0/src/check.py -sn list_excel
    """
    sys.exit()
else:
    for sn in list:
        print sn
        excel = pd.read_excel(sn)
        excel = append_sample_excel(excel)
        d1 = append_gender(excel)
        print d1
        if u'样本编号' in d1.columns:
            colname = u'样本编号'
        if u'样本ID' in d1.columns:
            colname = u'样本ID'
        d1.index = range(len(d1['sample']))

        strategy = os.getcwd().split('/')[-1].split('_')[-1]
        print 'strategy: '+strategy
        truth = raw_input('Is the strategy right? If not, please enter it (WES IDT-PANEL PANEL Agilent-wes ... ): ')
        if truth != 'yes':
            strategy = truth

        md5_check = raw_input('Do you need to check the md5 of the copied files? If not, please enter no:')
        if md5_check != 'no':
            check_md5_file = open('check_md5.sh', 'w')
            original_rawdata_dir_single_judgement = raw_input('Is the original rawdata in a single dir or in a set of subdirs? If single, please enter single, if not, please enter set:  ')
            if original_rawdata_dir_single_judgement == 'single':
                original_rawdata_dir = raw_input('Please enter the absolute path of the dir of the original rawdata:  ')
                os.system('ls '+original_rawdata_dir)
                pattern = raw_input('please enter the pattern of R1 and R2:  ')
                md5 = raw_input('please enter the name of the md5 file:  ')
                check_md5_file.write('for i in `cut -f 1 rename0`; do for j in '+pattern+'; do a=`md5sum '+original_rawdata_dir+'/*$i*$j |cut -d" " -f 1`; b=`grep $i '+original_rawdata_dir+'/'+md5+' |grep $j |cut -d" " -f 1`; if [ "$a" != "$b" ]; then echo -e $i"_"$j" is problematic"; fi; done; done\n')
                if not os.path.exists('cp_md5'+str(list.index(sn))):
                    check_md5_file.write('for i in `cat list`; do for j in R1 R2; do md5sum raw/$i\_$j.fastq.gz >> cp_md5; done; done\n')
                check_md5_file.write('original_rawdata_dir='+original_rawdata_dir+'; md5='+md5+'; pattern_original=('+pattern+'); pattern_copied=(R1 R2); for i in `cat rename0 |tr "\\t" "-"`; do ID=${i%-*}; sample=${i#*-}; for j in 0 1; do original=`grep $ID $original_rawdata_dir/$md5 |grep ${pattern_original[$j]} |cut -d" " -f 1`; copied=`grep $sample cp_md5 |grep ${pattern_copied[$j]} |cut -d" " -f 1`; if [ "$original" != "$copied" ]; then echo -e $ID"_"${pattern_original[$j]}" is not copied correctly"; fi; done; done\n')
                check_md5_file.write('original_rawdata_dir='+original_rawdata_dir+'; md5='+md5+'; pattern=('+pattern+'); for i in `cut -f 1 rename0`; do for j in 0 1; do original=`grep $i $original_rawdata_dir/$md5 |grep ${pattern[$j]} |cut -d" " -f 1`; backup=`md5sum /anaData/anaData004/children_hos_genetic/rawdata/2018/'+strategy+'/$i/*$i*${pattern[$j]} |cut -d" " -f 1`; if [ "$original" != "$backup" ]; then echo -e $i"_"${pattern[$j]}"is not backup correctly"; fi; done; done\n')
                check_md5_file.close()
            else:
                sys.exit()

        slurm_check = raw_input('Check slurm?')
        if slurm_check != 'no':
            ofile = open('check_slurm.sh', 'w')
            ofile.write('for i in `find . -name \"slurm*\"`; do a=`grep -ci error $i`; if [ $a -gt 0 ]; then echo -e $i has errors; grep -i error -n $i; fi; done\n')
            ofile.close()
            os.system('bash check_slurm.sh')
            #os.system('rm check.sh')
        
        qcsum_check = raw_input('Check qcsum?')
        if qcsum_check != 'no':
            if os.path.exists('all.qcsum'):
                d2 = pd.read_table('all.qcsum', dtype=str)
                d2.index = d2['sample']
            	#os.system('cat all.qcsum')
            else:
                loss = []
                for i in d1['sample']:
                    if not os.path.exists(i+'/2_mapping/'+i+'.qcsum'):
                        loss.append('qcsum of '+i+' does not exist')
                if len(loss) > 0:
                    for i in loss:
                        print i
                    sys.exit()
                d2 = pd.DataFrame(columns=['sample', '1Xcov', '20Xcov', 'AvgDepth', 'duplication'])
                for j in d1['sample']:
                    fn = j+'/2_mapping/'+j+'.qcsum'
                    qcsum = pd.read_table(fn)
                    d2 = pd.concat([d2,qcsum], ignore_index=True)
                d2.index = d2['sample']
                d2.to_csv('all.qcsum', sep="\t", index=None)
            with open('dbevn.sh') as f:
                for l in f:
                    if l.startswith('panel'):
                        l = l.strip('\n').split(' ')
                        strategy = l[2]
            strategy = os.getcwd().split('/')[-1].split('_')[-1]
            print 'strategy: '+strategy
            truth = raw_input('Is the strategy right? If not, please enter it (WES PANEL IDT-PANEL Agilent-PANEL Agilent-wes ... ): ')
            if truth != '':
                strategy = truth
            if strategy == 'PANEL' or strategy == 'Agilent-PANEL' or strategy == 'IDT-PANEL':
                for i in d2.index:
                    if float(d2.loc[i]['20Xcov']) < 95:
                        print 'the 20X coverage of '+i+' is less than 95%'
            else:
                for i in d2.index:
                    if float(d2.loc[i]['20Xcov']) < 90:
                        print 'the 20X coverage of '+i+' is less than 90%'
        
        sex_ped_check = raw_input('Check sex and pedigree?')
        if sex_ped_check != 'no':
            if os.path.exists('peddy/pedtrio/results/peddy.sex_check.csv'):
                sex = pd.read_csv('peddy/pedtrio/results/peddy.sex_check.csv', dtype=str)
                for i in sex.index:
                    if sex.loc[i]['error'] == 'True':
                        print 'the sex of '+sex.loc[i]['sample_id']+' is not the same with the predicted result'
            if os.path.exists('peddy/pedtrio/results/peddy.ped_check.csv'):
                ped = pd.read_csv('peddy/pedtrio/results/peddy.ped_check.csv', dtype=str)
                for i in ped.index:
                    if ped.loc[i]['parent_error'] == 'True':
                        print ped.loc[i]['sample_a']+' and '+ped.loc[i]['sample_b']+' are not child-parent'
            if os.path.exists('peddy/pedtrio/results/peddy.html'):
                os.system('mv peddy/pedtrio/results/peddy.html .')
        
        report_check = raw_input('Check the final reporter dir?')
        if report_check != 'no':
            rootDir = os.getcwd().split('/')[-1]
            print rootDir
            rootDir_truth = raw_input('Is the dir name correct? If not, please enter the name:')
            if rootDir_truth != '':
                rootDir = rootDir_truth
            if os.path.exists(rootDir):
                all_check_file = flatten_list(rootDir)
                for i in all_check_file:
                    print unicode(i,'utf-8'),
                wrong = []
                for i in all_check_file:
                    target_file = find_file('.', i)
                    if len(target_file) > 2:
                        print target_file
                        print 'more than 2 files are found'
                        sys.exit()
                    elif len(target_file) < 2:
                        print target_file
                        print 'only one file are found'
                    else:
                        if not cmp(target_file[0],target_file[1]):
                            wrong.extend(target_file)
                if len(wrong) > 0:
                    for i in wrong:
                        print wrong[i][0]+' is not the same with '+wrong[i][1]
                    sys.exit()

#        md5_check = raw_input('Do you need to check the md5 of the copied files? If not, please enter no:')
#        if md5_check != 'no':
#            cp_md5 = {}
#            if os.path.exists('raw/cp_md5'+str(list.index(sn))):
#                md5_cp = open('raw/cp_md5'+str(list.index(sn))).readlines()
#                for i in d1['sample']:
#                    for j in ['_R1', '_R2']:
#                        for k in md5_cp:
#                            if i in k and j in k:
#                                cp_md5[i+j] = k.split()[0]
#            else:
#                for i in d1.index:
#                    for j in ['_R1', '_R2']:
#                        if os.path.exists('raw/'+d1.loc[i]['sample']+j+'.fastq.gz'):
#                            v = os.popen('md5sum raw/'+d1.loc[i]['sample']+j+'.fastq.gz').readlines()[0].split()[0]
#                            cp_md5[d1.loc[i]['sample']+j] = v
#                print cp_md5
#                ofile = open('raw/cp_md5'+str(list.index(sn)), 'w')
#                for k in cp_md5:
#                    ofile.write(cp_md5[k]+' '+k+'\n')
#                ofile.close()
#            redo_judgement = raw_input('Has the raw data been copied to the backup space? If not, please enter no:')
#            if redo_judgement == 'no':
#                single_judgement = raw_input('Is the raw data in a single dir or in a set of subdirs? If single, please enter single, if not, please enter set:  ')
#                if single_judgement == 'single':
#                    rawdata_dir = raw_input('Please enter the absolute path of the raw data:  ')
#                    os.system('ls '+rawdata_dir)
#                    #print os.listdir(rawdata_dir)
#                    pattern = raw_input('please enter the pattern of R1 and R2:  ')
#                    origin_md5 = {}
#                    md5 = raw_input('please enter the absolute path of the md5 file:  ')
#                    os.system('cp '+md5+' raw/origin_md5'+str(list.index(sn)))
#                    md5_origin = open(md5).readlines()
#                    for i in d1.index:
#                        for j in pattern.split():
#                            if '1' in j:
#                                if 'R1' in j:
#                                    for k in md5_origin:
#                                        if d1.loc[i][colname] in k and j in k:
#                                            origin_md5[d1.loc[i][colname]+'_R1'] = k.split()[0]
#                                if 'R2' in j:
#                                    for k in md5_origin:
#                                        if d1.loc[i][colname] in k and j in k:
#                                            origin_md5[d1.loc[i][colname]+'_R2'] = k.split()[0]
#                                if '1.fq' in j:
#                                    for k in md5_origin:
#                                        if d1.loc[i][colname] in k and j in k:
#                                            origin_md5[d1.loc[i][colname]+'_R1'] = k.split()[0]
#                            else:
#                                for k in md5_origin:
#                                    if d1.loc[i][colname] in k and j in k:
#                                        origin_md5[d1.loc[i][colname]+'_R2'] = k.split()[0]
#                else:
#                    parent_raw_dir = raw_input('please enter the parent directory:  ')
#                    for i in d1[colname]:
#                        rawdata_dir = parent_raw_dir+'/*'+i+'*'
#                        os.system('ls '+rawdata_dir)
#                    origin = raw_input('please enter the pattern of R1 and R2:  ')
#                    md5 = raw_input('please enter the pattern of md5:  ')
#                    origin_md5 = {}
#                    for i in d1[colname]:
#                        if os.path.exists(rawdata_dir+'/*'+i+'*'):
#                            md5_origin = open(rawdata_dir+'/*'+i+'*/*'+md5+'*').readlines()
#                            for j in origin.split():
#                                if '1' in j:
#                                    for k in md5_origin:
#                                        if i in k and j in k:
#                                            origin_md5[i+'_R1'] = k.split()[0]
#                                else:
#                                    for k in md5_origin:
#                                        if i in k and j in k:
#                                            origin_md5[i+'_R2'] = k.split()[0]
#            else:
#                rawdata_dir = raw_input('Please enter the absolute path of the raw data:  ')
#                for i in d1[colname]:
#                    if os.path.exists(rawdata_dir+'/'+i):
#                        print os.listdir(rawdata_dir+i)
#                        os.system('cp '+rawdata_dir+'/'+i+'/'+i+'.md5 raw/')
#                pattern = raw_input('please enter the pattern of R1 and R2:  ')
#                origin_md5 = {}
#                for i in d1.index:
#                    if os.path.exists(rawdata_dir+'/'+i):
#                        md5_origin = open(rawdata_dir+'/'+d1.loc[i][colname]+'/'+d1.loc[i][colname]+'.md5').readlines()
#                        for j in pattern.split():
#                            if '1' in j:
#                                for k in md5_origin:
#                                    if d1.loc[i][colname] in k and j in k:
#                                        origin_md5[d1.loc[i][colname]+'_R1'] = k.split()[0]
#                            else:
#                                for k in md5_origin:
#                                    if d1.loc[i][colname] in k and j in k:
#                                        origin_md5[d1.loc[i][colname]+'_R2'] = k.split()[0]
#            wrong = []
#            for i in d1.index:
#                for j in ['_R1', '_R2']:
#                    if cp_md5.has_key(d1.loc[i]['sample']+j) and origin_md5.has_key(d1.loc[i][colname]+j):
#                        if origin_md5[d1.loc[i][colname]+j] != cp_md5[d1.loc[i]['sample']+j]:
#                            wrong.append(d1.loc[i][colname]+j+' '+d1.loc[i]['sample']+j)
#            if len(wrong) > 0:
#                for i in wrong:
#                    print i.split()[1]+' is not the same with '+i.split()[0]
#                sys.exit()
#            else:
#                print 'raw data is copied correctly'

