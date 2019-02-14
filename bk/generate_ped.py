#!/usr/bin/env python
# -*- coding: utf-8 -*-
#########################################################################
# File Name: generate_ped.py
# Author: sslyu
# mail: 1573077420@qq.com
# Created Time:Fri 28 Sep 2018 05:15:56 PM CST
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
kwargs={'-sn':'', '-strategy':'', '-filtmode':''}
for i in range(len(kwargs_raw)):
    for j in kwargs.keys():
        if(kwargs_raw[i]==j):
            kwargs.update({j:kwargs_raw[i+1]})

if len(kwargs_raw) != 6:
    print """
    examples:
    python generate_ped.py -sn filename -strategy (WES/IDT-PANEL...) -filtmode (hard/vqsr)
    """
    sys.exit()

#read sample list
if(kwargs.get('-sn')!=''):
    sn = kwargs.get('-sn')
    print "\nNote:sigle sample mode"

strategy = kwargs.get('-strategy')
filtmode = kwargs.get('-filtmode')

def generate_ped(k,df):
    df.index = df['sample']
    pedigree_k = df['sample'][df['familyname'] == k].tolist()
    ped_f = df.loc[pedigree_k][['familyname', 'sample', 'father', 'mother', 'gender', 'phenotype1']]
    ped_f = ped_f.loc[ped_f['familyname'] == k] # incase one parent have two sick children
    ped_f.index = ped_f['sample']
    ped_mendel = df.loc[pedigree_k][['familyname', 'sample', 'father', 'mother', 'gender', 'phenotype2']]
    ped_mendel = ped_mendel.loc[ped_mendel['familyname'] == k] 
    ped_mendel.index = ped_mendel['sample']
    if not os.path.exists('trio/ped'):
        os.makedirs('trio/ped')
    f1name = 'trio/ped/'+k+'.ped'
    f2name = 'trio/ped/'+k+'.mendel.ped'
    ped1 = open(f1name, 'w')
    ped2 = open(f2name, 'w')
    relation = ''
    for i in pedigree_k:
        relation = relation+df.loc[i]['relationship']
    c_f_m = []
    c_f_m.append(k)
    ped1.write('\t'.join(ped_f.loc[k])+'\n')
    ped2.write('\t'.join(ped_mendel.loc[k])+'\n')
    if u'父' in relation and u'母' in relation:
        for i in pedigree_k:
            if df.loc[i]['relationship'] == u'父':
                c_f_m.append(i)
                ped1.write('\t'.join(ped_f.loc[i])+'\n')
                ped2.write('\t'.join(ped_mendel.loc[i])+'\n')
                break
        for i in pedigree_k:
            if df.loc[i]['relationship'] == u'母':
                c_f_m.append(i)
                ped1.write('\t'.join(ped_f.loc[i])+'\n')
                ped2.write('\t'.join(ped_mendel.loc[i])+'\n')
                break
    elif u'父' in relation or u'母' in relation:
        parent = ''
        if u'父' in relation:
            parent = u'父'
        else:
            parent = u'母'
        for i in pedigree_k:
            if df.loc[i]['relationship'] == parent:
                c_f_m.append(i)
                ped1.write('\t'.join(ped_f.loc[i])+'\n')
                ped2.write('\t'.join(ped_mendel.loc[i])+'\n')
    n_c_f_m = [i for i in pedigree_k if i not in c_f_m]
    if len(n_c_f_m) > 0:
        for i in n_c_f_m:
            ped1.write('\t'.join(ped_f.loc[i])+'\n')
            ped2.write('\t'.join(ped_mendel.loc[i])+'\n')
    ped1.close()
    ped2.close()

def generate_single_cmd(path):
    single_cmd = open(path+'/cmd.sh', 'w')
    with open('/DATA/sslyu/trio_BWA-GATK_3.0/cmd.sh') as f:
        for l in f:
            single_cmd.write(l)
    single_cmd.write('cp 2_mapping/$sample.depth.sample_gene_summary ../annotation\n')
    single_cmd.write('cp 3_variants/$sample.{ann.hg19_multianno.txt,CADD,maf,link} ../annotation\n')
    single_cmd.write('cd ../annotation\n')
    single_cmd.write('echo '+path+' >> list\n')
    single_cmd.write('python /DATA/sslyu/trio_BWA-GATK_3.0/annotation/score_re.py -sn $sample\n')
    single_cmd.write('python /DATA/sslyu/trio_BWA-GATK_3.0/annotation/annotation_filt_ver3.0.py -sn $sample\n')
    single_cmd.close()

def generate_trio_cmd(path,filtmode):
    trio_cmd = open(path+'/cmd.sh', 'w')
    n = 0
    with open('/DATA/sslyu/trio_BWA-GATK_3.0/trio/cmd.sh') as f:
        for l in f:
            n = n+1
            if n<=22:
                if n == 12:
                    trio_cmd.write('filtmode=\''+filtmode+'\'\n')
                else:
                    trio_cmd.write(l)
    return trio_cmd

def generate_dbevn(strategy):
    dbevn = open('dbevn.sh', 'w')
    if strategy != 'none':
        with open('/DATA/sslyu/trio_BWA-GATK_3.0/dbevn.sh') as f:
            for l in f:
                if l.startswith('#panel'):
                    if strategy in l:
                        dbevn.write(l[1:])
                    else:
                        dbevn.write(l)
                else:
                    dbevn.write(l)
    else:
        with open('/DATA/sslyu/trio_BWA-GATK_3.0/dbevn.sh') as f:
            for l in f:
                dbevn.write(l)
        panel = raw_input('please enter the absolute path of the bed:')
        dbevn.write('panel=\"'+panel+'\"\n')
    dbevn.close()

def generate_run2_sentieon(path,str):
    run2_sentieon = open(path+'/run2_sentieon.sh', 'w')
    n = 0
    with open('/DATA/sslyu/trio_BWA-GATK_3.0/run2_sentieon.sh') as f:
        for l in f:
            n = n+1
            if n>0 and n<= 23:
                run2_sentieon.write(l)
            elif n>= 25:
                run2_sentieon.write(l)
            elif n==24:
                run2_sentieon.write('cmd="'+str+'"\n')
            else:
                continue
    run2_sentieon.close()            

os.system('mkdir -p raw/bk trio/ped trio/gvcf trio/peddy annotation CNV')
os.system('cp /DATA/sslyu/trio_BWA-GATK_3.0/sub{_sep.sh,_list.sh,_cp_single.sh,_cp_list.sh} .')
os.system('cp /DATA/sslyu/trio_BWA-GATK_3.0/trio/gvcf.py .')
os.system('cp /DATA/sslyu/trio_BWA-GATK_3.0/cnv/R* .')

if os.path.exists('dbevn.sh'):
    with open('dbevn.sh') as f:
        for l in f:
            if l.startswith('panel'):
                print l.strip('\n')
    dbevn_truth = raw_input('Is dbevn.sh correct?')
    if dbevn_truth == 'no':
        strategy = raw_input('strategy(WES PANEL Agilent-PANEL IDT-PANEL Agilent-wes brain blood-disease, if not exist, please enter none):')
        generate_dbevn(strategy)
else:
    strategy = raw_input('strategy(WES PANEL Agilent-PANEL IDT-PANEL Agilent-wes brain blood-disease, if not exist, please enter none):')
    generate_dbevn(strategy)
    with open('dbevn.sh') as f:
        for l in f:
            if l.startswith('panel'):
                print l.strip('\n')

df = pd.read_table(sn, dtype=str, encoding='utf-8')
df['sample'].to_csv('list',index=None,header=None)
os.system('cp /DATA/sslyu/trio_BWA-GATK_3.0/sub_{list,sep}.sh .')
trio = [item for item, count in collections.Counter(df['familyname']).items() if count > 1]
pedigree = {}
for i in trio:
    pedigree[i] = []
for i in trio:
    for j in df.index:
        if df.loc[j]['familyname'] == i:
            pedigree[i].append(df.loc[j]['sample'])
    pedigree[i].append
single = [item for item, count in collections.Counter(df['familyname']).items() if count == 1]
for i in df['sample']:
    if not os.path.exists(str(i)):
        os.mkdir(str(i))
truth = raw_input('Do you need to generate cmd files?')
if truth != 'no':
    n = 0
    with open('/DATA/sslyu/trio_BWA-GATK_3.0/run2_sentieon.sh') as f:
        for l in f:
            n = n+1
            if n>22 and n<=30:
                print l.strip('\n')
    if len(single) > 0:
        list_single = open('list_single', 'w')
        cmd_single = raw_input('select cmds from following for single (qc clean mapping precalling calling gvcf index gtgvcf vqsr/hardfilt left-normalize annotation depth qcsum qcvcf): ')
        for i in single:
            list_single.write(str(i)+'\n')
            generate_single_cmd(str(i))
            generate_ped(str(i),df)
            generate_run2_sentieon(str(i),cmd_single)
        list_single.close()
    if len(trio) > 0:
        cmd_trio = raw_input('select cmds from following for trio (qc clean mapping precalling gvcf index depth qcsum): ')
        list_trio = open('list_trio', 'w')
        list = open('trio/list', 'w')
        for k in pedigree:
            list_trio.write(str(k)+'\n')
            list.write(str(k)+'\n')
            generate_ped(str(k),df)
            path = 'trio/'+k
            if not os.path.exists(path):
                os.makedirs(path)
            os.system('cp /DATA/sslyu/trio_BWA-GATK_3.0/trio/{mendel_to_annovar.sh,sort_sample.sh,vep-mendelscan.sh,vfilt.sh} '+path)
            os.system('cp /DATA/sslyu/trio_BWA-GATK_3.0/trio/trio.sh '+path)
            trio_cmd = generate_trio_cmd(path, filtmode)
            trio_cmd.write('[ -e ../../$sample/2_mapping/$sample.depth.sample_gene_summary ] && cp ../../$sample/2_mapping/$sample.depth.sample_gene_summary ../../annotation\n')
            trio_cmd.write('cp segtrio/$sample.{ann.hg19_multianno.txt,CADD,maf,link} ../../annotation\n')
            if not os.path.exists('annotation'):
                os.makedirs('annotation')
            trio_cmd.write('cd ../../annotation\n')
            trio_cmd.write('echo '+k+' >> list\n')
            trio_cmd.write('python /DATA/sslyu/trio_BWA-GATK_3.0/annotation/score_re.py -sn $sample\n')
            relation = ''
            for i in pedigree[k]:
                relation = relation+df.loc[i]['relationship']
            if u'父' in relation and u'母' in relation:
                trio_cmd.write('python /DATA/sslyu/trio_BWA-GATK_3.0/annotation/annotation_filt_ver2.8.py -sn $sample -c_f_m c_f_m\n')
            else:
                if u'父' in relation:
                    trio_cmd.write('python /DATA/sslyu/trio_BWA-GATK_3.0/annotation/annotation_filt_ver2.8.py -sn $sample -c_f_m c_f\n')
                else:
                    trio_cmd.write('python /DATA/sslyu/trio_BWA-GATK_3.0/annotation/annotation_filt_ver2.8.py -sn $sample -c_f_m c_m\n')
            trio_cmd.write('fi\n')
            trio_cmd.close()

            trio_k_ID = []
            for i in pedigree[k]:
                trio_k_ID.append(i)
            for i in trio_k_ID:
                list_trio.write(str(i)+'\n')
                generate_cmd(str(i))
                generate_run2_sentieon(str(i),cmd_trio)
        list_trio.close()
        list.close()

os.system('sort list |uniq > tmp; mv tmp list')
if not os.path.exists('trio/peddy'):
    os.makedirs('trio/peddy')
if not os.path.exists('trio/gvcf'):
    os.makedirs('trio/gvcf')
os.system('cp /DATA/sslyu/trio_BWA-GATK_3.0/trio/{cmd.sh,trio_peddy.sh} trio/peddy')    
if os.path.exists('peddy_tmp'):
    os.system('rm peddy_tmp')
if df.shape[0] == 1:
    if strategy == 'WES':
        os.system('cp /DATA/sslyu/trio_BWA-GATK_3.0/single_sex_check/B180625105.g.vcf* trio/gvcf')
        os.system('cp /DATA/sslyu/trio_BWA-GATK_3.0/single_sex_check/B180625105.mendel.ped trio/ped')
    if strategy == 'IDT-PANEL':
        os.system('cp /DATA/sslyu/trio_BWA-GATK_3.0/single_sex_check/18P3855.g.vcf* trio/gvcf')
        os.system('cp /DATA/sslyu/trio_BWA-GATK_3.0/single_sex_check/18P3855.mendel.ped trio/ped')
os.system('for i in `ls trio/ped/*mendel.ped`; do cat $i >> peddy_tmp; done')
os.system('awk \'{OFS="\t"; print \'0\',$2,$3,$4,$5,$6}\' peddy_tmp > trio/ped/peddy.ped')
if os.path.exists('trio/list'):
    list = open('trio/list', 'a')
    list.write('peddy\n')                    
    list.close()
else:
    list = open('trio/list', 'w')
    list.write('peddy\n')
    list.close()
os.system('sort trio/list |uniq > tmp; mv tmp trio/list')
with open('trio/ped/peddy.ped') as f:
    for l in f:
        print l.strip('\n')



    






