#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
version:2.6
by:lss
"""

import numpy as np
import pandas as pd
import sys
reload(sys)
sys.setdefaultencoding("utf8")
import os
import re
import collections

def append_sample_excel(df):
    ID = []
    if u'样本编号' in df.columns and u'原始样本ID' in df.columns:
        for i in range(len(df[u'原始样本ID'])):
            if df[u'原始样本ID'][i] == 'nan':
                ID.append(str(df.iloc[i][u'样本编号']))
            else:
                ID.append(str(df.iloc[i][u'原始样本ID']))
    if u'样本ID' in df.columns and u'原始编码' in df.columns:
        for i in range(len(df[u'原始编码'])):
            if df[u'原始编码'][i] == 'nan':
                ID.append(str(df.iloc[i][u'样本ID']))
            else:
                ID.append(str(df.iloc[i][u'原始编码']))
    df['sample'] = ID
    return df

def append_sample_txt(df):
    ID = []
    if u'样本编号' in df.columns and u'原始样本ID' in df.columns:
        for i in range(len(df[u'原始样本ID'].isnull())):
            if df[u'原始样本ID'].isnull()[i]:                
                ID.append(df.iloc[i][u'样本编号'])
            else:
                ID.append(df.iloc[i][u'原始样本ID'])
    if u'样本ID' in df.columns and u'原始编码' in df.columns:
        for i in range(len(df[u'原始编码'].isnull())):
            if df[u'原始编码'].isnull()[i]:
                ID.append(df.iloc[i][u'样本ID'])
            else:
                ID.append(df.iloc[i][u'原始编码'])
    df['sample'] = ID
    return df

def append_gender(df):
    df['gender'] = df[u'性别'].apply(lambda x:'1' if u'男' in x else '2' if u'女' in x else 'unknown')
    return df

def generate_pedigree(df):
    pedigree = {}
    if 'familyname' not in df.columns:
        df.index = df[u'姓名']
        for i in df[u'姓名']:
            pedigree[i] = []
        for k in pedigree:
            for i in df[u'姓名']:
                if k in i: # the name of the child should be totally included in other family member's names
                    pedigree[k].append(i)
        list_keys = pedigree.keys()
        for k in list_keys:
            if len(pedigree[k]) == 1:
                pedigree.pop(k)
    elif 'familyname' in df.columns:
        import collections
        df.index = df['sample']
        trio = [item for item, count in collections.Counter(df['familyname'].tolist()).items() if count > 1]
        for i in trio:
            pedigree[df.loc[i][u'姓名']] = []
            for j in df['sample']:
                if df.loc[j]['familyname'] == i:
                    pedigree[df.loc[i][u'姓名']].append(df.loc[j][u'姓名'])
    df.index = df[u'姓名']
    return pedigree

def append_relation(df):
    df.index = df[u'姓名']
    if 'relationship' not in df.columns:
        df['relationship'] = None
    pedigree_dict = generate_pedigree(df)
    trio = pedigree_dict.keys()
    for i in df.index:
        for j in trio:
            if j in i:
                if j == i:
                    df.loc[i,'relationship'] = u'子'
                #elif re.match(j+u'父',i):    # the regular expression has a bug
                #    df.loc[i]['relationship'] == u'父'
                #elif re.match(j+u'母',i):
                #    df.loc[i]['relationship'] == u'母'
                #elif j+u'父' == i:
                #    df.loc[i,'relationship'] = u'父'
                #elif j+u'母' ==  i:
                #    df.loc[i,'relationship'] = u'母'
                elif u'父' in i:
                    df.loc[i,'relationship'] = u'父'  #there is no u'父' in i
                elif u'母' in  i:
                    df.loc[i,'relationship'] = u'母'  #there is no u'母' in i
                else:
                    df.loc[i,'relationship'] = 'other'
    df[u'家系关系'] = None
    df[u'家系关系'] = df['relationship']
    return df

def append_pedigree(df):
    pedigree_dict = generate_pedigree(df)
    df.index = df[u'姓名']
    if 'familyname' not in df.columns:
        df['familyname'] = df['sample']
        for i in df.index:
            for k in pedigree_dict:
                if k in i:  # the name of the child should be totally included in other family member's names
                    df.loc[i,'familyname'] = df.loc[k]['sample']
        df[u'样本间关系'] = None
        for i in df.index:
            for k in pedigree_dict:
                if k in i:
                    df.loc[i,u'样本间关系'] = 'fam'+df.loc[k]['sample']
        #df[u'样本间关系'] = df['familyname']
    else:
        df[u'样本间关系'] = None
        for i in df.index:
            for k in pedigree_dict:
                if df.loc[i]['pedigree'] == df.loc[k]['pedigree']:
                    df.loc[i,u'样本间关系'] = 'fam'+df.loc[k]['sample']
    return df                    

def flatten_list(nested):
    if isinstance(nested, list):
        for sublist in nested:
            for item in flatten_list(sublist):
                yield item
    else:
        yield nested

def append_phenotype(df):
    pedigree_dict = generate_pedigree(df)
    trio = pedigree_dict.keys()
    try:
        single = [i for i in df[u'姓名'] if i not in flatten_list(pedigree_dict.values())] # dict.values()
    except:
        pass
    df['phenotype1'] = '0'
    df['phenotype2'] = None
    #try:
    #    for i in df.index:
    #        if i in trio:
    #            df.loc[i,'phenotype2'] = '2'
    #        elif i in single:
    #            df.loc[i,'phenotype2'] = '2'
    #        else:
    #            df.loc[i,'phenotype2'] = '1'
    #except:
    #    for i in df.index:
    #        if i in trio:
    #            df.loc[i,'phenotype2'] = '2'
    #        else:
    #            df.loc[i,'phenotype2'] = '1'
    if len(trio) > 0 and len(single) > 0:
        for i in df.index: 
            if i in trio:
                df.loc[i,'phenotype2'] = '2'
            elif i in single:
                df.loc[i,'phenotype2'] = '2'
            else:
                df.loc[i,'phenotype2'] = '1'
    elif len(trio) > 0:
        for i in df.index:
            if i in trio:
                df.loc[i,'phenotype2'] = '2'
            else:
                df.loc[i,'phenotype2'] = '1'
    else:
        for i in df.index:
            df.loc[i,'phenotype2'] = '2'
    return df

def generate_trio(pedigree_dict):
    trio = []
    for k in pedigree_dict:
        trio.append(k)
    return trio
def generate_single(df):
    pedigree_dict = generate_pedigree(df)
    trio = generate_trio(pedigree_dict)
    trio_m = []
    for i in df[u'姓名']:
        for k in pedigree_dict:
            trio_m.extend(pedigree_dict[k])
    single = [i for i in df[u'姓名'] if i not in trio_m]
    return single

def append_father_mother(df):
    df['father'] = '0'
    df['mother'] = '0'
    pedigree_dict = generate_pedigree(df)
    trio = pedigree_dict.keys()
    try:
        single = [i for i in df[u'姓名'] if i not in pedigree_dict.values()[0]]
    except:
        single = df[u'姓名'].tolist()
    df.index = df[u'姓名']
    for i in trio:
        relation = ''
        for j in pedigree_dict[i]:
            relation = relation+df.loc[j]['relationship']
        for j in df.index:
            if j in pedigree_dict[i]:
                if u'父' in relation and u'母' in relation:
                    if j == i:
                        for k in pedigree_dict[i]:
                            if u'父' in k:
                                father_name = k
                            if u'母' in k:
                                mother_name = k
                        df.loc[j,'father'] = df.loc[father_name]['sample']
                        df.loc[j,'mother'] = df.loc[mother_name]['sample']
                    else:
                        if len(pedigree_dict[i]) == 3:
                            df.loc[j,'father'] = '0'
                            df.loc[j,'mother'] = '0'
                        else:
                            if u'姐' in j or u'妹' in j or u'哥' in j or u'弟' in j:
                                for k in pedigree_dict[i]:
                                    if u'父' in k:
                                        father_name = k
                                    if u'母' in k:
                                        mother_name = k
                                df.loc[j,'father'] = df.loc[father_name]['sample']
                                df.loc[j,'mother'] = df.loc[mother_name]['sample']
                            elif u'爷' not in relation and u'奶' not in relation and u'外婆' not in relation and u'外公' not in relation:
                                if u'父' in j:
                                    df.loc[j,'father'] = '0'
                                    df.loc[j,'mother'] = '0'
                                if u'母' in j:
                                    df.loc[j,'father'] = '0'
                                    df.loc[j,'mother'] = '0'
                            else:
                                father = raw_input('father ID of '+j+'(if not exist, please enter 0): ')
                                df.loc[j,'father'] = father
                                #try:
                                #    df.loc[j]['father'] = father
                                #except:
                                #    df.loc[j]['father'] = '0'
                                mother = raw_input('mother ID of '+j+'(if not exist, please enter 0):')
                                df.loc[j,'mother'] = mother
                                #try:
                                #    df.loc[j]['mother'] = mother
                                #except:
                                #    df.loc[j]['mother'] = '0'
                else:
                    if j == i:
                        if u'父' in relation:
                            for k in pedigree_dict[i]:
                                if u'父' in k:
                                    father_name = k
                            df.loc[j,'father'] = df.loc[father_name]['sample']
                        if u'母' in relation:
                            for k in pedigree_dict[i]:
                                if u'母' in k:
                                    mother_name = k
                            df.loc[j,'mother'] = df.loc[mother_name]['sample']
                    elif u'父' in j:
                        if u'爷' not in relation and u'奶' not in relation and u'外婆' not in relation and u'外公' not in relation:
                            df.loc[j,'father'] = '0'
                            df.loc[j,'mother'] = '0'
                    elif u'母' in j:
                        if u'爷' not in relation and u'奶' not in relation and u'外婆' not in relation and u'外公' not in relation:
                            df.loc[j,'father'] = '0'
                            df.loc[j,'mother'] = '0'
                    else:
                        father = raw_input('father ID of '+j+'(if not exist, please enter 0):')
                        df.loc[j,'father'] = father
                        mother = raw_input('mother ID of '+j+'(if not exist, please enter 0):')
                        df.loc[j,'mother'] = mother
    return df                        

def generate_sinle_ped(k,d):
    d.index = d['sample']
    ped_mendel = d.loc[k][['sample', 'father', 'mother', 'gender', 'phenotype2']]
    fname = d.loc[k]['sample']+'/'+d.loc[k]['sample']+'.mendel.ped'
    fname = 'ped/'+d.loc[k]['sample']+'.mendel.ped'
    ped = open(fname, 'w')
    ped.write(str(d.loc[k]['sample'])+'\t'+'\t'.join(ped_mendel)+'\n')
    ped.close()
    d.index = d[u'姓名']
def generate_ped(k,d):
    pedigree_dict = generate_pedigree(d)
    ped_f = d.loc[pedigree_dict[k]][['sample', 'father', 'mother', 'gender', 'phenotype1']]
    ped_f.index = ped_f['sample']
    ped_mendel = d.loc[pedigree_dict[k]][['sample', 'father', 'mother', 'gender', 'phenotype2']]
    ped_mendel.index = ped_mendel['sample']
    f1name = 'trio/ped/'+d.loc[k]['sample']+'.ped'
    f2name = 'trio/ped/'+d.loc[k]['sample']+'.mendel.ped'
    f1name = 'ped/'+d.loc[k]['sample']+'.ped'
    f2name = 'ped/'+d.loc[k]['sample']+'.mendel.ped'
    ped1 = open(f1name, 'w')
    ped2 = open(f2name, 'w')
    relation = []
    for i in pedigree_dict[k]:
        relation.append(d.loc[i]['relationship'])
    a = ' '.join(relation)
    key1 = d.loc[k]['sample']
    ped1.write(key1+'\t'+'\t'.join(ped_f.loc[key1])+'\n')
    ped2.write(key1+'\t'+'\t'.join(ped_mendel.loc[key1])+'\n')
    if u'父' in a and u'母' in a:
        #for i in pedigree_dict[k]:
        #    for j in [u'父', u'母']:
        for j in [u'父', u'母']:
            for i in pedigree_dict[k]:
                if i != k and j in i:
                    key2 = d.loc[i]['sample']
                    ped1.write(key1+'\t'+'\t'.join(ped_f.loc[key2])+'\n')
                    ped2.write(key1+'\t'+'\t'.join(ped_mendel.loc[key2])+'\n')
    	n_c_f_m = [i for i in pedigree_dict[k] if u'父' not in i and u'母' not in i and i != k]
        for i in n_c_f_m:
            key2 = d.loc[i]['sample']
            ped1.write(key1+'\t'+'\t'.join(ped_f.loc[key2])+'\n')
            ped2.write(key1+'\t'+'\t'.join(ped_mendel.loc[key2])+'\n')
    elif u'父' in a or u'母' in a:
        parent = ''
        if u'父' in a:
            parent = u'父'
        else:
            parent = u'母'
        for i in pedigree_dict[k]:
            if parent in i and i != k:
                key2 = d.loc[i]['sample']
                ped1.write(key1+'\t'+'\t'.join(ped_f.loc[key2])+'\n')
                ped2.write(key1+'\t'+'\t'.join(ped_mendel.loc[key2])+'\n')
        n_c_f_m = [i for i in pedigree_dict[k] if parent not in i and i != k]
        for i in n_c_f_m:
            key2 = d.loc[i]['sample']
            ped1.write(key1+'\t'+'\t'.join(ped_f.loc[key2])+'\n')
            ped2.write(key1+'\t'+'\t'.join(ped_mendel.loc[key2])+'\n')
    else:
    	n_c = [i for i in pedigree_dict[k] if i != k]
    	for i in n_c:
    	    key2 = d.loc[i]['sample']
    	    ped1.write(key1+'\t'+'\t'.join(ped_f.loc[key2])+'\n')
    	    ped2.write(key1+'\t'+'\t'.join(ped_mendel.loc[key2])+'\n')
    ped1.close()
    ped2.close()

def generate_peddy(trio,d):
    peddy = open('trio/ped/peddy.ped', 'a')
    peddy_ID = '0'
    for i in trio:
        trio_ID = d.loc[i]['sample']
        with open('trio/ped/'+trio_ID+'.mendel.ped') as f:
            for l in f:
                l = l.strip('\n').split(' ')
                peddy.write(peddy_ID+'\t'+'\t'.join(l[1:])+'\n')
    peddy.close()

def generate_trio_list(trio, df):
    child = []
    for i in trio:
        child.append(df.loc[i]['sample'])
    ofile = open('trio/list', 'w')
    for i in child:
    	ofile.write(i+'\n')
    #ofile.write('peddy\n')
    ofile.close()

def generate_cmd(path):
    os.system('cp /DATA/sslyu/trio_BWA-GATK_3.0/src/cmd.sh '+path)

def generate_single_cmd_step1(path):
    single_cmd_step1 = open(path+'/cmd_step1.sh', 'w')
    single_cmd_step1.write('#!/bin/sh\n\n#SBATCH -J jobname\n#SBATCH -p BP10\n#SBATCH -N 1\n#SBATCH -n 4\n\nsample=list\n\nsh run2_sentieon_step1.sh $sample $SLURM_NPROCS\n')
    single_cmd_step1.write('sed -n \'2p\' peddy.sex_check.csv |awk -F\',\' \'{if($1=="True"){print "sex of '+path+' is wrong"}}\' >> ../sex_check_log\n')
    single_cmd_step1.close()
def generate_single_cmd_step2(path):
    single_cmd_step2 = open(path+'/cmd_step2.sh', 'w')
    single_cmd_step2.write('#!/bin/sh\n\n#SBATCH -J jobname\n#SBATCH -p BP10\n#SBATCH -N 1\n#SBATCH -n 4\n\nsample=list\n\nsh run2_sentieon_step2.sh $sample $SLURM_NPROCS\n')
    single_cmd_step2.write('cp 2_mapping/$sample.depth.sample_gene_summary ../annotation\n')
    single_cmd_step2.write('cp 3_variants/$sample.{ann.hg19_multianno.txt,CADD,maf,link} ../annotation\n')
    single_cmd_step2.write('cd ../annotation\n')
    single_cmd_step2.write('echo '+path+' >> list\n')
    single_cmd_step2.write('python /DATA/sslyu/trio_BWA-GATK_3.0/src/score_re.py -sn $sample\n')
    single_cmd_step2.write('python /DATA/sslyu/trio_BWA-GATK_3.0/src/annotation_filt_ver3.0.py -sn $sample\n')
    single_cmd_step2.close()

def generate_single_cmd(path):
    single_cmd = open(path+'/cmd.sh', 'w')
    with open('/DATA/sslyu/trio_BWA-GATK_3.0/src/cmd.sh') as f:
        for l in f:
            single_cmd.write(l)
    single_cmd.write('cp 2_mapping/$sample.depth.sample_gene_summary ../annotation\n')
    single_cmd.write('cp 3_variants/$sample.{ann.hg19_multianno.txt,CADD,maf,link} ../annotation\n')
    single_cmd.write('cd ../annotation\n')
    single_cmd.write('echo $sample >> list\n')
    single_cmd.write('python /DATA/sslyu/trio_BWA-GATK_3.0/src/score_re.py -sn $sample\n')
    single_cmd.write('python /DATA/sslyu/trio_BWA-GATK_3.0/src/annotation_filt_ver3.0.py -sn $sample -c_f_m c\n')
    single_cmd.write('echo "end `date`" >> ../$sample/finished\n')
    if not os.path.exists('annotation'):
        os.mkdir('annotation')
    ofile = open('annotation/annotation.sh', 'a')
    ofile.write('python /DATA/sslyu/trio_BWA-GATK_3.0/src/annotation_filt_ver3.0.py -sn '+path+' -c_f_m c\n')
    ofile.close()
    single_cmd.close()
	
def generate_trio_cmd(path,filtmode):
    trio_cmd = open(path+'/trio_cmd.sh', 'w')
    n = 0
    with open('/DATA/sslyu/trio_BWA-GATK_3.0/src/trio_cmd.sh') as f:
        for l in f:
            n = n+1
            if n<=22:
                if n == 12:
                    trio_cmd.write('filtmode=\''+filtmode+'\'\n')
                else:
                    trio_cmd.write(l)
    return trio_cmd

def generate_trio_cmd_merge(path,filtmode):
    trio_cmd_merge = open(path+'/trio_cmd_merge.sh', 'w')
    n = 0
    with open('/DATA/sslyu/trio_BWA-GATK_3.0/src/trio_cmd.sh') as f:
        for l in f:
            n = n+1
            if n<=22:
                if n == 12:
                    trio_cmd_merge.write('filtmode=\''+filtmode+'\'\n')
                elif n == 20:
                    trio_cmd_merge.write('sh trio_peddy_merge.sh $sample $SLURM_NPROCS\n')
                elif n == 22:
                    trio_cmd_merge.write('sh trio_merge.sh $sample $SLURM_NPROCS $filtmode\n')
                else:
                    trio_cmd_merge.write(l)
    return trio_cmd_merge

def generate_dbevn(strategy):
    dbevn = open('dbevn.sh', 'w')
    if strategy != 'none':
        with open('/DATA/sslyu/trio_BWA-GATK_3.0/src/dbevn.sh') as f:
            for l in f:
                if l.startswith('#panel'):
                    if strategy in l:
                        dbevn.write(l[1:])
                    else:
                        dbevn.write(l)
                else:
                	dbevn.write(l)
    else:
        with open('/DATA/sslyu/trio_BWA-GATK_3.0/src/dbevn.sh') as f:
            for l in f:
                dbevn.write(l)
        panel = raw_input('please enter the absolute path of the bed:')
        dbevn.write('panel=\"'+panel+'\"\n')
    dbevn.close()

def generate_run2_sentieon(path,str):
	run2_sentieon = open(path+'/run2_sentieon.sh', 'w')
	n = 0
	with open('/DATA/sslyu/trio_BWA-GATK_3.0/src/run2_sentieon.sh') as f:
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

def generate_run2_sentieon_step1(path,str):
	run2_sentieon_step1 = open(path+'/run2_sentieon_step1.sh', 'w')
	n = 0
	with open('/DATA/sslyu/trio_BWA-GATK_3.0/src/run2_sentieon.sh') as f:
		for l in f:
			n = n+1
			if n>0 and n<= 23:
				run2_sentieon_step1.write(l)
			elif n>= 25:
				run2_sentieon_step1.write(l)
			elif n==24:
				run2_sentieon_step1.write('cmd="'+str+'"\n')
			else:
				continue
	run2_sentieon_step1.close()
def generate_run2_sentieon_step2(path,str):
	run2_sentieon_step2 = open(path+'/run2_sentieon_step2.sh', 'w')
	n = 0
	with open('/DATA/sslyu/trio_BWA-GATK_3.0/src/run2_sentieon.sh') as f:
		for l in f:
			n = n+1
			if n>0 and n<= 23:
				run2_sentieon_step2.write(l)
			elif n>= 25:
				run2_sentieon_step2.write(l)
			elif n==24:
				run2_sentieon_step2.write('cmd="'+str+'"\n')
			else:
				continue
	run2_sentieon_step2.close()

#def generate_annotation(df):
#    pedigree_dict = generate_pedigree(df)
#    trio = generate_trio(pedigree_dict)
#    single = generate_single(df)
#    if not os.path.exists('annotation'):
#        os.makedirs('annotation')
#    ofile = open('step3_annotation.sh', 'a')
#    if len(trio) > 0:
#        ofile1 = open('annotation/list', 'a')
#        ofile2 = open('annotation/list_n_c_f_m', 'a')
#        for i in trio:
#            ID = df.loc[i]['sample']
#            ofile.write('mv '+ID+'/2_mapping/'+ID+'.depth.sample_gene_summary annotation\n')
#            ofile.write('mv trio/'+ID+'/segtrio/'+ID+'.{ann*,CADD,link,maf} annotation\n')
#            relation = []
#            for j in pedigree_dict[i]:
#                relation.append(df.loc[j]['relationship'])
#            a = ' '.join(relation)
#            if u'父' in a and u'母' in a:
#                ofile1.write(ID+'\n')
#            else:
#                ofile2.write(ID+'\n')
#        ofile1.close()
#        ofile2.close()
#    if os.path.exists('annotation/list_n_c_f_m'):
#        if os.path.getsize('annotation/list_n_c_f_m') == 0:
#            os.system('rm annotation/list_n_c_f_m')
#    if os.path.exists('annotation/list'):
#        if os.path.getsize('annotation/list') == 0:
#            os.system('rm annotation/list')
#    if len(single) > 0:
#        if os.path.exist('annotation/list_n_c_f_m'):
#            ofile3 = open('annotation/list_n_c_f_m', 'a')
#        else:
#            ofile3 = open('annotation/list', 'a')
#        for i in single:
#            ID = df.loc[i]['sample']
#            ofile.write('mv '+ID+'/2_mapping/'+ID+'.depth.sample_gene_summary annotation\n')
#            ofile.write('mv '+ID+'/3_variants/'+ID+'.{ann*,CADD,link,maf} annotation\n')
#            ofile3.write(ID+'\n')
#        ofile3.close()
#    ofile.write('cd annotation\n')
#    if os.path.exists('annotation/list_n_c_f_m'):
#        ofile.write('python /DATA/sslyu/trio_BWA-GATK_3.0/annotation/score_re.py -l list_n_c_f_m\n')
#        ofile.write('nohup python /DATA/sslyu/trio_BWA-GATK_3.0/annotation/annotation_filt_ver3.0.py -l list_n_c_f_m -c_f_m no &\n')
#    if os.path.exists('annotation/list'):
#        ofile.write('python /DATA/sslyu/trio_BWA-GATK_3.0/annotation/score_re.py -l list\n')
#        ofile.write('nohup python /DATA/sslyu/trio_BWA-GATK_3.0/annotation/annotation_filt_ver3.0.py -l list &\n')
#    ofile.write('cd ..\n')
#    ofile.close()

