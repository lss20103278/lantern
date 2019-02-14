#!/usr/bin/env python
# -*- coding: utf-8 -*-
#########################################################################
# File Name: lib.py
# Author: sslyu
# mail: 1573077420@qq.com
# Created Time:Sat 29 Sep 2018 04:18:00 PM CST
#########################################################################

import numpy as np
import pandas as pd
import sys
reload(sys)
sys.setdefaultencoding("utf8")
import os
import re
import collections

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

def generate_pedigree(df):
    pedigree = {}
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
    return pedigree

def append_relation(df):
    df.index = df[u'姓名']
    df['relationship'] = u'子'
    pedigree_dict = generate_pedigree(df)
    trio = pedigree_dict.keys()
    for i in df.index:
        for j in trio:
            if j in i:
                if j == i:
                    df.loc[i,'relationship'] = u'子'
                elif u'父' in i:
                    df.loc[i,'relationship'] = u'父'  #there is no u'父' in i
                elif u'母' in  i:
                    df.loc[i,'relationship'] = u'母'  #there is no u'母' in i
                else:
                    df.loc[i,'relationship'] = 'other'
    return df

def append_pedigree(df):
    pedigree_dict = generate_pedigree(df)
    df.index = df[u'姓名']
    df['familyname'] = df['sample']
    for i in df.index:
        for k in pedigree_dict:
            if k in i: # the name of the child should be totally included in other family member's names
                df.loc[i]['familyname'] = df.loc[k]['sample']
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
    print pedigree_dict.values()
    trio = pedigree_dict.keys()
    try:
        single = [i for i in df[u'姓名'] if i not in flatten_list(pedigree_dict.values())] # dict.values()
    except:
        #single = df[u'姓名'].tolist()
        pass
    #for i in single:
    #    print unicode(i, 'utf-8')
    df['phenotype1'] = '0'
    df['phenotype2'] = None
    try:
        for i in df.index:
            if i in trio:
                df.loc[i,'phenotype2'] = '2'
            elif i in single:
                df.loc[i,'phenotype2'] = '2'
            else:
                df.loc[i,'phenotype2'] = '1'
    except:
        for i in df.index:
            if i in trio:
                df.loc[i,'phenotype2'] = '2'
            else:
                df.loc[i,'phenotype2'] = '1'
    return df

def append_father_mother(df):
    df['father'] = '0'
    df['mother'] = '0'
    pedigree_dict = generate_pedigree(df)
    trio = pedigree_dict.keys()
    try:
        single = [i for i in df[u'姓名'] if i not in pedigree_dict.values()[0]]
    except:
        single = df[u'姓名'].tolist()
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
                        df.loc[j]['father'] = df.loc[father_name]['sample']
                        df.loc[j]['mother'] = df.loc[mother_name]['sample']
                    else:
                        if len(pedigree_dict[i]) == 3:
                            df.loc[j]['father'] = '0'
                            df.loc[j]['mother'] = '0'
                        else:
                            if u'姐' in j or u'妹' in j or u'哥' in j or u'弟' in j:
                                for k in pedigree_dict[i]:
                                    if u'父' in k:
                                        father_name = k
                                    if u'母' in k:
                                        mother_name = k
                                df.loc[j]['father'] = df.loc[father_name]['sample']
                                df.loc[j]['mother'] = df.loc[mother_name]['sample']
                            elif u'爷' not in relation and u'奶' not in relation and u'外婆' not in relation and u'外公' not in relation:
                                if u'父' in j:
                                    df.loc[j]['father'] = '0'
                                    df.loc[j]['mother'] = '0'
                                if u'母' in j:
                                    df.loc[j]['father'] = '0'
                                    df.loc[j]['mother'] = '0'
                            else:
                                father = raw_input('father ID of '+j+'(if not exist, please enter 0):')
                                try:
                                    df.loc[j]['father'] = father
                                except:
                                    df.loc[j]['father'] = '0'
                                mother = raw_input('mother ID of '+j+'(if not exist, please enter 0):')
                                try:
                                    df.loc[j]['mother'] = mother
                                except:
                                    df.loc[j]['mother'] = '0'
                else:
                    if j == i:
                        if u'父' in relation:
                            for k in pedigree_dict[i]:
                                if u'父' in k:
                                    father_name = k
                            df.loc[j]['father'] = df.loc[father_name]['sample']
                        if u'母' in relation:
                            for k in pedigree_dict[i]:
                                if u'母' in k:
                                    mother_name = k
                            df.loc[j]['mother'] = df.loc[mother_name]['sample']
                    elif u'父' in j:
                        if u'爷' not in relation and u'奶' not in relation and u'外婆' not in relation and u'外公' not in relation:
                            df.loc[j]['father'] = '0'
                            df.loc[j]['mother'] = '0'
                    elif u'母' in j:
                        if u'爷' not in relation and u'奶' not in relation and u'外婆' not in relation and u'外公' not in relation:
                            df.loc[j]['father'] = '0'
                            df.loc[j]['mother'] = '0'
                    else:
                        father = raw_input('father ID of '+j+'(if not exist, please enter 0):')
                        try:
                            df.loc[j]['father'] = father
                        except:
                            df.loc[j]['father'] = '0'
                        mother = raw_input('mother ID of '+j+'(if not exist, please enter 0):')
                        try:
                            df.loc[j]['mother'] = mother
                        except:
                            df.loc[j]['mother'] = '0'
    return df                            

