#!/usr/bin/env python
# -*- coding: utf-8 -*-
#########################################################################
# File Name: format_web.py
# Author: sslyu
# mail: 1573077420@qq.com
# Created Time:Sat 15 Sep 2018 10:51:06 AM CST
# Usage: This script is for transforming the annotation excel to the appropriate form for the online system
#########################################################################

import numpy as np
import pandas as pd
import sys
reload(sys)
sys.setdefaultencoding("utf8")
import os
import re

infile = sys.argv[1]
outfile = sys.argv[2]

d = pd.ExcelFile('/DATA/sslyu/trio_BWA-GATK_3.0/annotation/示例文件.xlsx')
standard = d.parse(u'Sheet2', header=None)

d1 = pd.ExcelFile(infile)
selected_mutations = d1.parse(u'selected_mutations')
standard_colname = standard[0]

selected_mutations.rename(columns={'avsnp150':'avsnp147', 'gnomAD_genome_ALL':'gnomAD_exome_ALL', 'gnomAD_genome_EAS':'gnomAD_exome_EAS','CLNSIG':'CLINSIG','CLNDN':'CLNDBN','CLNALLELEID':'CLNACC','CLNDISDB':'CLNDSDB','CLNREVSTAT':'CLNDSDBID','GTinfo':'Gtinfo','InterVar_automated':'InterVar(automated)'}, inplace=True)
selected_mutations_columns = selected_mutations.columns.tolist()
index = []
for i in standard_colname:
    if i in selected_mutations_columns:
        index.append(selected_mutations_columns.index(i))
    else:
        index.append(i)
adjust_columns = []
for i in index:
    adjust_columns.append(selected_mutations_columns[i])
out = selected_mutations[adjust_columns]
out.to_excel(outfile,sheet_name="selected_mutations",index=False,header=True)
