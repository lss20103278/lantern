#!/opt/anaconda2/bin/python
import sys

import pandas as pd

sn = str(sys.argv[1])
df = pd.read_csv(sn+'.tmp.vcf',sep = '\t')
df_tmp = df.copy()

with open(sn+'.gvcf.list') as f:
	sample_order = f.readlines()
sample_order = map(lambda x:x.split('gvcf/')[1].split('.g.')[0],sample_order)
sample_num = len(sample_order)
columns = df.columns.tolist()
OrderInVCF = df.columns.tolist()[-sample_num:]
if sample_order != OrderInVCF:
	for i in range(sample_num):
		df[OrderInVCF[i]]=df_tmp[sample_order[i]]
	new_columns = columns[:-sample_num]
	new_columns.extend(sample_order)
	df.columns = new_columns
	df.to_csv(sn+'.sort.tmp.vcf',sep = '\t',index = False)


