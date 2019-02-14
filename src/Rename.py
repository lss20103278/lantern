import os
import pandas as pd
#with open('list_old') as f_name_old:
	#lst_name_old = map(lambda x:x.strip('\n'),f_name_old.readlines())
#with open('list') as f_name_new:
	#lst_name_new=map(lambda x:x.strip('\n'),f_name_new.readlines())
#dict_name = dict(zip(lst_name_old,lst_name_new))
df = pd.read_csv('rename',sep = '\t',header = None,index_col=0)

path = os.getcwd()
for file in os.listdir(path):
		for i in df.index:
			if file.startswith(i):
				os.rename(file,df.ix[i,1]+'_'+'_'.join(file.split('_')[1:]))
