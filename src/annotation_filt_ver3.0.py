#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
version:3.0
by:lss
inherited from ykl
update:2018/9/15
"""
ThisVer='v3.0'
Thisreleasetime='2018/9/15'

import os
import sys
import pandas as pd
import numpy as np


sys.path.append("/DATA/sslyu/trio_BWA-GATK_3.0/")
from lib.GAF import *
from lib.xlsx import *


##verion confirm
print "\nGenetic annotation, %s, update:%s, by:sslyu" % (ThisVer,Thisreleasetime)
##


##read parameters of cmd line
kwargs_raw = sys.argv[1:]
kwargs={'-edge':'ucsc_hg19_refGene.csv','-l':'list','-AF':"combine_run1.snp.AF",'-sn':'', '-c_f_m':'', '-thisminedge':'15'}
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
	list=pd.read_table(kwargs.get('-l'),header=None)[0].tolist()
	print "\nNote:list samples mode"

c_f_m = kwargs.get('-c_f_m') ##deal with trio that doesn't have a child_father_mother relatonship

#database confirm:
database = [kwargs.get('-edge'),kwargs.get('-AF'),"OMIM_selected_ch.csv","OMIM_selected_en.csv","phenotype_annotation_hpoteam_include_chpo_doid_xref.csv","morbidmap","DIDGene.csv","OMIM_Gene-Map-Search_20180329.csv","Phenotype20180329.csv"]
missing = []
for i in database:
    if not os.path.exists("/DATA/sslyu/trio_BWA-GATK_3.0/doc/"+i):
        missing.append(i)
if len(missing) > 0:
    print ','.join(missing)+' are missing'
    sys.exit()
else:
    print """
    Database confirming:
    -edge: /DATA/sslyu/trio_BWA-GATK_3.0/%s
    -list:%s
    -locMAF: /DATA/sslyu/trio_BWA-GATK_3.0/%s 
    OMIM_ch: /DATA/sslyu/trio_BWA-GATK_3.0/%s
    OMIM_en: /DATA/sslyu/trio_BWA-GATK_3.0/%s
    HPO: /DATA/sslyu/trio_BWA-GATK_3.0/%s
    OMIM_morbidmap: /DATA/sslyu/trio_BWA-GATK_3.0/%s
    DIDGene: /DATA/sslyu/trio_BWA-GATK_3.0/%s
    OMIM_Gene-Map-Search_20180329: /DATA/sslyu/trio_BWA-GATK_3.0/%s
    OMIM_Phenotype: /DATA/sslyu/trio_BWA-GATK_3.0/%s
    """ % ("doc/"+kwargs.get('-edge'),kwargs.get('-l'),"doc/"+kwargs.get('-AF'),"doc/OMIM_selected_ch.csv","doc/OMIM_selected_en.csv","doc/phenotype_annotation_hpoteam_include_chpo_doid_xref.csv","doc/morbidmap","doc/DIDGene.csv","doc/OMIM_Gene-Map-Search_20180329.csv","doc/Phenotype20180329.csv")

AF=pd.read_table("/DATA/sslyu/trio_BWA-GATK_3.0/doc/"+kwargs.get('-AF'),sep="\t",header = None,dtype = str)
AF = AF[[0,2]]
AF.columns = [0,'locAF']

def maf_process(sn, ds):
    df_maf = pd.read_csv(sn+'.maf',sep = '\t',header = 0,dtype = str)
    df_maf = df_maf[['Chromosome','Start_Position','End_Position','Reference_Allele','Tumor_Seq_Allele2','HGVSc','HGVSp']]
    df_maf.columns = ['Chr','Start','End','Ref','Alt','HGVSc','HGVSp'] # HGVSc and HGVSp
    df_maf_tmp = df_maf['Ref'].apply(lambda x:1 if x=='-' else 0)
    df_maf['End'] = df_maf['End'].astype(int)-df_maf_tmp.astype(int)
    df_maf['End']=df_maf['End'].astype(str)
    ds = pd.merge(ds,df_maf,how = 'left',on = ['Chr','Start','End','Ref','Alt'])
    return ds

def p_ID_l_g(ds):    
    p_ID_l = []
    for i in ds.columns:
    	if 'GTinfo_' in i:
    		p_ID_l.append(i[7:])
    return p_ID_l

def percent(ds,i,p_ID_l):            
    perc = []
    for j in ds['GTinfo_'+i]:
        try:
            c = j.split(':')[1].split(',')
            d = float(c[0])+float(c[1])
            perc.append('{:.2f}'.format(float(c[1])/d*100)+'%')
        except:
            perc.append(np.nan)
    return perc

def link_process(sn, ds):
    df_link = pd.read_csv(sn+'.link',sep = '\t',header = None,dtype = str)
    df_link[0] = df_link[0].apply(lambda x:x if x.startswith('chr') else 'chr'+x)
    df_link['link'] = df_link[0]+'-'+df_link[6]+'-'+df_link[8]+'-'+df_link[9]
    df_link = df_link.drop([6,8,9],axis = 1)		
    p_ID_l = p_ID_l_g(ds)   #discriminate single or pedigree
    if len(p_ID_l) == 0:
        p_ID_l.append('GTinfo_'+sn) # sn should be the same with the name in sn.link
        df_link = df_link[[0,1,2,3,4,14,'link']]
        df_link.columns = ['Chr','Start','End','Ref','Alt','GTinfo_'+sn,'link']
        perc = percent(df_link, sn, p_ID_l)
        ds.insert(10,sn+'_allele_percentage',perc)
        df_link = df_link.drop('GTinfo_'+sn,axis = 1)
    else:
        df_link = df_link[[0,1,2,3,4,'link']]
        df_link.columns = ['Chr','Start','End','Ref','Alt','link']
        for i in p_ID_l:
            perc = percent(ds,i,p_ID_l)
            ds.insert(10+p_ID_l.index(i),i+'_allele_percentage',perc)
    ds = pd.merge(ds,df_link,how = 'left',on = [u'Chr',u'Start',u'End',u'Ref',u'Alt'])
    return ds

def cadd_process(sn, ds):
    df_cadd = pd.read_csv(sn+'.CADD',sep = '\t',header = 0,dtype = str)
    df_cadd['link'] = df_cadd['#CHROM'].apply(lambda x:'chr'+x)+'-'+df_cadd['POS']+'-'+df_cadd['REF']+'-'+df_cadd['ALT']
    df_cadd = df_cadd.drop(['#CHROM','POS','REF','ALT'],axis = 1)
    df_cadd.columns = ['CADD_RawScore','CADD_PHRED','link']
    ds = pd.merge(ds,df_cadd,how = 'left',on = ['link'])
    return ds

def del_process(sn, ds):
####### deal with the deletion problem to see if the deletion of nuclei acids affects amnio acids #########
    ds_del = ds[ds['Alt']=='-']  #ds_del: Alt == '-' what situation?

    value = ds_del['Ref'].apply(lambda x:len(x))
    ds_del.loc[:,'conf'] = value  #ds_del['conf'] = ds_del['Ref'].apply(lambda x:len(x)) #SettingWithCopyWarning:Try using .loc[row_indexer,col_indexer] = value instead

    ds_del_length = ds_del[ds_del['conf']>6] # Alt == '-', Ref of length more than 6 nt ?
    ds_del_length = ds_del_length.fillna('-')
    
    ds_del_length_fs = ds_del_length[ds_del_length['AAChange.refGene'].str.contains('[\d]fs')] # number+fs; fs:frameshift; loci whose ref is more than 6 nt length contains frameshift
   
    ds_del_length_fs.loc[:,'Deletion_Problem'] = 'note' #ds_del_length_fs['Deletion_Problem'] = 'note' # add column 'Deletion_Problem' SettingWithCopyWarning:Try using .loc[row_indexer,col_indexer] = value instead

    ds_del_length = ds_del_length[~ds_del_length['AAChange.refGene'].str.contains('[\d]fs')] #  loci whose ref is more than 6 nt length doesn't contain frameshift
    ds_del_length = ds_del_length[ds_del_length['AAChange.refGene'].str.contains('p.*del')] # loci doesn't contain frameshipt but contains p.*del
    AA_num = map(lambda x:int(x[1])-int(x[0]),ds_del_length['AAChange.refGene'].apply(lambda x:x.split(',')[0].split('p.')[1].split('del')[0].split('_'))) # extract the number between p. and del in p.*del Amino acid number
    ds_del_length.index = range(len(ds_del_length))
    Deletion_Problem = []
    for i in ds_del_length.index:
    	if int(ds_del_length.ix[i,'conf']) != AA_num[i]*3 and int(ds_del_length.ix[i,'conf']) != (AA_num[i]+1)*3: # deletion affects the coding? 3 nucleic acids code 1 amino acid
    		Deletion_Problem.append('note')
    	else:
    		Deletion_Problem.append('-')
    ds_del_length['Deletion_Problem'] = pd.Series(Deletion_Problem,index = ds_del_length.index)
    ds_del_length_all = pd.concat([ds_del_length_fs,ds_del_length])
    ds_del_length_all = ds_del_length_all[['Chr','Start','End','Ref','Alt','Deletion_Problem']]
    ds = pd.merge(ds,ds_del_length_all,how = 'left',on = ['Chr','Start','End','Ref','Alt'])
    return ds

def prediction_process(sn, ds):
    if 'Dead_Zone' in ds.columns.tolist() and 'Problem_High' in ds.columns.tolist() and 'Problem_Low' in ds.columns.tolist():
        ds[['Dead_Zone','Problem_High','Problem_Low']]=ds[['Dead_Zone','Problem_High','Problem_Low']].replace('Name=Y','Y') #
        ds['HI_Predictions'] = ds['HI_Predictions'].fillna('-').apply(lambda x:x if x=='-' else x.split('=')[1]) # 
    
    flag = True
    for i in ['InterVar_automated', 'PVS1', 'PS2', 'PS1', 'PS3', 'PS4', 'PM1', 'PM2', 'PM3', 'PM4', 'PM5', 'PM6', 'PP1', 'PP2', 'PP3', 'PP4', 'PP5', 'BA1', 'BS1', 'BS2', 'BS4', 'BP1','BP2','BP3','BP4','BP5','BP6','BP7']:
        if i not in ds.columns.tolist():
            flag = False
            break
    if flag:
        ds.insert(ds.columns.tolist().index('InterVar_automated')+1,'PVS1/PS1/PS2/PS3/PS4/PM1/PM2/PM3/PM4/PM5/PM6/PP1/PP2/PP3/PP4/PP5/BA1/BS1/BS2/BS3/BS4/BP1/BP2/BP3/BP4/BP5/BP6/BP7',np.nan)  #ds.insert(ds.columns.tolist().index('InterVar(automated)')+1,'PVS1/PS1/PS2/PS3/PS4/PM1/PM2/PM3/PM4/PM5/PM6/PP1/PP2/PP3/PP4/PP5/BA1/BS1/BS2/BS3/BS4/BP1/BP2/BP3/BP4/BP5/BP6/BP7',np.nan)# annovar 20180416 InterVar(automated)->InterVar_automated  
        ds['PVS1/PS1/PS2/PS3/PS4/PM1/PM2/PM3/PM4/PM5/PM6/PP1/PP2/PP3/PP4/PP5/BA1/BS1/BS2/BS3/BS4/BP1/BP2/BP3/BP4/BP5/BP6/BP7'] = \
    	ds['PVS1'].apply(lambda x:str(x))+'/'+ds['PS1'].apply(lambda x:str(x))+'/'+\
    	ds['PS2'].apply(lambda x: str(x)) + '/' + ds['PS3'].apply(lambda x: str(x)) + '/' + ds['PS4'].apply(lambda x: str(x)) + '/'+ \
    	ds['PM1'].apply(lambda x: str(x)) + '/' + ds['PM2'].apply(lambda x: str(x)) + '/' + ds['PM3'].apply(lambda x: str(x)) + '/'+ \
    	ds['PM4'].apply(lambda x: str(x)) + '/' + ds['PM5'].apply(lambda x: str(x)) + '/' + ds['PM6'].apply(lambda x: str(x)) + '/'+ \
    	ds['PP1'].apply(lambda x: str(x)) + '/' + ds['PP2'].apply(lambda x: str(x)) + '/' + ds['PP3'].apply(lambda x: str(x)) + '/'+ \
    	ds['PP4'].apply(lambda x: str(x)) + '/' + ds['PP5'].apply(lambda x: str(x)) + '/' + ds['BA1'].apply(lambda x: str(x)) + '/'+ \
    	ds['BS1'].apply(lambda x: str(x)) + '/' + ds['BS2'].apply(lambda x: str(x)) + '/' + ds['BS3'].apply(lambda x: str(x)) + '/'+ \
    	ds['BS4'].apply(lambda x: str(x)) + '/' + ds['BP1'].apply(lambda x: str(x)) + '/' + ds['BP2'].apply(lambda x: str(x)) + '/'+ \
    	ds['BP3'].apply(lambda x: str(x)) + '/' + ds['BP4'].apply(lambda x: str(x)) + '/' + ds['BP5'].apply(lambda x: str(x)) + '/'+ds['BP6'].apply(lambda x: str(x)) + '/' + ds['BP7'].apply(lambda x: str(x))
    if 'PVS1' in ds.columns.tolist() and 'CLINSIG' in ds.columns.tolist():
        ds = ds.drop(ds.columns.tolist()[ds.columns.tolist().index('PVS1'):ds.columns.tolist().index('CLINSIG')],axis=1)      
        #ds = ds.drop(ds.columns.tolist()[ds.columns.tolist().index('PVS1'):ds.columns.tolist().index('CLNALLELEID')],axis=1) # annovar 20180416 the order of clinvar columns changes

    ds['Start']=ds['Start'].astype(long)
    return ds

def check_chr(ds):
    #chr tag confirm
    if("chr" in ds.iloc[0]['Chr']):
    	print "chrid_check : ok"
    else:
    	ds['Chr']="chr"+ds['Chr']
    	print "chrid_rewrite : ok"
    """
    check=[]
    for i in ds.index:
    	line=ds.loc[i]['Gene.refGene'].split(";")
    	tmp=[id for id in line if id in genelist]
    	if(tmp):
    		check.append(True)
    	else:
    		check.append(False)
    
    ds=ds[check]
    """

def file_check(sn):
###### check if all files are appropriately prepared ###############
    for i in ["maf", "link", "CADD"]:
        if not os.path.exists(sn+'.'+i):
            print sn+'.'+i+" is missing"
            sys.exit()

def filt1(sn, ds, kwargs):
    ##step1:filt1+annotation_AF+annotation_OMIM
    print "\nStep1:\n--filting by snp/indel location"
    
    #read edge unkown genes:
    ex_edge=pd.read_csv("/DATA/sslyu/trio_BWA-GATK_3.0/doc/"+kwargs.get('-edge'))
    ex_edge=ex_edge.sort_values(['chrom'])
    
    #build edge dic:  used to accelerate below
    chroms=[]
    for i in range(22):
    	chroms.append("chr"+str(i+1))
    
    chroms.append('chrX')
    chroms.append('chrY')
    ex_init={}
    ex_end={}
    for i in chroms:
    	sub=ex_edge[ex_edge['chrom']==i]
    	ex_init.update({i:sub.index[0]})
    	ex_end.update({i:sub.index[-1]})   #end

    #determine the edge of exon
    if(kwargs.get('-thisminedge')!=''):
    	thisminedge_set = int(kwargs.get('-thisminedge'))
    	print "\nNote:the edge of exon has been set to "+kwargs.get('-thisminedge')

    ds.insert(7,'min_edge',[np.nan]*len(ds))
    tar=['exonic','splicing','exonic;splicing']
    ds_a=ds[ds['Func.refGene'].isin(tar)].copy() # .copy() modifications to the data or indices of the copy will not be reflected in the original object

    ds_b=ds[~ds['Func.refGene'].isin(tar)].copy() # Func.refGene is neither exonic nor splicing 
    ds_b.index=range(len(ds_b))
    for i in ds_b.index:
    	line=ds_b.loc[i]
    	chrom=line['Chr']
    	pos=int(line['Start'])
    	genes=line['Gene.refGene'].replace(",",";").split(";")
    	
    	#sub of edge by chr
    	""" bugs maybe: the end should +1? """
    	ex_sub=ex_edge.loc[ex_init.get(chrom):ex_end.get(chrom)] #??		
    	ex_sub=ex_edge[ex_edge['symbol'].isin(genes)] #??
    	refedge='::edge::'
    	#thisminedge=11
    	thisminedge=thisminedge_set
    	confheader=''
    	for j in ex_sub.index:						  
    		ex_subline=ex_sub.loc[j]     
    		eCount=ex_subline.exonCount
    		eStart=ex_subline.exonStarts.split()
    		for k in range(len(eStart)):eStart[k]=int(eStart[k])+1
    		eEnd=ex_subline.exonEnds.split()
    		for k in range(len(eEnd)):eEnd[k]=int(eEnd[k])
    		if(ex_subline.strand=='+'):
    			for k in range(eCount):
    				if(eStart[k] > pos and eStart[k]-pos <= thisminedge_set):
    					if(eStart[k]-pos < thisminedge_set+1): thisminedge=eStart[k]-pos
    					confheader='-' # in plus strand, '-' means it is located ahead of the start of the exon
    					refedge=refedge + ex_subline['name']+":"+ex_subline.symbol+":exon"+str(k+1)+":"+str(eStart[k])+":-"+str(eStart[k]-pos)+";"
    				if(eEnd[k] < pos and pos-eEnd[k] <= thisminedge_set):
    					if(pos-eEnd[k] < thisminedge_set+1): thisminedge=pos-eEnd[k]
    					confheader='+' # in plus strand, '+' means it is located behind the end of the exon
    					refedge=refedge + ex_subline['name']+":"+ex_subline.symbol+":exon"+str(k+1)+":"+str(eEnd[k])+":+"+str(pos-eEnd[k])+";"
    		elif(ex_subline.strand=='-'):
    			for k in range(eCount):
    				if(eEnd[k] < pos and pos-eEnd[k] <= thisminedge_set):
    					if(pos-eEnd[k] < thisminedge_set+1): thisminedge=pos-eEnd[k]
    					confheader='-' # in minus strand, '-' means it is located behind the end of the exon
    					refedge=refedge + ex_subline['name'] + ":" + ex_subline.symbol + ":exon" + str(eCount-k) +":"+str(eEnd[k])+ ":-" + str(pos-eEnd[k]) + ";"
    				if(eStart[k]>pos and eStart[k]-pos <= thisminedge_set):
    					if(eStart[k]-pos < thisminedge_set+1): thisminedge=eStart[k]-pos
    					confheader='+' # in minus strand, '+' means it is located ahead of the start of the exon
    					refedge=refedge + ex_subline['name'] + ":" + ex_subline.symbol + ":exon" + str(eCount-k) + ":"+str(eStart[k])+ ":+" + str(eStart[k]-pos) + ";"
    	ds_b.loc[i,'min_edge'] = confheader+str(thisminedge)
    	ds_b.loc[i,'GeneDetail.refGene'] = str(ds_b.loc[i]['GeneDetail.refGene']) + refedge.strip(";")		  
    	if(i%1000==0 or i==ds_b.index[-1]):
    	   view_bar(i+1,len(ds_b))
    conf=[]
    for i in ds_b.index:
    	if(ds_b.loc[i]['GeneDetail.refGene'].split("::edge::")[1]==''):
    		conf.append(False)
    	else:
    		conf.append(True)
    ds_c=ds_b[conf] # it is located in the min_edge of the exons
    ds=pd.concat([ds_a,ds_c])
    ds.index=range(len(ds))
    return ds

def add_AF(sn, ds, AF):    
    #annotation:local MAF
    print "\n--adding local AF"
    ds[0] = ds['Chr'].astype(str)+'-'+ds['Start'].astype(str)+'-'+ds['End'].astype(str)+'-'+ds['Ref'].astype(str)+'-'+ds['Alt'].astype(str)
    ds = pd.merge(ds,AF,how = 'left',on=[0])
    ds = ds.drop(0,axis = 1)
    return ds

def add_other(sn, ds):    
    #annotation:OMIM_ch/OMIM_en/chpo/morbidmap/Phenotype
    print "--adding OMIM_ch/OMIM_en/CHPO/morbidmap/Phenotype"
    gene_disease_9606=pd.read_csv("/DATA/sslyu/trio_BWA-GATK_3.0/doc/gene_disease.9606.tsv",encoding='gbk',dtype=str,sep="\t")
    gene_disease_9606=gene_disease_9606.fillna("NA")
    OMIM_selected_ch=pd.read_csv("/DATA/sslyu/trio_BWA-GATK_3.0/doc/OMIM_selected_ch.csv",encoding='gbk',dtype=str)
    OMIM_selected_en=pd.read_csv("/DATA/sslyu/trio_BWA-GATK_3.0/doc/OMIM_selected_en20180330.csv",encoding='gbk',dtype=str)
    phenon=pd.read_csv("/DATA/sslyu/trio_BWA-GATK_3.0/doc/phenotype_annotation_hpoteam_include_chpo_doid_xref.csv",encoding='gbk',dtype=str)
    phenon=phenon[phenon.infosource=='OMIM']
    phenon=phenon.fillna("NA")
    morbidmap = pd.read_csv("/DATA/sslyu/trio_BWA-GATK_3.0/doc/morbidmap",sep = '\t',names = ['Phenotype','Gene_Symbols','MIM_Number','Cyto_Location','Gene'],encoding = 'gbk',dtype = str) #update
    DIDGene = pd.read_csv('/DATA/sslyu/trio_BWA-GATK_3.0/doc/DIDGene.csv',sep = '\t',names = ['Gene'])
    OMIM_Gene_Map_Search = pd.read_csv('/DATA/sslyu/trio_BWA-GATK_3.0/doc/OMIM_Gene-Map-Search_20180329.csv',header = 0,encoding = 'gbk',dtype = str)
    OMIM_Phenotype = pd.read_csv('/DATA/sslyu/trio_BWA-GATK_3.0/doc/Phenotype20180329.csv',header = 0,encoding = 'gbk',dtype = str)
    OMIM_Phenotype = OMIM_Phenotype[['OMIM_id','Phenotype']]
    OMIM_Phenotype = OMIM_Phenotype.fillna('-')
    OMIM_Phenotype['OMIM_id'] = OMIM_Phenotype['OMIM_id'].apply(lambda x:x.strip())
    OMIM_Phenotype['Phenotype'] = OMIM_Phenotype['Phenotype'].apply(lambda x:x.strip())
    addgene_disease_9606=[]
    addOMIM_id=[]
    addOMIM_gene=[]
    addOMIM_disease_ch=[]
    addOMIM_disease_en=[]
    add_chpo_name_CH1=[]
    add_chpo_define_EN1=[]
    add_chpo_define_CH1=[]
    # update
    IP = []
    DI = []
    Inheritance = []
    Comments = []
    Gene_Locus_MIM_number = []
    Phenotype = []
    Phenotype_MIM_number = []
    add_mor_Gene = []
    add_mor_Phenotype = []
    add_mor_MIM_Number = []
    add_mor_Cyto_Location = []
    omim_phenotype = []
    ###########
    for i in ds.index:
        if(i%1000==0 or i==ds.index[-1]):
     	   view_bar(i+1,len(ds.index))
        
        targene=ds.loc[i]['Gene.refGene'].split(',')[0]
        # update 
        info_gene_disease_9606 = gene_disease_9606[gene_disease_9606.subject_label == targene]
        if len(info_gene_disease_9606) == 0:
            addgene_disease_9606.append('')
        else:
            addgene_disease_9606.append(','.join(info_gene_disease_9606['object_label'])+";  "+','.join(info_gene_disease_9606['relation_label'])+';  '+','.join(info_gene_disease_9606['source']))
        info_gene = DIDGene[DIDGene.Gene==targene]
        if(len(info_gene)==0):
     	   DI.append('')
        else:
     	   DI.append('YES')
        info_OMIM_Gene = OMIM_Gene_Map_Search[OMIM_Gene_Map_Search['Gene']==targene]
        if(len(info_OMIM_Gene)==0):
     	   Inheritance.append('')
     	   Comments.append('')
     	   Gene_Locus_MIM_number.append('')
     	   Phenotype.append('')
     	   Phenotype_MIM_number.append('')
        else:
     	   Inheritance.append(info_OMIM_Gene['Inheritance'].str.cat(sep = '|'))
     	   Comments.append(info_OMIM_Gene['Comments'].str.cat(sep = '|'))
     	   Gene_Locus_MIM_number.append(info_OMIM_Gene['Gene/Locus_MIM_number'].drop_duplicates().str.cat(sep = '|'))
     	   Phenotype.append(info_OMIM_Gene['Phenotype'].str.cat(sep = '||'))
     	   Phenotype_MIM_number.append(info_OMIM_Gene['Phenotype_MIM_number'].str.cat(sep='|'))		   
        info_mor = morbidmap[morbidmap.Gene == targene]
        if(len(info_mor) == 0):
     	   add_mor_Gene.append('')
     	   add_mor_Phenotype.append('')
     	   add_mor_MIM_Number.append('')
     	   add_mor_Cyto_Location.append('')
        else:
     	   add_mor_Gene.append(info_mor.iloc[0].Gene)
     	   add_mor_Phenotype.append(info_mor.iloc[0].Phenotype)
     	   add_mor_MIM_Number.append(info_mor.iloc[0].MIM_Number)
     	   add_mor_Cyto_Location.append(info_mor.iloc[0].Cyto_Location)
        ############################# 
        info_ch=OMIM_selected_ch[OMIM_selected_ch.OMIM_gene==targene]
        info_en=OMIM_selected_en[OMIM_selected_en.OMIM_gene==targene]
        if(len(info_en)==0):
     	   IP.append('')
     	   addOMIM_disease_en.append('')
        else:
     	   IP.append(info_en['IP'].str.cat(sep = '|'))
     	   addOMIM_disease_en.append(info_en['OMIM_disease'].str.cat(sep = '||'))		 
        if(len(info_ch)==0):
     	   addOMIM_id.append('')
     	   addOMIM_gene.append('')
     	   addOMIM_disease_ch.append('')
     	   add_chpo_name_CH1.append('')
     	   add_chpo_define_EN1.append('')
     	   add_chpo_define_CH1.append('')
        else:
     	   addOMIM_gene.append(info_ch.iloc[0].OMIM_gene)
     	   addOMIM_disease_ch.append(info_ch.iloc[0].OMIM_disease_ch)
     	   addOMIM_id.append(info_ch.iloc[0].OMIM_id)
     	   subphenon=phenon[phenon.infosourceid==info_ch.iloc[0].OMIM_id]
     	   if(len(subphenon)==0):
     		   add_chpo_name_CH1.append('')
     		   add_chpo_define_EN1.append('')
     		   add_chpo_define_CH1.append('')
     	   else:
     		   chpo_name=u";"
     		   chpo_define_EN=u";"
     		   chpo_define_CH=u";"
     		   for j in subphenon.index:
     			   chpo_name=chpo_name+subphenon.loc[j].chpo_name_CH1+u";"
     			   chpo_define_EN=chpo_define_EN + subphenon.loc[j].chpo_define_EN1 + u';'	
     			   chpo_define_CH=chpo_define_CH + subphenon.loc[j].chpo_define_CH1 + u";"
     		   add_chpo_name_CH1.append(chpo_name.strip(';'))
     		   add_chpo_define_EN1.append(chpo_define_EN.strip(";").strip("NA;"))
     		   add_chpo_define_CH1.append(chpo_define_CH.strip(";").strip("NA;"))
    
        #ename=''
        #for j in info_en.OMIM_disease.tolist():
     	   #ename = ename + j + "||"
        #addOMIM_disease_en.append(ename.strip('||'))
    ds.insert(ds.columns.tolist().index('HGVSc')+1,'monarch', addgene_disease_9606)
    #ds['subject_label/relation_label/source'] = pd.Series(addgene_disease_9606, index = ds.index)
    ds['IP'] = pd.Series(IP,index = ds.index)
    ds['DI'] = pd.Series(DI,index = ds.index)
    ds['Inheritance'] = pd.Series(Inheritance,index = ds.index)
    ds['Comments'] = pd.Series(Comments,index = ds.index)
    ds['Gene_Locus_MIM_number'] = pd.Series(Gene_Locus_MIM_number,index = ds.index)
    ds['Phenotype'] = pd.Series(Phenotype,index = ds.index)
    ds['Phenotype_MIM_number'] = pd.Series(Phenotype_MIM_number,index = ds.index)
    #####		
    ds['OMIM_id']=pd.Series(addOMIM_id,index=ds.index)
    ds['OMIM_gene']=pd.Series(addOMIM_gene,index=ds.index)
    ds['OMIM_disease_ch']=pd.Series(addOMIM_disease_ch,index=ds.index)
    ds['OMIM_disease_en']=pd.Series(addOMIM_disease_en,index=ds.index)
    ds['chpo_name_CH']=pd.Series(add_chpo_name_CH1,index=ds.index)
    ds['chpo_define_EN']=pd.Series(add_chpo_define_EN1,index=ds.index)
    ds['chpo_define_CH']=pd.Series(add_chpo_define_CH1,index=ds.index)
    # update
    ds['morbidmap_Gene']=pd.Series(add_mor_Gene,index=ds.index)
    ds['morbidmap_Phenotype']=pd.Series(add_mor_Phenotype,index=ds.index)
    ds['MIM_Number']=pd.Series(add_mor_MIM_Number,index=ds.index)
    ds['Cyto_Location']=pd.Series(add_mor_Cyto_Location,index=ds.index)
    #############################
    for i in ds.index:
        taromim = ds.loc[i]['Phenotype_MIM_number']
        taromim_lst = set(taromim.split('|')) # remove repeated elements
        tmp_lst = []
        for each_id in taromim_lst:
     	   if each_id=='':
     		   tmp_lst.append('-')
     	   else:
     		   tmp = OMIM_Phenotype[OMIM_Phenotype['OMIM_id'].str.contains(each_id)]
     		   tmp = tmp['Phenotype'].str.cat(sep = '#').replace('\r','').replace('\n','') # there are line breaks in text columns which cause extra lines
     		   tmp_lst.append(tmp)
        omim_phenotype.append('||'.join(tmp_lst))
    ds['OMIM_Phenotype'] = pd.Series(omim_phenotype,index = ds.index)
    return ds

def filt2(sn, ds):    
    print "\n\nStep2:filting by MAF(public)"
    """ old maf cutoff till 2017/2/28
    conf=ds.esp6500_all.apply(pd.to_numeric,errors='coerce').fillna(0)<0.1
    ds=ds[conf]
    conf=ds.ExAC_ALL.apply(pd.to_numeric,errors='coerce').fillna(0)<0.1
    ds=ds[conf]
    conf=ds.ExAC_EAS.apply(pd.to_numeric,errors='coerce').fillna(0)<0.05
    ds=ds[conf]
    """
    #conf=ds['1000g2015aug_all'].apply(pd.to_numeric,errors='coerce').fillna(0)<0.05
    #ds=ds[conf]
    if '1000g2015aug_all' in ds.columns.tolist() and 'CLINSIG' in ds.columns.tolist() and 'gnomAD_genome_ALL' in ds.columns.tolist() and 'Hom_EAS' in ds.columns.tolist():
        ds = ds[(ds['1000g2015aug_all'].apply(pd.to_numeric,errors='coerce').fillna(0)>=0.05)&(ds['CLINSIG'].str.contains('Pathogenic'))|(ds['1000g2015aug_all'].apply(pd.to_numeric,errors='coerce').fillna(0)<0.05)]
        #ds = ds[(ds['1000g2015aug_all'].apply(pd.to_numeric,errors='coerce').fillna(0)>=0.05)&(ds['CLNSIG'].str.contains('Pathogenic'))|(ds['1000g2015aug_all'].apply(pd.to_numeric,errors='coerce').fillna(0)<0.05)] # annovar 20180416 CLINSIG->CLNSIG
        ds = ds[~((ds['gnomAD_genome_ALL'].apply(pd.to_numeric,errors = 'coerce').fillna(0)>0.5)&(ds['Hom_EAS'].apply(pd.to_numeric,errors = 'coerce').fillna(0)>500))]
    return ds

def filt3(sn, ds):    
    ##step3:filt3_mis/sence/1/1 add gene coverage
    print "\nStep3:filting by mis/sence/1/1 and add gene coverage"
    
    ds_b=ds[~ds['Func.refGene'].isin(['exonic'])]
    ds_a=ds[ds['Func.refGene'].isin(['exonic'])]
    ds_a=ds_a[ds_a['ExonicFunc.refGene']!="synonymous SNV"]
    
    #ds=pd.concat([ds_a2,ds_a1,ds_b])
    ds=pd.concat([ds_a,ds_b])
    if 'Trio_GTinfo' in ds.columns:
    	ds=ds[ds['Depth_'+sn].apply(pd.to_numeric,errors = 'coerce')>7]  # child depth
    	ds=ds[ds['Qual'].apply(pd.to_numeric,errors = 'coerce')>200]     # tiro qual
    	ds=ds[ds['Trio_GTinfo'].apply(lambda x:"".join(set(x.split(":")))) != '1/1'] # filter trio 1/1  why?
    ds.index=range(len(ds))
    if '1000g2015aug_eas' in ds.columns.tolist():
        ds.insert(ds.columns.tolist().index('1000g2015aug_eas')+1,'High_Frequency_Pathogenic',ds['1000g2015aug_all'])
        ds['High_Frequency_Pathogenic'] = ds['High_Frequency_Pathogenic'].apply(pd.to_numeric,errors='coerce').fillna(0).apply(lambda x:'-' if x<0.05 else 'YES')

    ds.insert(ds.columns.tolist().index('Qual')+1,'Gene_20X_cvg',np.nan) # add Gene_average_cvg
    if os.path.exists(sn+'.depth.sample_gene_summary'):
        gene_cvg = pd.read_csv(sn+'.depth.sample_gene_summary',sep = '\t',encoding = 'gbk',header = 0) # read gene_covrage
        gene_cvg = gene_cvg.set_index('Gene') # set index
    
        for i in ds.index:
        	if ds.loc[i]['Gene.refGene'].find(',') != -1:
        		Gene_lst = ds.loc[i]['Gene.refGene'].split(',')
        		average_cvg = 0
        		for gene in Gene_lst:
        			if gene not in gene_cvg.index:continue
        			average_cvg += gene_cvg.loc[gene][sn+'_%_above_20']
        		ds.loc[i,'Gene_20X_cvg'] = average_cvg/float(len(Gene_lst))
        	if ds.loc[i]['Gene.refGene'].find(';') != -1:
        		Gene_lst = ds.loc[i]['Gene.refGene'].split(';')
        		average_cvg = 0
        		for gene in Gene_lst:
        			if gene not in gene_cvg.index:continue
        			average_cvg += gene_cvg.loc[gene][sn+'_%_above_20']
        		ds.loc[i,'Gene_20X_cvg'] = average_cvg/float(len(Gene_lst))
        	if ds.loc[i]['Gene.refGene'] in gene_cvg.index:
        		ds.loc[i,'Gene_20X_cvg'] = gene_cvg.loc[ds.loc[i]['Gene.refGene']][sn+'_%_above_20']

    flag = True                
    for i in ['SiPhy_29way_logOdds','SiPhy_29way_logOdds_rankscore','score_total','esp6500_all','CLINSIG']:
        if i not in ds.columns.tolist():
            flag = False
            break
    if flag:
        ds_tmp3 = ds[['SiPhy_29way_logOdds','SiPhy_29way_logOdds_rankscore','score_total']]
        for i in ds_tmp3.columns:
            ds = ds.drop(i,axis = 1)
        ds_tmp = ds.drop(ds.columns.tolist()[ds.columns.tolist().index('esp6500_all'):ds.columns.tolist().index('CLINSIG')],axis=1)
        #ds_tmp = ds.drop(ds.columns.tolist()[ds.columns.tolist().index('esp6500_all'):ds.columns.tolist().index('CLNALLELEID')],axis=1) # clinvar_20180603, the order and names of clinvar columns changes
        for i in ds_tmp3.columns:
            ds_tmp.insert(ds_tmp.columns.tolist().index('AAChange.refGene')+1,i,ds_tmp3[i]) # adjust the order of columns
        ds_tmp2 = ds.loc[:,'esp6500_all':'PVS1/PS1/PS2/PS3/PS4/PM1/PM2/PM3/PM4/PM5/PM6/PP1/PP2/PP3/PP4/PP5/BA1/BS1/BS2/BS3/BS4/BP1/BP2/BP3/BP4/BP5/BP6/BP7']
        ds = pd.concat([ds_tmp,ds_tmp2],axis = 1)

    if 'Dead_Zone' in ds.columns.tolist():
        ds_tmp = ds[['Inheritance','Comments','Gene_Locus_MIM_number','Phenotype','Phenotype_MIM_number','OMIM_Phenotype']]
        ds = ds.drop(['Inheritance','Comments','Gene_Locus_MIM_number','Phenotype','Phenotype_MIM_number','OMIM_Phenotype'],axis = 1)
        for i in ds_tmp.columns:
            ds.insert(ds.columns.tolist().index('Dead_Zone'),i,ds_tmp[i]) # adjust the order of columns
    if 'HGVSc' in ds.columns.tolist() and 'HGVSp' in ds.columns.tolist():
        ds_tmp = ds[['HGVSc','HGVSp']]
        ds = ds.drop(['HGVSc','HGVSp'],axis = 1)
        for i in ds_tmp.columns:
            ds.insert(ds.columns.tolist().index('AAChange.refGene'),i,ds_tmp[i]) # adjust the order of columns
    if 'link' in ds.columns.tolist():    
        ds_tmp = ds['link']
        ds = ds.drop('link',axis = 1)
        ds.insert(ds.columns.tolist().index('ExAC_ALL'),'gnomAD_link',ds_tmp)  # adjust the order of columns
    if 'CADD_RawScore' in ds.columns.tolist() and 'CADD_PHRED' in ds.columns.tolist():    
        ds_tmp = ds[['CADD_RawScore','CADD_PHRED']]
        ds = ds.drop(['CADD_RawScore','CADD_PHRED'],axis = 1)
        for i in ds_tmp.columns:
            ds.insert(ds.columns.tolist().index('Inheritance'),i,ds_tmp[i]) # adjust the order of columns
    return ds 

#__main__:
for sn in list:
    sn=str(sn) 
    #file_check(sn)
    print "##Sample:%s" % sn   
    ds=pd.read_table(sn+".score",sep="\t",dtype=str)
    if os.path.exists(sn+'.maf'):
        ds = maf_process(sn, ds)
    if os.path.exists(sn+'.link'):
        ds = link_process(sn, ds)
    if os.path.exists(sn+'.CADD'):
        ds = cadd_process(sn, ds)
    ds = del_process(sn, ds)
    ds_pre = prediction_process(sn, ds)
    check_chr(ds_pre)
    
    ##step1:filt1+annotation_AF+annotation_OMIM
    ds_filt1 = filt1(sn, ds_pre, kwargs)
    #annotation:local MAF
    ds_AF = add_AF(sn, ds_filt1, AF)
    #annotation:OMIM_ch/OMIM_en/chpo/morbidmap/Phenotype
    ds_other = add_other(sn, ds_AF)
    ds_other.to_csv(sn+".sheet1.txt",sep="\t",index=False,encoding='utf-8')
    
    ##step2:filting by MAF(public)
    ds_filt2 = filt2(sn, ds_other)
    ds_filt2.to_csv(sn+".sheet2.txt",sep="\t",index=False,encoding='utf-8')
    
    ##step3:filt3_mis/sence/1/1 add gene coverage
    ds_filt3 = filt3(sn, ds_filt2)
    ds_filt3.to_csv(sn+".sheet3.txt",sep="\t",index=False,encoding='utf-8')
    
    ##step4:create xlsx
    print "\nStep4:creating xlsx"
    xlsxgen(sn,ThisVer,c_f_m)
    xlsxgen_children_hospital(sn,ThisVer,c_f_m)
