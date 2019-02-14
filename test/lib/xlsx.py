import pandas as pd
import sys
import os
import copy
reload(sys)
sys.setdefaultencoding('utf8')

def GT_Flag(lst):
    result= ['0','0','0']
    if lst[0] ==lst[1]:
        result[0]='1'
    if lst[0] == lst[2]:
        result[1]='1'
    if lst[1] == lst[2]:
        result[2] ='1'
    return ':'.join(result)

def add_tmp(ds):
    GTinfo = ds['Trio_GTinfo'].str.split(':')
    GT_flag = []
    for i in GTinfo.values:
        GT_flag.append(GT_Flag(i))
		# GT_flag is a list containing three elements, the first is whether the genotype of child is the same as that of father, the second is whether the genotype of child is the same as that of mother, the third is whether the genotype of father is the same as that of mother, 0 stands for different, 1 stands for the same
    ds['tmp'] = pd.Series(GT_flag, index=ds.index)
    return ds

def filt3_process(sn,c_f_m):
    filt3=pd.read_table(sn.strip()+".sheet3.txt",sep="\t",dtype=str,encoding='utf-8')
    a = []
    for i in filt3.columns:
    	if 'GTinfo_' in i:
    		a.append(i[7:])
    if len(a) == 0:
        a.append(sn) # sn should be the same with the name in sn.link
    if c_f_m == 'c':
        list_name = [sn,'father','mother']
        col_name = ['GTinfo_', 'GT_', 'Qual_']
        for i in range(len(list_name)):
            filt3.insert(filt3.columns.tolist().index('HI_Predictions')+i*3+1,col_name[i]+list_name[i],'-')
        filt3.rename(columns={'depth':'Depth_'+a[0]}, inplace=True)
        list_name2 = ['father', 'mother']
        for i in list_name2:
            filt3.insert(filt3.columns.tolist().index('Depth_'+a[0])+1+list_name2.index(i),'Depth_'+i,'-')
        list_name3 = ['Trio_GTinfo','SEGSCORE','DENOVO','DENOVO_Heterozygosity','lowGQ']
        for i in list_name3:
            filt3.insert(filt3.columns.tolist().index('Qual')+list_name3.index(i),i,'-')
    elif c_f_m == 'c_f':
        for i in ['GTinfo_','GT_','Qual_','Depth']:
            filt3.insert(filt3.columns.tolist().index(i+a[0])+2,i+'_mother','-')
    elif c_f_m == 'c_m':
        for i in ['GTinfo_','GT_','Qual_','Depth']:
            filt3.insert(filt3.columns.tolist().index(i+a[0])+2,i+'_father','-')
    return filt3            

def filt3_process2(sn,filt3):    
    filt3.to_csv(sn+'.sheet3.txt.tmp', index=None, sep="\t", encoding='utf-8')
    os.system('/DATA/ypliu/opt/R-3.4.3/bin/Rscript /DATA/sslyu/trio_BWA-GATK_3.0/src/Data_transform_v3.R '+sn.strip()+'.sheet3.txt.tmp > '+sn+'.sheet3.tsv')
    filt3=pd.read_table(sn.strip()+".sheet3.tsv",sep="\t",dtype=str,encoding='utf-8')
    filt3.insert(filt3.columns.tolist().index('AAChange.refGene')+1,'SIFT/Polyphen2_HVAR/MutationTaster/FATHMM/M-CAP/REVEL',filt3['SIFT_pred']+"/"+filt3['Polyphen2_HVAR_pred']+"/"+filt3['MutationTaster_pred']+"/"+filt3['FATHMM_pred']+"/"+filt3['M-CAP_pred']+"/"+filt3['REVEL'].replace('NA',0).apply(lambda x:'D' if float(x)>0.5 else '.')) # the source of .
    filt3.score_total=filt3.score_total.apply(pd.to_numeric,errors='coerce')
    filt3=filt3.sort_values(['score_total'],ascending=0)
    return filt3

def filt4_process(filt3):    
    filt4 = filt3[(filt3['Func.refGene'].apply(lambda x : 'splicing' in x))|(filt3['Heterozygosity']=='hom')| \
        (filt3['ExonicFunc.refGene'].fillna('nothing').apply(lambda x: 'stoploss' in x or 'stopgain' in x)) | \
        (filt3['ExonicFunc.refGene']=='frameshift deletion')| \
        (filt3['ExonicFunc.refGene']=='frameshift insertion')| \
        (filt3['SIFT_pred']=='D') | (filt3['Polyphen2_HVAR_pred'].apply(lambda x: x=="D" or x=="P")) | \
        (filt3['min_edge'].notnull())]
    filt4 = filt4[filt4['Func.refGene'] != 'ncRNA_intronic']
    return filt4

def filt_het_1_process(filt_het, gene_count_dict):        
    filt_het_1 = filt_het.copy()
    gene_count_lst_1 = [k for k, v in gene_count_dict.items() if v == 1]
    for i in filt_het_1.index:
        if filt_het_1.loc[i,'Gene.refGene'] not in gene_count_lst_1:
            filt_het_1 = filt_het_1.drop(i)
    filt_het_1 = filt_het_1[
        (filt_het_1['InheritanceStatus'].str.contains('AD')) | (filt_het_1['InheritanceStatus'].str.contains('CX'))] 
    filt_het_1.index = range(len(filt_het_1))
    filt_het_1 = add_tmp(filt_het_1)
    filt_het_1 = filt_het_1[filt_het_1['SIFT/Polyphen2_HVAR/MutationTaster/FATHMM/M-CAP/REVEL'].fillna('-').str.startswith(('D/D', 'D/P', 'D/B', 'T/D', './D', 'D/.'))] # only one loci occurs in this Gene.refGene and its InheritanceStatus is AD or CX and its SIFT/Polyphen2_HVAR/MutationTaster/FATHMM/M-CAP/REVEL is D/D or D/P or D/B or T/D or ./D or D/.
    return filt_het_1

def filt_het_2_process(filt_het, gene_count_dict):
    filt_het_2 = filt_het.copy()
    gene_count_lst_2_3 = [k for k, v in gene_count_dict.items() if 1 < v <= 3]
    filt_het_2.index = range(len(filt_het_2))
    for i in filt_het_2.index:
        if filt_het_2.loc[i,'Gene.refGene'] not in gene_count_lst_2_3:
            filt_het_2 = filt_het_2.drop(i)
    filt_het_2 = add_tmp(filt_het_2)
    gene_tmp = []
    for each_gene in gene_count_lst_2_3:
        filt_het_2_tmp = filt_het_2[filt_het_2['Gene.refGene'] == each_gene]
        gene_tmp_flag = filt_het_2_tmp['tmp'].tolist()
        if '0:1:0' in gene_tmp_flag and '1:0:0' in gene_tmp_flag: # 0:1:0 the genotype of child is different from that of the father but the same with that of the mother; 1:0:0 the genotype of the child is the same with that of the father but different from that of the mother
            if filt_het_2_tmp['SIFT/Polyphen2_HVAR/MutationTaster/FATHMM/M-CAP/REVEL'].fillna('-').str.startswith(('D/D', 'D/P', 'D/B', 'T/D', './D', 'D/.')).tolist().count(True):
                gene_tmp.append(each_gene)
    for i in filt_het_2.index:
        if filt_het_2.loc[i,'Gene.refGene'] not in gene_tmp:
            filt_het_2 = filt_het_2.drop(i)
    return filt_het_2

def filt5_process(filt3):        
    filt_tmp = filt3[filt3['Heterozygosity'] != 'unexact']

    filt_hom = filt_tmp[filt_tmp['Heterozygosity'] == 'hom']
    filt_hom.index = range(len(filt_hom))
    filt_hom = add_tmp(filt_hom)
    filt_hom = filt_hom[(filt_hom['tmp'] == '0:0:1') | (filt_hom['tmp'] == '0:0:0') | (filt_hom['DENOVO'] != '-')]
    filt_hom = filt_hom[filt_hom['SIFT/Polyphen2_HVAR/MutationTaster/FATHMM/M-CAP/REVEL'].fillna('-').str.startswith(('D/D', 'D/P', 'D/B', 'T/D', './D', 'D/.'))]

	#################### het
    filt_het = filt_tmp[filt_tmp['Heterozygosity'] == 'het']
    gene_count_dict = dict(filt_het['Gene.refGene'].value_counts())
    gene_count_lst = [k for k, v in gene_count_dict.items() if v > 3]
    for each_gene in gene_count_lst:
        filt_het = filt_het[filt_het['Gene.refGene'] != each_gene]

    filt_het_1 = filt_het_1_process(filt_het, gene_count_dict)
    filt_het_2 = filt_het_2_process(filt_het, gene_count_dict)
    filt_het = pd.concat([filt_het_1, filt_het_2])

    filt5 = pd.concat([filt_hom, filt_het])
    del filt5['tmp']
    filt5 = filt5.sort_values(['score_total'], ascending=0)
    return filt5

def xlsxgen(sn,ver,c_f_m):
    print "creating xlsx for sample : %s" %(sn)
    wb=pd.ExcelWriter(sn.strip()+"_"+ver+".xlsx",engine='openpyxl')
    readme = pd.read_excel('/DATA/sslyu/trio_BWA-GATK_3.0/doc/readme.xlsx',sheetname = 'readme',header = 0)
    filt1=pd.read_table(sn.strip()+".sheet1.txt",sep="\t",dtype=str,encoding='utf-8')
    filt2=pd.read_table(sn.strip()+".sheet2.txt",sep="\t",dtype=str,encoding='utf-8')
    filt3 = filt3_process(sn,c_f_m)
    filt3 = filt3_process2(sn,filt3)
    filt4 = filt4_process(filt3)

    filt1.to_excel(wb,"filt1",na_rep="-",index=False)
    filt2.to_excel(wb,"filt2",na_rep="-",index=False)
    filt3.to_excel(wb,"selected_mutations",na_rep="-",index=False)
    filt4.to_excel(wb,"splicing,stop gain loss,DD,DP",na_rep="-",index=False)
    if c_f_m == 'c_f_m':
        filt5 = filt5_process(filt3)
        filt5.to_excel(wb,'filt_selected_mutations',na_rep = "-",index = False)
    readme.to_excel(wb,'readme',index = False)

    ws = wb.sheets['selected_mutations']
    for i in range(2,len(filt3)+1):
        i = str(i)
        ws['AA'+i].hyperlink = "http://gnomad.broadinstitute.org/variant/"+ws['AA'+i].value    
    wb.save()

def xlsxgen_children_hospital(sn,ver,c_f_m):
    print "creating xlsx for sample : %s" %(sn)
    wb=pd.ExcelWriter(sn.strip()+"_children_hospital_"+ver+".xlsx",engine='openpyxl')
    readme = pd.read_excel('/DATA/sslyu/trio_BWA-GATK_3.0/doc/readme.xlsx',sheetname = 'readme',header = 0)
    filt1=pd.read_table(sn.strip()+".sheet1.txt",sep="\t",dtype=str,encoding='utf-8')
    filt2=pd.read_table(sn.strip()+".sheet2.txt",sep="\t",dtype=str,encoding='utf-8')

    filt3=pd.read_table(sn.strip()+".sheet3.txt",sep="\t",dtype=str,encoding='utf-8')
    filt3.insert(filt3.columns.tolist().index('AAChange.refGene')+1,'SIFT/Polyphen2_HVAR/MutationTaster/FATHMM/M-CAP/REVEL',filt3['SIFT_pred']+"/"+filt3['Polyphen2_HVAR_pred']+"/"+filt3['MutationTaster_pred']+"/"+filt3['FATHMM_pred']+"/"+filt3['M-CAP_pred']+"/"+filt3['REVEL'].replace('NA',0).apply(lambda x:'D' if float(x)>0.5 else '.')) # the source of .
    filt3.score_total=filt3.score_total.apply(pd.to_numeric,errors='coerce')
    filt3=filt3.sort_values(['score_total'],ascending=0)
    filt4 = filt4_process(filt3)


    filt1.to_excel(wb,"filt1",na_rep="-",index=False) # the source of -
    filt2.to_excel(wb,"filt2",na_rep="-",index=False) # the source of -
    filt3.to_excel(wb,"selected_mutations",na_rep="-",index=False) # the source of -
    filt4.to_excel(wb,"splicing,stop gain loss,DD,DP",na_rep="-",index=False) # the source of -

    if c_f_m == 'c_f_m':
        filt5 = filt5_process(filt3)
        filt5.to_excel(wb,'filt_selected_mutations',na_rep = "-",index = False)
    readme.to_excel(wb,'readme',index = False)
    ws = wb.sheets['selected_mutations']
    for i in range(2,len(filt3)+1):
        i = str(i)
        ws['AA'+i].hyperlink = "http://gnomad.broadinstitute.org/variant/"+ws['AA'+i].value    
    wb.save()

