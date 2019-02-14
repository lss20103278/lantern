"""
version:2.4
by:gp
update:2017/3/28
"""
ThisVer='v2.4'
Thisreleasetime='2017/3/28'

import sys
import pandas as pd
import numpy as np
from lib.GAF import *
from lib.xlsx import *



##verion confirm
print "\nGenetic annotation, %s, update:%s, by:gp" % (ThisVer,Thisreleasetime)
##


##read parameters of cmd line
kwargs_raw = sys.argv[1:]
kwargs={'-edge':'ucsc_hg19_refGene.csv','-l':'list','-AF':"combine_run1.snp.AF",'-sn':''}
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


#read edge unkown genes:
ex_edge=pd.read_csv("db/"+kwargs.get('-edge'))
ex_edge=ex_edge.sort_values(['chrom'])

#build edge dic:  used to accelerate bellow
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

AF=pd.read_table("db/"+kwargs.get('-AF'),sep="\t")

OMIM_selected_ch=pd.read_csv("db/OMIM_selected_ch.csv",encoding='gbk',dtype=str)
OMIM_selected_en=pd.read_csv("db/OMIM_selected_en.csv",encoding='gbk',dtype=str)
phenon=pd.read_csv("db/phenotype_annotation_hpoteam_include_chpo_doid_xref.csv",encoding='gbk',dtype=str)
phenon=phenon[phenon.infosource=='OMIM']
phenon=phenon.fillna("NA")

#database confirm:
print """
Database confirming:
-edge:%s
-list:%s
-locMAF:%s
OMIM_ch:%s
OMIM_en:%s
HPO:%s
""" % ("db/"+kwargs.get('-edge'),kwargs.get('-l'),"db/"+kwargs.get('-AF'),"db/OMIM_selected_ch.csv","db/OMIM_selected_en.csv","db/phenotype_annotation_hpoteam_include_chpo_doid_xref.csv")



#__main__:
for sn in list:
    sn=str(sn) 
    print "##Sample:%s" % sn   
    ds=pd.read_table(sn+".score",sep="\t",dtype=str)
    ds['Start']=ds['Start'].astype(long)
    
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
    
    ##step1:filt1+annotation_AF+annotation_OMIM
    print "\nStep1:\n--filting by snp/indel location"
    
    ds.insert(7,'min_edge',[np.nan]*len(ds))
    tar=['exonic','splicing','exonic;splicing']
    ds_a=ds[ds['Func.refGene'].isin(tar)].copy() 
    ds_b=ds[~ds['Func.refGene'].isin(tar)].copy()
    ds_b.index=range(len(ds_b))
    for i in ds_b.index:
        line=ds_b.loc[i]
        chrom=line['Chr']
        pos=int(line['Start'])
        genes=line['Gene.refGene'].replace(",",";").split(";")
        
        #sub of edge by chr
        """ bugs maybe: the end should +1? """
        ex_sub=ex_edge.loc[ex_init.get(chrom):ex_end.get(chrom)]        
        ex_sub=ex_edge[ex_edge['symbol'].isin(genes)]
        refedge='::edge::'
        thisminedge=11
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
                    if(eStart[k] > pos and eStart[k]-pos <= 10):
                        if(eStart[k]-pos < thisminedge): thisminedge=eStart[k]-pos
                        confheader='-'
                        refedge=refedge + ex_subline['name']+":"+ex_subline.symbol+":exon"+str(k+1)+":-"+str(eStart[k]-pos)+";"
                    if(eEnd[k] < pos and pos-eEnd[k] <= 10):
                        if(pos-eEnd[k] < thisminedge): thisminedge=pos-eEnd[k]
                        confheader='+'
                        refedge=refedge + ex_subline['name']+":"+ex_subline.symbol+":exon"+str(k+1)+":+"+str(pos-eEnd[k])+";"                
            elif(ex_subline.strand=='-'):                                
                for k in range(eCount):
                    if(eEnd[k] < pos and pos-eEnd[k] <= 10):
                        if(pos-eEnd[k] < thisminedge): thisminedge=pos-eEnd[k]
                        confheader='-'
                        refedge=refedge + ex_subline['name'] + ":" + ex_subline.symbol + ":exon" + str(eCount-k) + ":-" + str(pos-eEnd[k]) + ";"
                    if(eStart[k]>pos and eStart[k]-pos <= 10):
                        if(eStart[k]-pos < thisminedge): thisminedge=eStart[k]-pos
                        confheader='+'
                        refedge=refedge + ex_subline['name'] + ":" + ex_subline.symbol + ":exon" + str(eCount-k) + ":+" + str(eStart[k]-pos) + ";"
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
    ds_c=ds_b[conf]
    ds=pd.concat([ds_a,ds_c])
    ds.index=range(len(ds))
    
    
    
    #annotation:local MAF
    print "\n--adding local AF"
    af=[]
    for i in ds.index:
        line=ds.loc[i]
        tag=ds.loc[i].Chr+":"+str(ds.loc[i].Start)+":"+ds.loc[i].Ref+">"+ds.loc[i].Alt
        flt=AF[AF['tag']==tag]
        if(len(flt)==1):
            af.append(flt.iloc[0].AF)
        else:
            af.append('.')
    ds['locAF']=pd.Series(af,index=ds.index)		#end
    
    #annotation:OMIM_ch/OMIM_en/chpo
    print "--adding OMIM_ch/OMIM_en/CHPO"
    addOMIM_id=[]
    addOMIM_gene=[]
    addOMIM_disease_ch=[]
    addOMIM_disease_en=[]
    add_chpo_name_CH1=[]
    add_chpo_define_EN1=[]
    add_chpo_define_CH1=[]
    for i in ds.index:
        if(i%1000==0 or i==ds.index[-1]):
            view_bar(i+1,len(ds.index))
        
        targene=ds.loc[i]['Gene.refGene'].split(',')[0]
        info_ch=OMIM_selected_ch[OMIM_selected_ch.OMIM_gene==targene]
        info_en=OMIM_selected_en[OMIM_selected_en.OMIM_gene==targene]
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

        ename=''
        for j in info_en.OMIM_disease.tolist():
            ename = ename + j + "||"
        addOMIM_disease_en.append(ename.strip('||'))
    ds['OMIM_id']=pd.Series(addOMIM_id,index=ds.index)
    ds['OMIM_gene']=pd.Series(addOMIM_gene,index=ds.index)
    ds['OMIM_disease_ch']=pd.Series(addOMIM_disease_ch,index=ds.index)
    ds['OMIM_disease_en']=pd.Series(addOMIM_disease_en,index=ds.index)
    ds['chpo_name_CH']=pd.Series(add_chpo_name_CH1,index=ds.index)
    ds['chpo_define_EN']=pd.Series(add_chpo_define_EN1,index=ds.index)
    ds['chpo_define_CH']=pd.Series(add_chpo_define_CH1,index=ds.index)
    
    ds.to_csv(sn+".sheet1.txt",sep="\t",index=False,encoding='gbk')
    
    
    ##step2:filt_MAF
    print "\n\nStep2:filting by MAF(public)"
    """ old maf cutoff till 2017/2/28
    conf=ds.esp6500_all.apply(pd.to_numeric,errors='coerce').fillna(0)<0.1
    ds=ds[conf]
    conf=ds.ExAC_ALL.apply(pd.to_numeric,errors='coerce').fillna(0)<0.1
    ds=ds[conf]
    conf=ds.ExAC_EAS.apply(pd.to_numeric,errors='coerce').fillna(0)<0.05
    ds=ds[conf]
    """
    conf=ds['1000g2014oct_all'].apply(pd.to_numeric,errors='coerce').fillna(0)<0.05
    ds=ds[conf]
    ds.to_csv(sn+".sheet2.txt",sep="\t",index=False,encoding='gbk')
    
    
    ##step3:filt3_mis/sence
    print "\nStep3:filting by mis/sence"
    ds_b=ds[~ds['Func.refGene'].isin(['exonic'])]
    ds_a=ds[ds['Func.refGene'].isin(['exonic'])]
    ds_a=ds_a[ds_a['ExonicFunc.refGene']!="synonymous SNV"]
    
    #ds=pd.concat([ds_a2,ds_a1,ds_b])
    ds=pd.concat([ds_a,ds_b])
    ds.index=range(len(ds))
    ds.to_csv(sn+".sheet3.txt",sep="\t",index=False,encoding='gbk')
    
    
    ##step4:create xlsx
    print "\nStep4:creating xlsx"
    xlsxgen(sn,ThisVer)
