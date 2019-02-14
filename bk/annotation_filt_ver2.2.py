"""
version:2.2
by:gp
update:2017/2/25
"""


import sys
import pandas as pd
import numpy as np
from lib.GAF import *
from lib.xlsx import *

#reload(sys)
#sys.setdefaultencoding('utf8')

##
print "\nGenetic annotation, ver2.2, update:2017/2/25, by:gp"


##read parameters of cmd line
kwargs_raw = sys.argv[1:]
kwargs={'-edge':'hg19_wxs_edge.bed','-l':'list','-AF':"combine_run1.snp.AF",'-sn':''}
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


#genelist=pd.read_table("genelist_uniq.txt",header=None)[0].tolist()
#ex_bound=pd.read_table("chihos_edge_20.bed",sep="\t")
#ex_exonin=pd.read_table("chihos_exon_in_5.bed",sep="\t")
#read database:
ex_edge=pd.read_table("db/"+kwargs.get('-edge'),sep="\t")
ex_edge=ex_edge.sort_values(['chrom','pos'])

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
#phenon=pd.read_table("db/phenotype_annotation_hpoteam_include_chpo_doid_xref.txt",encoding='gbk',dtype=str)
phenon=phenon[phenon.infosource=='OMIM']
phenon=phenon.fillna("NA")

#database confirm:
print "Database confirming:\n-edge:%s\n-list:%s\n-MAF:%s\nOMIM_ch:%s\nOMIM_en:%s\nHPO:%s\n\n" % ("db/"+kwargs.get('-edge'),kwargs.get('-l'),"db/"+kwargs.get('-AF'),"db/OMIM_selected_ch.csv","db/OMIM_selected_en.csv","db/phenotype_annotation_hpoteam_include_chpo_doid_xref.csv")




for sn in list:
    sn=str(sn) 
    print "##Sample:%s" % sn   
    ds=pd.read_table(sn+".score",sep="\t",dtype=str)
    ds['Start']=ds['Start'].astype(long)
    #ds=ds.sort_values(['Chr','Start']) #dealed by GATK
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
    
    tar=['exonic','splicing','exonic;splicing']
    ds_a=ds[ds['Func.refGene'].isin(tar)].copy() 
    ds_b=ds[~ds['Func.refGene'].isin(tar)].copy()
    ds_b.index=range(len(ds_b))
    ds_b['GeneDetail.refGene']='NA'
    for i in ds_b.index:
        line=ds_b.loc[i]
        chrom=line['Chr']
        pos=line['Start']
        
        #ex_sub=ex_edge[ex_edge['chrom']==chrom]
        ex_sub=ex_edge.loc[ex_init.get(chrom):ex_end.get(chrom)]
        
        """slow:
        ex_pos=ex_sub.pos.tolist()
        ex_index = range(len(ex_pos))
        posmin=ex_index[0]
        posmax=ex_index[-1]
        while(1):
            if(posmax-posmin<=2):
                break
            posmid=(posmin+posmax)/2
            if(pos>ex_pos[posmid]):
                posmin=posmid
            else:
                posmax=posmid
        
        ex_sub=ex_sub.loc[posmin:posmax]
        
        """
        ex_sub=ex_sub[ex_sub['pos']>(pos-16)]
        ex_sub=ex_sub[ex_sub['pos']<(pos+100)]
        ex_mindis=16
        
        for j in ex_sub.index:
            ex_pos=ex_sub.loc[j].pos
            if(abs(pos-ex_pos)<=15):
                ex_mindis=abs(pos-ex_pos)
        if(ex_mindis<16):
            ds_b.loc[i,['GeneDetail.refGene']] = "exon_edge_"+str(ex_mindis)
        
        if(i%1000==0 or i==ds_b.index[-1]):
           view_bar(i+1,len(ds_b))
            
        """ slow & bugs:
        ex_sub_left=ex_sub[ex_sub.side=="left"]
        ex_sub_right=ex_sub[ex_sub.side=="right"]
        
        Bugs: gene on different strands should be dealed differently when confirm distance to exons; here use both sides 
        
        for j in range(15):
            j=j+1
            pos_re=pos-j
            if(pos_re in ex_sub_right.pos.tolist() or pos_re in ex_sub_left.pos.tolist()):
                ds_b.loc[i]['GeneDetail.refGene']="exon_edge_"+str(j)
            pos_re=pos+j
            if(pos_re in ex_sub_left.pos.tolist() or pos_re in ex_sub_right.pos.tolist()):
                ds_b.loc[i]['GeneDetail.refGene']="exon_edge_"+str(j)
            
        
        for j in range(5):
            j=j+1
            pos_re=pos-j
            if(pos_re in ex_sub_right.pos.tolist()):
                ds_b.loc[i]['GeneDetail.refGene']="exonin_+"+str(j)
            pos_re=pos+j
            if(pos_re in ex_sub_left.pos.tolist()):
                ds_b.loc[i]['GeneDetail.refGene']="exonin_-"+str(j)
        """
    
    ds_c=ds_b[ds_b['GeneDetail.refGene']!='NA']
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
    conf=ds.1000g2014oct_all.apply(pd.to_numeric,errors='coerce').fillna(0)<0.05
    ds=ds[conf]
    ds.to_csv(sn+".sheet2.txt",sep="\t",index=False,encoding='gbk')
    
    
    ##step3:filt3_mis/sence
    print "\nStep3:filting by mis/sence"
    """ old
    tags_exonin=ex_exonin['tags'].tolist()
    
    for i in ds.index:
        dsline=ds.loc[i]
        dstag=dsline['Chr']+":"+str(dsline['Start'])
        addinfo="NA"
        if(dstag in tags_exonin):
            addinfo=ex_exonin[ex_exonin['tags']==dstag].iloc[0].pos
        elif(dstag in tags_bound):
            addinfo=ex_bound[ex_bound['tags']==dstag].iloc[0].pos
        ds.loc[i,'GeneDetail.refGene']=addinfo

    ds_a=ds[ds['Func.refGene'].isin(['exonic'])]
    ds_a1=ds_a[ds_a['GeneDetail.refGene']!='NA']
    ds_a2=ds_a[ds_a['GeneDetail.refGene']=='NA']
    ds_a2=ds_a2[ds_a2['ExonicFunc.refGene']!="synonymous SNV"]
    """
    ds_b=ds[~ds['Func.refGene'].isin(['exonic'])]
    ds_a=ds[ds['Func.refGene'].isin(['exonic'])]
    ds_a=ds_a[ds_a['ExonicFunc.refGene']!="synonymous SNV"]
    
    #ds=pd.concat([ds_a2,ds_a1,ds_b])
    ds=pd.concat([ds_a,ds_b])
    ds.index=range(len(ds))
    ds.to_csv(sn+".sheet3.txt",sep="\t",index=False,encoding='gbk')
    
    
    ##step4:create xlsx
    print "\nStep4:creating xlsx"
    xlsxgen(sn)
