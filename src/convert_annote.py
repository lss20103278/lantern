import sys
import pandas as pd

sn=sys.argv[1]


data=pd.read_table(sn+".txt",sep="\t")

"""
dadd=pd.read_table(sn+".short.txt",sep="\t")

data['gt_seq_father']=dadd.gt_seq_father
data['gt_seq_mother']=dadd.gt_seq_mother
data['gt_seq_child']=dadd.gt_seq_child


del data['CONDEL']
del data['ENSEMBL_GENE']
del data['OVERALL_SCORE']
del data['ANNOTATION_SCORE']
"""



#REF/ALT types
def convert_annovar_type(ref,alt):
    if(len(ref)==len(alt)):
        return([ref,alt])
    elif(len(ref)>len(alt)):
        cutlen=len(alt)
        ref=ref[cutlen:]
        alt='-'
        return([ref,alt])
    else:
        cutlen=len(ref)
        ref='-'
        alt=alt[cutlen:]
        return([ref,alt])

#Start/End types
def convert_annovar_pos(ref,alt,pos):
    if(len(ref)==len(alt)):
        return([pos,pos])
    elif(len(ref)>len(alt)):
        return([pos+len(alt),pos + len(ref) - len(alt)])
    else:
        return([pos+len(ref),pos + len(ref)])


#extract depth
def extract_depth(x):
    if(x=="./."):
        return(".")
    else:
        return((x.split(":")[2]))
    
    
#filt lowGQ
def lowGQ_child(x):
    if("lowGQ" in x):
        return(False)
    else:
        return(True)

#extract_gt
def extract_gt(ref,alt,x):
    if(x=="./."):
        return("./.")
    else:
        gt=x.split(":")[0]
        if(gt=="0/0"):
            return(ref+"/"+ref)
        if(gt=="0/1"):
            return(ref+"/"+alt)
        if(gt=="1/1"):
            return(alt+"/"+alt)
        return("./.")    


#het
def extract_het(x):
    gt=x.split(":")[0]
    if(gt=="0/0"):
        return('hom')
    if(gt=="0/1"):
        return('het')
    if(gt=="1/1"):
        return('hom')
    return("unexact")




conf=[]
for i in data.index:
    conf.append(lowGQ_child(data.loc[i].gt_child))

data=data[conf]


gt_seq_child=[]
gt_seq_father=[]
gt_seq_mother=[]
het=[]
for i in data.index:
    gt_seq_child.append(extract_gt(data.loc[i].REF,data.loc[i].ALT,data.loc[i].gt_child))
    gt_seq_father.append(extract_gt(data.loc[i].REF,data.loc[i].ALT,data.loc[i].gt_father))
    gt_seq_mother.append(extract_gt(data.loc[i].REF,data.loc[i].ALT,data.loc[i].gt_mother))
    het.append(extract_het(data.loc[i].gt_child))

data['gt_seq_child']=pd.Series(gt_seq_child,index=data.index)
data['gt_seq_father']=pd.Series(gt_seq_father,index=data.index)
data['gt_seq_mother']=pd.Series(gt_seq_mother,index=data.index)
data['heterozygosity']=pd.Series(het,index=data.index)



#data.rename(columns={'X1000g2014oct_all':'MAF_1000G_all'})
conf2=(data.gt_seq_child!=data.REF+"/"+data.REF) & (data.gt_seq_child!='./.')
data=data[conf2]



depth=[]
for i in data.index:
    gt=data.loc[i].gt_child
    depth.append(extract_depth(gt))


data['depth']=pd.Series(depth,index=data.index)
data.to_csv(sn+".csv",index=False)








reann=pd.DataFrame()
reann['chr']=data.CHROM
reann['pos']=data.POSITION
reann['ref']=data.REF
reann['alt']=data.ALT


ofile=open(sn+".avinput",'w')
for i in reann.index:
    newformat=convert_annovar_type(reann.loc[i].ref,reann.loc[i].alt)
    newpos=convert_annovar_pos(reann.loc[i].ref,reann.loc[i].alt,reann.loc[i].pos)
    ofile.write("%s\t%d\t%d\t%s\t%s\n" % (reann.loc[i].chr,newpos[0],newpos[1],newformat[0],newformat[1]))
ofile.close()








