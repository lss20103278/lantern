# Usage:?
import pandas as pd

samples=pd.read_table("list",header=None)[0].tolist()

for sn in samples:
    data=pd.read_table(sn+".ann.hg19_multianno.txt",sep="\t")
    dinfo=pd.read_csv(sn+".csv")
    data['SEGREGATION_SCORE']=dinfo['SEGREGATION_SCORE']
    data['QUAL']=dinfo['QUAL']
    data['DENOVO']=dinfo['DENOVO']
    data['depth']=dinfo['depth']
    data['gt_father']=dinfo['gt_father']
    data['gt_mother']=dinfo['gt_mother']
    data['gt_child']=dinfo['gt_child']
    data['gt_seq_father']=dinfo['gt_seq_father']
    data['gt_seq_mother']=dinfo['gt_seq_mother']
    data['gt_seq_child']=dinfo['gt_seq_child']
    data['heterozygosity']=dinfo['heterozygosity']
    data.to_csv(sn+".ann.txt",index=False,sep="\t")
