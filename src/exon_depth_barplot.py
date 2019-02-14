#!/usr/bin/env python
# -*- coding: utf-8 -*-
#########################################################################
# File Name: exon_depth_barplot.py
# Author: sslyu
# mail: 1573077420@qq.com
# Created Time:Fri 31 Aug 2018 01:22:16 PM CST
#########################################################################

import numpy as np
import pandas as pd
import sys
reload(sys)
sys.setdefaultencoding("utf8")
import os
import re
#import matplotlib.pyplot as plt

kwargs_raw = sys.argv[1:]
kwargs = {'-gene':'', '-sn':'', '-l':'list'}
for i in range(len(kwargs_raw)):
    for j in kwargs.keys():
        if kwargs_raw[i] == j:
            kwargs.update({j:kwargs_raw[i+1]})

if len(kwargs_raw) < 1:
    print """
    Examples:
    python /DATA/sslyu/trio_BWA-GATK_3.0/src/exon_depth_barplot.py -gene genename -sn sampleID
    python /DATA/sslyu/trio_BWA-GATK_3.0/src/exon_depth_barplot.py -gene genename -l filename of files containing list of sampleIDs
    """
    sys.exit()

#read sample list
if(kwargs.get('-sn')!=''):
    list = []
    list.append(kwargs.get('-sn'))
    print "\nNote:single sample mode"
else:
    if kwargs.get('-l')!='':
        list = pd.read_table(kwargs.get('-l'), header=None)[0].tolist()
        print "\nNote:list samples mode"
            
# read gene name
gene = kwargs.get('-gene')

##########step1: prepare gene bed file ####################
if not os.path.exists(gene+'.bed'):
    print """
    grep -w gene /DATA/sslyu/soft/annovar/humandb/refseq.sorted.txt |cut -f 3,10-11 |awk '{l=split($2,a,","); split($3,b,","); for (i=1;i<l;i++){OFS="\t"; print $1,a[i],b[i]}}' |sort |uniq > gene.bed
    """
    sys.exit()
    os.system('grep -w '+gene+' /DATA/sslyu/soft/annovar/humandb/refseq.sorted.txt |cut -f 3,10-11 |awk \'{l=split($2,a,","); split($3,b,","); for (i=1;i<l;i++){OFS="\\t"; print $1,a[i],b[i]}}\' |sort |uniq > '+gene+'.bed')
    #os.system('grep -w '+gene+' /DATA/sslyu/soft/annovar/humandb/refseq.sorted.txt |cut -f 3,10-11,13')
    #chrom = raw_input('please enter the chromosome(e.g. chr1):')
    #start = raw_input('please enter the start coordinates of the exons of '+gene+':')
    #end = raw_input('please enter the end coordinates of the exons of '+gene+':')
    #bed_file = open(gene+'.bed', 'w')
    #start = start.split(',')
    #end = end.split(',')
    #for i in range(len(start)):
    #    bed_file.write(chrom+'\t'+start[i]+'\t'+end[i]+'\n')
    #bed_file.close()
###########step2: calculate depth and coverage #################
for sn in list:    
    if not os.path.exists(sn+'_'+gene+'.bed.depth'):
        print """
        samtools depth -b gene.bed sn.sort.mkdup.bam > sn_gene.bed.depth
        """
        sys.exit()
        #os.system('samtools depth -b '+gene+'.bed bam/'+sn+'.sort.mkdup.bam > '+sn+'_'+gene+'.bed.depth')
    exon_depth = pd.DataFrame(columns=['count'])
    exon_depth_v = []
    ix = []
    with open(gene+'.bed') as f:
        for l in f:
            l = l.strip('\n').split('\t')
            ix.append(l[1]+'-'+l[2])
            with open(sn+'_'+gene+'.bed.depth') as f1:
                depth=0
                for l1 in f1:
                    l1 = l1.strip('\n').split('\t')
                    if int(l1[1])>=int(l[1]) and int(l1[1])<=int(l[2]):
                        depth = depth+int(l1[2])
            cdslen = int(l[2])-int(l[1])+1
            exon_depth_v.append(float(depth)/float(cdslen))
    exon_depth['count'] = exon_depth_v
    exon_depth.index = ix
    exon_depth.to_csv(sn+'_'+gene+'_exon.depth.tsv', header=None, sep="\t")

    ##########step3: plotting #######################################
    os.system('Rscript /DATA/sslyu/trio_BWA-GATK_3.0/src/exon_depth_barplot.R '+sn+'_'+gene+'_exon.depth.tsv '+gene+' '+sn)
    #plt.figure()
    #n = exon_depth.shape[0]
    #ix = [i+1 for i in range(n)]
    #exon_depth.index = ix
    #exon_depth.plot(kind="bar")
    #plt.title(gene+' exon region count of '+sn)
    #plt.savefig(sn+'_'+gene+'_exon_count.png')
