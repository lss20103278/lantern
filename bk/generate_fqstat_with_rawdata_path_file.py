#!/usr/bin/env python
# -*- coding: utf-8 -*-
#########################################################################
# File Name: sub_fqstat.py
# Author: sslyu
# mail: 1573077420@qq.com
# Created Time:Thu 11 Oct 2018 03:36:20 PM CST
#########################################################################

import numpy as np
import pandas as pd
import pandas as pd

sn = pd.read_table('list_46_highrisk_infant_3rd_group_reads_miss_ID', header=None)[0].tolist()

for i in sn:
    ofile = open(i+'_fqstat.sh', 'w')
    rawdata_path = open(i+'_rawdata_path').readlines()
    R1 = rawdata_path[0].strip('\n')
    R2 = rawdata_path[1].strip('\n')
    ofile.write('#!/bin/bash\n\n#SBATCH -J '+i+'\n#SBATCH -p BP10\n#SBATCH -N 1\n#SBATCH -n 4\n\nsample='+i+'\n\nR1='+R1+'\nR2='+R2+'\n\n/home/ana005/anaconda2/bin/iTools Fqtools stat -InFq $R1 -InFq $R2 -OutStat $sample.info\n')
    ofile.close()
