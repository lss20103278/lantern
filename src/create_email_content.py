#!/usr/bin/env python
# -*- coding: utf-8 -*-
#########################################################################
# File Name: /DATA/sslyu/trio_BWA-GATK_3.0/src/create_email_content.py
# Author: sslyu
# mail: 1573077420@qq.com
# Created Time:Tue 04 Dec 2018 10:26:22 AM CST
#########################################################################

import numpy as np
import pandas as pd
import sys
import os

if os.getcwd().startswith('/DATA'):
    dirname = "/".join(os.getcwd().split()[-2:])
    dirname = "/DATA/sslyu/Project/Genetics_children_hospital/2018/"+dirname
else:
    dirname = os.getcwd()

for f in os.listdir('.'):
    if 'eml' in f:
        sn = f.split(']')[1].split('.eml')[0]
        #sn = f[1:-18]

ofile = open('email.txt','w')
ofile.write('以下为儿童遗传病分析结果：\n批次：'+sn+'结果\n\n位置\n'+dirname+'\n\n祝好，\n吕珊珊\n\n')
ofile.close()

