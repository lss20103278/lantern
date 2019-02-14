#!/usr/bin/env python
# -*- coding: utf-8 -*-
#########################################################################
# File Name: txt2excel.py
# Author: sslyu
# mail: 1573077420@qq.com
# Created Time:Mon 15 Oct 2018 04:20:20 PM CST
#########################################################################

import numpy as np
import pandas as pd
import sys

infile = sys.argv[1]
df = pd.read_table(infile, dtype=str)
ofile = infile+'.xlsx'
wb=pd.ExcelWriter(ofile, engine='openpyxl')
df.to_excel(wb,index=False)
wb.save()
