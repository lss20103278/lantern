#!/usr/bin/env python
# -*- coding: utf-8 -*-
#########################################################################
# File Name: /home/sslyu/trio_BWA-GATK_2.7/read_excel.py
# Author: sslyu
# mail: 1573077420@qq.com
# Created Time:Tue 04 Sep 2018 03:59:26 PM CST
#########################################################################

import numpy as np
import pandas as pd
import sys
reload(sys)
sys.setdefaultencoding("utf8")
import os
import re

dataframe = pd.read_excel(sys.argv[1])
dataframe.to_csv(sys.argv[1].split('.xlsx')[0]+'.csv', index=None, sep="\t")
