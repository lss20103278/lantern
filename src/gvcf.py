#!/opt/anaconda2/bin/python
import sys
import os


sn = sys.argv[1]
#with open('gvcf/'+sn+'.g.vcf','w') as t:
#    with open('gvcf/'+sn+'.raw.g.vcf') as f:
with open(sn+'.g.vcf','w') as t:
    with open(sn+'.raw.g.vcf') as f:
        lst = f.readlines()
        lst.insert(1,'##FORMAT=<ID=Q,Number=1,Type=Float,Description="Quality">\n')
        for i in lst:
            if i.startswith('#'):
                t.write(i)
                continue
            line = i.strip().split('\t')
            if line[5] != '.':
                line[-2] = line[-2]+':'+'Q'
                line[-1] = line[-1]+':'+line[5]
	    t.write('\t'.join(line)+'\n')
