#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
version:3.0
by:lss
"""

import numpy as np
import pandas as pd
import sys
reload(sys)
sys.setdefaultencoding("utf8")
import os

sys.path.append("/DATA/sslyu/trio_BWA-GATK_3.0/")
from lib.prepare import *
from lib.GAF import *
from lib.xlsx import *

kwargs_raw = sys.argv[1:]
kwargs = {'-sn':'', '-l':''}
for i in range(len(kwargs_raw)):
    for j in kwargs.keys():
        if kwargs_raw[i] == j:
            kwargs.update({j:kwargs_raw[i+1]})
            
#read sample list
if(kwargs.get('-sn')!=''):
    list=[]
    list.append(kwargs.get('-sn'))
    print "\nNote:sigle sample mode"
    print unicode(list[0], 'utf-8')
    #sep = raw_input('please enter the split characters:')
    if u'分析任务单' in kwargs.get('-sn'):
        sep = u'分析任务单'
    else:
        sep = u'产品信息表'
    seperator = []
    seperator.append(sep)
else:
    if kwargs.get('-l')!='':
        list=pd.read_table(kwargs.get('-l'),header=None)[0].tolist()
        print "\nNote:list samples mode"
        print list
        sep = raw_input('please enter the split characters of each file, seperated by space:')
        seperator = sep.split(' ')

if len(list) > 0:
    for sn in list:
        excel = pd.read_excel(sn)
        excel = append_sample_txt(excel)
        merge_columns = excel.columns.tolist()
        d1 = append_gender(excel)
        d1 = append_relation(d1)
        d1 = append_pedigree(d1)
        d1 = append_phenotype(d1)
        print d1[[u'姓名', 'sample']]
        d1 = append_father_mother(d1)

        strategy = os.getcwd().split('/')[-1].split('_')[-1]
        #print 'strategy: '+strategy
        #truth = raw_input('Is the strategy right? If not, please enter it (WES IDT-PANEL PANEL Agilent-wes ... ): ')
        #if truth != 'yes':
        #    strategy = truth

        pedigree_dict = generate_pedigree(d1)
        import collections
        trio = [item for item, count in collections.Counter(d1['familyname'].tolist()).items() if count > 1]
        single = [item for item, count in collections.Counter(d1['familyname'].tolist()).items() if count == 1]
        ped_v = []
        for j in d1['sample']:
        	if j in single:
        		ped_v.append(0)
        	else:
        		ped_v.append(1)
        d1['pedigree'] = ped_v
        d1[u'分析流程'] = d1['pedigree'].map(lambda x:'single' if x==0 or x=='nan' else 'trio')
        
        if os.path.exists('all.qcsum'):
            qcsum = pd.read_table('all.qcsum', header=None)
            qcsum.columns = ['sample', '1Xcov', '20Xcov', 'AvgDepth', 'duplication']
            qcsum.index = qcsum['sample']
        else:
            qcsum = pd.DataFrame(columns=['sample', '1Xcov', '20Xcov', 'AvgDepth', 'duplication'])
            for j in d1['sample']:
            	fn = j+'/2_mapping/'+j+'.qcsum'
            	df = pd.read_table(fn, dtype=str)
            	qcsum = pd.concat([qcsum,df], ignore_index=True)
            qcsum.index = qcsum['sample']
            qcsum.to_csv('all.qcsum', sep="\t", index=None)
        
        if os.path.exists('all.qcvcf'):
            qcvcf= pd.read_table('all.qcvcf', header=None)
            qcvcf.columns = ['sample', 'indel_genomic', u'snp_genomic', u'ts/tv_genomic',u'snp_exonic', u'ts/tv_exonic']
            qcvcf.index = qcvcf['sample']
        else:
            qcvcf = pd.DataFrame(columns=['sample', 'indel_genomic', u'snp_genomic', u'ts/tv_genomic',u'snp_exonic', u'ts/tv_exonic'])
            qcvcf_list = []
            d1.index = d1['sample']
            for j in d1['sample'].tolist():
                if j in single:
                    qcvcf_list.append(j)
                elif j in trio:
                    qcvcf_list.append(j)
                else:
                    continue
            for j in qcvcf_list:
            	if j in single:
            		fn = j+'/3_variants/'+j+'.qcvcf'
            	else:
            		fn = 'trio/'+j+'/triotmp/'+j+'.qcvcf'
            	df = pd.read_table(fn, dtype=str)
            	qcvcf = pd.concat([qcvcf,df], ignore_index=True)
            qcvcf.index = qcvcf['sample']
        
        #print d1.columns
        #wrong_header = raw_input('Please enter the colnames of reads, bases and Q30: ')
        #wrong_header = wrong_header.split(', ')
        #print d1[wrong_header]
        #right_header = ['Reads','bases(Mb)','Q30']
        #for i in range(len(wrong_header)):
        #    d1.rename(columns={wrong_header[i]:right_header[i]}, inplace=True)
        #    merge_columns[merge_columns.index(wrong_header[i])] = right_header[i]
        #if os.path.exists('all_fqstat'):
        #    d4 = pd.read_table('all_fqstat', header=None, index_col=0)
        #    d1[[u'Reads', u'bases(Mb)', u'Q30']] = d4[[1,2,3]]
        #    print d1[[u'Reads', u'bases(Mb)', u'Q30']]
        d1.index = d1['sample']	

        #merge_columns = excel.columns.tolist()
        d1_merge = d1[merge_columns]
        qcsum_merge = qcsum.drop(u'duplication', axis=1)
        jiaofuxinxi_table = pd.merge(d1_merge, qcsum_merge, on="sample") # notice the type of the value
        jiaofuxinxi_full = pd.merge(d1, qcsum_merge, on="sample")

        l2 = [u'策略', 'vcf', 'QC']
        l3 = [strategy, 'available', u'合格']
        for j in range(len(l2)):
            if l2[j] not in jiaofuxinxi_table.columns:
                jiaofuxinxi_table.insert(len(jiaofuxinxi_table.columns.tolist()), l2[j], l3[j])
                jiaofuxinxi_full.insert(len(jiaofuxinxi_full.columns.tolist()), l2[j], l3[j])
        jiaofuxinxi_table = jiaofuxinxi_table.drop('sample', axis=1)
        print sn.split(seperator[list.index(sn)])[0]
        print sn.split(seperator[list.index(sn)])[1]
        #ofile = sn.split(seperator[list.index(sn)])[0]+'交付信息表'+sn.split(seperator[list.index(sn)])[1]
        tmp = sn.split(seperator[list.index(sn)])[1].split(')')
        tmp1 = tmp[0]+u'）'+tmp[1]
        ofile = sn.split(seperator[list.index(sn)])[0]+'交付信息表'+tmp1
        wb=pd.ExcelWriter(ofile,engine='openpyxl')
        jiaofuxinxi_table.to_excel(wb,index=False)
        wb.save()
        
        #navicat_truth = raw_input('Generate navicat excel?  ')
        navicat_truth = "yes"
        if navicat_truth != 'no':
            header = pd.read_excel('/DATA/sslyu/trio_BWA-GATK_3.0/doc/navicat_header.xlsx')
            total_table = pd.read_excel('/DATA/sslyu/trio_BWA-GATK_3.0/bk/navicat_bak.xlsx')
            colname3 = header.columns
            sample_n = d1.shape[0]
            #start = raw_input('please enter the start number of this batch of samples in the database:')
            start = int(total_table['order'].tolist()[-1])+1
            header[u'order'] = range(int(start), int(start)+sample_n)
            jiaofuxinxi_full.rename(columns={u'原始样本ID':u'样本ID', u'科室':u'申请科室'}, inplace=True)
            header[[u'批次', u'样本编号', u'样本ID', u'策略', u'姓名', u'性别', u'年龄',u'样本间关系', u'申请科室',u'exp_reads',u'exp_bases(Mb)', u'exp_Q30', u'exp_浓度(ng/ul)', u'exp_体积(ul)',u'exp_质量(ug)', u'exp_OD260/OD280', u'exp_OD260/OD230', u'exp_QC结论',u'ana_qc_comment',u'ana_pepline', u'ana_vcf']] = jiaofuxinxi_full[[u'批次', u'样本编号', u'样本ID', u'策略', u'姓名', u'性别', u'年龄',u'样本间关系', u'申请科室',u'reads',u'bases(Mb)',u'Q30',u'浓度(ng/ul)',u'体积(ul)',u'质量(ug)',u'OD260/OD280',u'OD260/OD230',u'QC',u'QC',u'分析流程',u'vcf']]
            header.index = d1.index # very important
            header[u'家系关系'] = None
            header[u'家系关系'] = d1[u'家系关系']
            header[[u'process_送样日期', u'process_上机日期', u'process_下机日期']] = None
            header[[u'process_送样日期', u'process_上机日期', u'process_下机日期']] = d1[[u'收样日期', u'上机日期',u'下机日期']]
            header[[u'ana_coverage1X(%)',u'ana_coverage20X(%)', u'ana_avg_depth(X)', u'ana_duplication(%)']] = qcsum[[u'1Xcov', u'20Xcov', u'AvgDepth', u'duplication']]
            header[[u'ana_indel(genomic)', u'ana_snp(genomic)', u'ana_ts/tv(genomic)',u'ana_snp(exonic)', u'ana_ts/tv(exonic)']] = qcvcf[[u'indel_genomic', u'snp_genomic', u'ts/tv_genomic',u'snp_exonic', u'ts/tv_exonic']]
            #annotation_ver = raw_input('Please enter the version of annotation: ')
            annotation_ver = "3.0"
            header[u'ana_annotation'] = annotation_ver
            header[u'ana_reporter'] = u'吕珊珊'
            bk_rawdata = u'/anaData/anaData004/children_hos_genetic/rawdata/2018/'+strategy+'/'
            for i in d1[u'样本编号']:
                os.system('ls '+bk_rawdata+i)
            print bk_rawdata
            #bk_rawdata_truth = raw_input('Is the rawdata backup dirname correct? If not, please enter the correct dirname of the bk_rawdata: ')
            bk_rawdata_truth = "yes"
            if bk_rawdata_truth != 'yes':
                bk_rawdata = bk_rawdata_truth
            header[u'bk_rawdata'] = bk_rawdata
            bk_reports = u'/DATA/BPshare/opm/遗传病-交大附属儿医/2018/'+'/'.join(os.getcwd().split('/')[-2:])
            #os.system('ls '+bk_reports)
            #print bk_reports
            #bk_reprots_truth = raw_input('Is the dirname correct? If not, please enter the correct dirname of the bk_reports:')
            #if bk_reprots_truth != 'yes':
            #    bk_reports = bk_reports_truth
            header[u'bk_reports'] = bk_reports
            result = total_table.append(header)
            
            #ofile = sn.split(seperator[list.index(sn)])[0]+'navicat'+sn.split(seperator[list.index(sn)])[1]
            tmp = sn.split(seperator[list.index(sn)])[1].split(')')
            tmp1 = tmp[0]+u'）'+tmp[1]
            ofile = sn.split(seperator[list.index(sn)])[0]+'navicat'+tmp1
            wb=pd.ExcelWriter(ofile, engine='openpyxl')
            header.to_excel(wb,index=False)
            wb.save()
            ofile1 = '/DATA/sslyu/trio_BWA-GATK_3.0/doc/navicat.xlsx'
            wb1=pd.ExcelWriter(ofile1, engine='openpyxl')
            result.to_excel(wb1,index=False)
            wb1.save()
            os.system('cp /DATA/sslyu/trio_BWA-GATK_3.0/doc/navicat.xlsx /DATA/sslyu/trio_BWA-GATK_3.0/bk/navicat_bak.xlsx')
