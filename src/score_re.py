import sys
import re
import os
import pandas as pd
##read parameters of cmd line
kwargs_raw = sys.argv[1:]
kwargs={'-edge':'ucsc_hg19_refGene.csv','-l':'list','-AF':"combine_run1.snp.AF",'-sn':'', '-c_f_m':'', '-thisminedge':'15'}
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

p1=re.compile("Pathogenic")  
#p2=re.compile("Likely pathogenic") #previous version annovar 2015-12-14
##p3=re.compile("risk factor") 
##p4=re.compile("association") 
#p3=re.compile("Likely benign") #previous version annovar 2015-12-14
p4=re.compile("Benign") #previous version annovar 2015-12-14
#p4_2 = re.compile("Uncertain significance") #previous version annovar 2015-12-14
p5=re.compile("CLIN_pathogenic")
p6=re.compile("CLIN_likely_pathogenic")
p7=re.compile("CLIN_risk_factor")
p8=re.compile("CLIN_association")

p2 = re.compile("Likely_pathogenic") # annovar 2018
p3 = re.compile("Likely_benign") # annovar 2018
p4_2 = re.compile("Uncertain_significance") # annovar 2018

s1=[4,3,2,1]  #coefficient?
#inhgmd=open("hgmd.vcf","r")
#hgmd={}
#for line in inhgmd.readlines():
#  hgmd["chr"+line.split('\t')[0]+','+line.split('\t')[1]]=1

######## related to DiseaseID, Disease, Disease type and InheritanceStatus #########
diseaseid={} # column DiseaseID
disease={} # column Disease
typedisease={} # column Disease type
inheri={} # column InheritanceStatus
indisease=open("/DATA/sslyu/trio_BWA-GATK_3.0/doc/gene_disease.txt","r")
for line in indisease.readlines():
  if line.split('\t')[3] in disease.keys():
    disease[line.split('\t')[3]]=disease[line.split('\t')[3]]+','+line.split('\t')[1]
  else:
    disease[line.split('\t')[3]]=line.split('\t')[1]

  if line.split('\t')[3] in diseaseid.keys():
    diseaseid[line.split('\t')[3]]=diseaseid[line.split('\t')[3]]+','+line.split('\t')[0]
  else:
    diseaseid[line.split('\t')[3]]=line.split('\t')[0]

  if line.split('\t')[3] in typedisease.keys():
    typedisease[line.split('\t')[3]]=typedisease[line.split('\t')[3]]+','+line.split('\t')[2]
  else:
    typedisease[line.split('\t')[3]]=line.split('\t')[2]
  
  if line.split('\t')[3] in inheri.keys():
    inheri[line.split('\t')[3]]=inheri[line.split('\t')[3]]+','+line.strip().split('\t')[4]
  else:
    inheri[line.split('\t')[3]]=line.strip().split('\t')[4]

for sn in list:
  sn = str(sn)
  print sn
  inline=open(sn+".ann.hg19_multianno.txt","r")
  outfile=open(sn+".score","w")
  outfile.write(open(sn+".ann.hg19_multianno.txt","r").readline().strip()+'\t'+"score_clinvar"+'\t'+"score_ensembl"+'\t'+"score_pre"+'\t'+"score_HGMD"+'\t'+"score_total"+'\t'+"DiseaseID"+'\t'+"Disease"+'\t'+"Disease type"+'\t'+"InheritanceStatus"+'\n')
  for line in inline.readlines()[1:]:
    s2=[0,0,0,0]
    s3=[0,0,0,0]
    scorepre=0
    scorehgmd=0
    
####### calculate score_clinvar ###########################################
    #if p1.search(line.strip().split('\t')[119]): #Pathogenic, annovar 2015
    #  s2[0]=1
    #if p2.search(line.strip().split('\t')[119]): #Likely pathogenic
    #  s2[1]=1
    #if p3.search(line.strip().split('\t')[119]): #Likely benign
    #  s2[2]=1
    #if p4.search(line.strip().split('\t')[119]): #Benign
    #  s2[3]=1
    #if p4_2.search(line.strip().split('\t')[119]): #Uncertain significance
    #  s2[3]=1
    if p1.search(line.strip().split('\t')[123]): #Pathogenic, annovar 2018 line[123]: column CLINSIG
      s2[0]=1
    if p2.search(line.strip().split('\t')[123]): #Likely pathogenic
      s2[1]=1
    if p3.search(line.strip().split('\t')[123]): #Likely benign
      s2[2]=1
    if p4.search(line.strip().split('\t')[123]): #Benign
      s2[3]=1
    if p4_2.search(line.strip().split('\t')[123]): #Uncertain significance
      s2[3]=1
    scoreclin=max(map(lambda x, y: x*y,  s1, s2))

####### calculate score_ensembl ############################################ 
    if p5.search(line.strip().split('\t')[124]): #CLIN_pathogenic line[124]: column ensembl
      s3[0]=1
    if p6.search(line.strip().split('\t')[124]): #CLIN_likely_pathogenic
      s3[1]=1
    if p7.search(line.strip().split('\t')[124]): #CLIN_risk_factor
      s3[2]=1
    if p8.search(line.strip().split('\t')[124]): #CLIN_association
      s3[3]=1
    scoreensembl=max(map(lambda x, y: x*y,  s1, s3))

###### calculate score_pre ##################################################
    if line.strip().split('\t')[25]=="D": # SIFT_pred
      scorepre=scorepre+1
    if line.strip().split('\t')[31]=="D": # Polyphen2_HVAR_pred
      scorepre=scorepre+1
    if line.strip().split('\t')[37]=="D": # MutationTaster_pred
      scorepre=scorepre+1
    if line.strip().split('\t')[43]=="D": # FATHMM_pred
      scorepre=scorepre+1
    if line.strip().split('\t')[57]=="D": # M-CAP_pred
    	scorepre=scorepre+1
    if line.strip().split('\t')[89] != 'NA': # REVEL
      if float(line.strip().split('\t')[89])>0.5:
        scorepre=scorepre+1

###### calculate score_HGMD #################################################
#    if line.split('\t')[0]+','+line.split('\t')[1] in hgmd.keys():
    confhgmd=['NA','','.','-']
    if (line.split('\t')[125] not in confhgmd): # hgmd
      scorehgmd=5
    
    if line.split('\t')[6].split(',')[0] in diseaseid.keys():
      outfile.write(line.strip()+'\t'+str(scoreclin)+'\t'+str(scoreensembl)+'\t'+str(scorepre)+'\t'+str(scorehgmd)+'\t'+str(scoreclin+scoreensembl+scorepre+scorehgmd)+'\t'+diseaseid[line.split('\t')[6].split(',')[0]]+'\t'+disease[line.split('\t')[6].split(',')[0]]+'\t'+typedisease[line.split('\t')[6].split(',')[0]]+'\t'+inheri[line.split('\t')[6].split(',')[0]]+'\n') # str(scoreclin+scoreensembl+scorepre+scorehgmd) calculate score_total
    else:
      outfile.write(line.strip()+'\t'+str(scoreclin)+'\t'+str(scoreensembl)+'\t'+str(scorepre)+'\t'+str(scorehgmd)+'\t'+str(scoreclin+scoreensembl+scorepre+scorehgmd)+'\t'+"NA"+'\t'+"NA"+'\t'+"NA"+'\t'+"NA"+'\n') # str(scoreclin+scoreensembl+scorepre+scorehgmd) calculate score_total
      
