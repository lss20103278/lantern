#!/opt/anaconda2/bin/python

import sys
import re


"""

"""




""" Functions """
#annovar format
def convert_2_annovar(pos,ref,alt):
    if(len(ref)>len(alt)):		#del
        newref = ref[len(alt):]
        newpos = pos + len(alt)
        newalt = "-"
        pos2 = newpos + len(newref) - 1
    elif(len(ref)<len(alt)):		#insert
        newalt = alt[len(ref):]
        #newpos = pos + len(ref) 
        newref = "-"
        newpos = pos
        pos2 = newpos
    else:
        newref = ref[-1]
        newalt = alt[-1]
        newpos = pos + len(ref) - 1
        pos2 = newpos
    return([newpos,pos2,newref,newalt])

#extract_gt
def extract_gt(ref,alt,x):
    gt = x.split(":")[0]
    if(gt=="./."):
        return("./.")
    elif(gt=="0/0"):
        return(ref+"/"+ref)
    elif(gt=="0/1"):
        return(ref+"/"+alt)
    elif(gt=="1/1"):
        return(alt+"/"+alt)
    elif(gt=="1/0"):
    	return(alt+"/"+ref)
    else:
        return(gt)

#extract_het
def extract_het(x):
    gt = x.split(":")[0]
    if(gt=="0/0" or gt=="1/1"):
        return("hom")
    elif(gt=="0/1"):
        return("het")
    else:
        return("unexact")

#extract_depth
def extract_depth(x):
    if(x=="./."):
        return(".")
    else:
        return(x.split(":")[2])
        
#extract_qual
def extract_qual(x):
    if(x=="./."):
        return(".")
    else:
        return(x.split(":")[-1])

""" Func END """



""" main: """
sn = str(sys.argv[1])
infile = open("segtrio/"+sn+".mendel.vcf",'r').readlines()
ofile = open("segtrio/"+sn+".avinput",'w')
header = open("segtrio/"+sn+".avinput.info","w")
for line in infile:
    if(re.match("##",line)):
        continue
    if(re.match("#CHROM",line)):
        samples=line.split("FORMAT")[1].split()
        headerinfo=""
        for i in samples: headerinfo += "\tGTinfo_" + i
        for i in samples: headerinfo += "\tGT_" + i
        for i in samples: headerinfo += "\tQual_"+i  # add trio qual
        for i in samples: headerinfo += "\tDepth_"+i # add trio depth 
        #headerinfo += "\tHeterozygosity\tQual\tDepth\tSEGSCORE\tDENOVO\tDENOVO_Heterozygosity\tlowGQ"
        headerinfo += "\tHeterozygosity\tQual\tTrio_GTinfo\tSEGSCORE\tDENOVO\tDENOVO_Heterozygosity\tlowGQ"
        header.write("%s" % headerinfo)
        continue
    
    line = line.strip().split("\t")
    gtinfo = ""
    gt = "" 
    depths=""
    quals=""
    Trio_GTinfo = []
    #if len(set(map(lambda x:x.split(':')[0],line[9:12]))) == 1:continue    #filter the same GT with trio
    #if set(map(lambda x:x.split(':')[0],line[9:12])) == set(['1/0']):continue  #filter trio 1/1
    for i in line[9:] : gtinfo += "\t" + i ; gt += "\t" + extract_gt(line[3],line[4],i);quals +="\t"+extract_qual(i);depths +="\t"+extract_depth(i)
    for i in line[9:]: Trio_GTinfo.append(i.split(':')[0])  # Trio_GTinfo
    
    for j in range(len(samples)):
        if(samples[j]==sn):break
    
    tar_gtinfo = gtinfo.strip().split('\t')[j]
    
    newpos = convert_2_annovar(int(line[1]),line[3],line[4])
    if(newpos[3]=="*") : continue   #fix: ALT==*
    if(line[9].split(":")[0]=="0/0") : continue #fix: gt='0/0'
        
    if(line[7].find("loConfDeNovo")!=-1):
        denovo="loConf"
        #denovo_heterozygosity=extract_het(line[9]) # Denovo heterozygosity
        denovo_heterozygosity=extract_het(tar_gtinfo)
    elif(line[7].find("hiConfDeNovo")!=-1):
        denovo="hiConf"
        #denovo_heterozygosity=extract_het(line[9]) # Denovo heterozygosity
        denovo_heterozygosity=extract_het(tar_gtinfo)
    else:
        denovo="-"
        denovo_heterozygosity="-" # Denovo heterozygosity
    
    if(gtinfo.find("lowGQ")!=-1):
        lowGQ="True"
    else:
        lowGQ="-"
    
    newline = "chr" + str(line[0]).upper().replace("CHR","",1) + "\t" + str(newpos[0]) +"\t" + str(newpos[1]) + "\t" + newpos[2] + "\t" + newpos[3] + "\t" + gtinfo + gt + quals+depths+"\t" +extract_het(tar_gtinfo)+"\t"+line[5]+"\t"+":".join(Trio_GTinfo)+"\t" + line[7].split("SEGSCORE=")[1].split(";")[0] + "\t" + denovo + "\t"+denovo_heterozygosity+"\t"+lowGQ
    ofile.write("%s\n" % newline)
ofile.close()
header.close()





