    ds=pd.read_table(sn+".score",sep="\t",dtype=str)
    df_maf = pd.read_csv(sn+'.maf',sep = '\t',header = 0,dtype = str)
    df_maf = df_maf[['Chromosome','Start_Position','End_Position','Reference_Allele','Tumor_Seq_Allele2','HGVSc','HGVSp']]
    df_maf.columns = ['Chr','Start','End','Ref','Alt','HGVSc','HGVSp']
    df_maf_tmp = df_maf['Ref'].apply(lambda x:1 if x=='-' else 0)
    df_maf['End'] = df_maf['End'].astype(int)-df_maf_tmp.astype(int)
    df_maf['End']=df_maf['End'].astype(str)
    ds = pd.merge(ds,df_maf,how = 'left',on = ['Chr','Start','End','Ref','Alt'])
    ############
    df_link = pd.read_csv(sn+'.link',sep = '\t',header = None,dtype = str)
    #discriminate single or pedigree
    a = []
    for i in ds.columns:
    	if 'GTinfo_' in i:
    		a.append(i[7:])
    if len(a) == 0:
        a.append(sn) # sn should be the same with the name in sn.link
    	#name = raw_input( 'Is '+sn+' is the name of the sample? If not, please enter the correct name; if yes, please press enter:')
    	#if name == '':
    	#	name = sn
    	#a.append(name)
    #df_link = df_link[[0,1,2,3,4,6,8,9]] # We need the 14th column to calculate the percentage of the reads containing alleles
    df_link[0] = df_link[0].apply(lambda x:x if x.startswith('chr') else 'chr'+x)
    df_link['link'] = df_link[0]+'-'+df_link[6]+'-'+df_link[8]+'-'+df_link[9]
    df_link = df_link.drop([6,8,9],axis = 1)		
    if len(a) == 1:
        b = []
        df_link_columns = [0,1,2,3,4,14,'link']
        df_link = df_link[df_link_columns]
        df_link.columns = ['Chr','Start','End','Ref','Alt',a[0],'link']
        ds = pd.merge(ds,df_link,how = 'left',on = ['Chr','Start','End','Ref','Alt'])
        for j in ds[a[0]]:
    		try:
    			c = j.split(':')[1].split(',')
    			d = float(c[0])+float(c[1])
    			b.append('{:.2f}'.format(float(c[1])/d*100)+'%')
    		except:
    			b.append(np.nan)
        ds.insert(10,a[0]+'_allele_percentage',b)
        ds.drop(a[0],axis=1)

    else:
        for i in a:
            b = []
            for j in ds['GTinfo_'+i]:
                try:
                    c = j.split(':')[1].split(',')
                    d = float(c[0])+float(c[1])
                    b.append('{:.2f}'.format(float(c[1])/d*100)+'%')
                except:
                    b.append(np.nan)
            ds.insert(10+a.index(i),i+'_allele_percentage',b)
    	#df_link_columns = [0,1,2,3,4,u'link']
    	#df_link = df_link[df_link_columns] # do not work, don't know why?
        df_link = df_link[[0,1,2,3,4,u'link']]
        df_link.columns = [u'Chr',u'Start',u'End',u'Ref',u'Alt',u'link']					
        ds = pd.merge(ds,df_link,how = 'left',on = [u'Chr',u'Start',u'End',u'Ref',u'Alt'])
    
    df_cadd = pd.read_csv(sn+'.CADD',sep = '\t',header = 0,dtype = str)
    df_cadd['link'] = df_cadd['#CHROM'].apply(lambda x:'chr'+x)+'-'+df_cadd['POS']+'-'+df_cadd['REF']+'-'+df_cadd['ALT']
    df_cadd = df_cadd.drop(['#CHROM','POS','REF','ALT'],axis = 1)
    df_cadd.columns = ['CADD_RawScore','CADD_PHRED','link']
    ds = pd.merge(ds,df_cadd,how = 'left',on = ['link'])
    
################ deal with the deletion problem to see if the deletion of nuclei acids affects amnio acids ##########################
    ds_del = ds[ds['Alt']=='-']  #ds_del: Alt == '-' what situation?
    #ds_del['conf'] = ds_del['Ref'].apply(lambda x:len(x)) #SettingWithCopyWarning:Try using .loc[row_indexer,col_indexer] = value instead
    value = ds_del['Ref'].apply(lambda x:len(x))
    ds_del.loc[:,'conf'] = value 
    ds_del_length = ds_del[ds_del['conf']>6] # Alt == '-', Ref has more than 6 nt ?
    ds_del_length = ds_del_length.fillna('-')
    
    ds_del_length_fs = ds_del_length[ds_del_length['AAChange.refGene'].str.contains('[\d]fs')] # number+fs what does it mean? fs:frameshift
    #ds_del_length_fs['Deletion_Problem'] = 'note' # add column 'Deletion_Problem' SettingWithCopyWarning:Try using .loc[row_indexer,col_indexer] = value instead
    ds_del_length_fs.loc[:,'Deletion_Problem'] = 'note'
    ds_del_length = ds_del_length[~ds_del_length['AAChange.refGene'].str.contains('[\d]fs')] #  AAChange.refGene of rows don't contain [\d]fs
    #ds_del_length = ds_del_length[ds_del_length['AAChange.refGene']!='-']
    ds_del_length = ds_del_length[ds_del_length['AAChange.refGene'].str.contains('p.*del')] # AAChagne.refGene of rows doesn't contain [\d]fs but contains p.*del
    AA_num = map(lambda x:int(x[1])-int(x[0]),ds_del_length['AAChange.refGene'].apply(lambda x:x.split(',')[0].split('p.')[1].split('del')[0].split('_'))) # extract the number between p. and del in p.*del Amino acid number
    ds_del_length.index = range(len(ds_del_length))
    Deletion_Problem = []
    for i in ds_del_length.index:
    	if int(ds_del_length.ix[i,'conf']) != AA_num[i]*3 and int(ds_del_length.ix[i,'conf']) != (AA_num[i]+1)*3: # deletion affect the codon (AA_num[i]+1)*3 ? 3 nucleic acids code 1 amino acid
    		Deletion_Problem.append('note')
    	else:
    		Deletion_Problem.append('-')
    ds_del_length['Deletion_Problem'] = pd.Series(Deletion_Problem,index = ds_del_length.index)
    ds_del_length_all = pd.concat([ds_del_length_fs,ds_del_length])
    ds_del_length_all = ds_del_length_all[['Chr','Start','End','Ref','Alt','Deletion_Problem']]
    ds = pd.merge(ds,ds_del_length_all,how = 'left',on = ['Chr','Start','End','Ref','Alt'])

    #############################
    ds[['Dead_Zone','Problem_High','Problem_Low']]=ds[['Dead_Zone','Problem_High','Problem_Low']].replace('Name=Y','Y') #
    ds['HI_Predictions'] = ds['HI_Predictions'].fillna('-').apply(lambda x:x if x=='-' else x.split('=')[1]) # 
    #ds.insert(ds.columns.tolist().index('InterVar(automated)')+1,'PVS1/PS1/PS2/PS3/PS4/PM1/PM2/PM3/PM4/PM5/PM6/PP1/PP2/PP3/PP4/PP5/BA1/BS1/BS2/BS3/BS4/BP1/BP2/BP3/BP4/BP5/BP6/BP7',np.nan)
    ds.insert(ds.columns.tolist().index('InterVar_automated')+1,'PVS1/PS1/PS2/PS3/PS4/PM1/PM2/PM3/PM4/PM5/PM6/PP1/PP2/PP3/PP4/PP5/BA1/BS1/BS2/BS3/BS4/BP1/BP2/BP3/BP4/BP5/BP6/BP7',np.nan) # annovar 20180416 InterVar(automated)->InterVar_automated
    ds['PVS1/PS1/PS2/PS3/PS4/PM1/PM2/PM3/PM4/PM5/PM6/PP1/PP2/PP3/PP4/PP5/BA1/BS1/BS2/BS3/BS4/BP1/BP2/BP3/BP4/BP5/BP6/BP7'] = \
    	ds['PVS1'].apply(lambda x:str(x))+'/'+ds['PS1'].apply(lambda x:str(x))+'/'+\
    	ds['PS2'].apply(lambda x: str(x)) + '/' + ds['PS3'].apply(lambda x: str(x)) + '/' + ds['PS4'].apply(lambda x: str(x)) + '/'+ \
    	ds['PM1'].apply(lambda x: str(x)) + '/' + ds['PM2'].apply(lambda x: str(x)) + '/' + ds['PM3'].apply(lambda x: str(x)) + '/'+ \
    	ds['PM4'].apply(lambda x: str(x)) + '/' + ds['PM5'].apply(lambda x: str(x)) + '/' + ds['PM6'].apply(lambda x: str(x)) + '/'+ \
    	ds['PP1'].apply(lambda x: str(x)) + '/' + ds['PP2'].apply(lambda x: str(x)) + '/' + ds['PP3'].apply(lambda x: str(x)) + '/'+ \
    	ds['PP4'].apply(lambda x: str(x)) + '/' + ds['PP5'].apply(lambda x: str(x)) + '/' + ds['BA1'].apply(lambda x: str(x)) + '/'+ \
    	ds['BS1'].apply(lambda x: str(x)) + '/' + ds['BS2'].apply(lambda x: str(x)) + '/' + ds['BS3'].apply(lambda x: str(x)) + '/'+ \
    	ds['BS4'].apply(lambda x: str(x)) + '/' + ds['BP1'].apply(lambda x: str(x)) + '/' + ds['BP2'].apply(lambda x: str(x)) + '/'+ \
    	ds['BP3'].apply(lambda x: str(x)) + '/' + ds['BP4'].apply(lambda x: str(x)) + '/' + ds['BP5'].apply(lambda x: str(x)) + '/'+ds['BP6'].apply(lambda x: str(x)) + '/' + ds['BP7'].apply(lambda x: str(x))
    
    #ds = ds.drop(ds.columns.tolist()[ds.columns.tolist().index('PVS1'):ds.columns.tolist().index('CLINSIG')],axis=1)
    ds = ds.drop(ds.columns.tolist()[ds.columns.tolist().index('PVS1'):ds.columns.tolist().index('CLNALLELEID')],axis=1) # clinvar_20180603 the order of clinvar columns changes
    ds['Start']=ds['Start'].astype(long)
    
    #chr tag confirm
    if("chr" in ds.iloc[0]['Chr']):
    	print "chrid_check : ok"
    else:
    	ds['Chr']="chr"+ds['Chr']
    	print "chrid_rewrite : ok"
    """
    check=[]
    for i in ds.index:
    	line=ds.loc[i]['Gene.refGene'].split(";")
    	tmp=[id for id in line if id in genelist]
    	if(tmp):
    		check.append(True)
    	else:
    		check.append(False)
    
    ds=ds[check]
    """
    
    ##step1:filt1+annotation_AF+annotation_OMIM
    print "\nStep1:\n--filting by snp/indel location"
    
    ds.insert(7,'min_edge',[np.nan]*len(ds))
    tar=['exonic','splicing','exonic;splicing']
    ds_a=ds[ds['Func.refGene'].isin(tar)].copy() # .copy() modifications to the data or indices of the copy will not be reflected in the original object
    ds_b=ds[~ds['Func.refGene'].isin(tar)].copy() # Func.refGene is neither exonic nor splicing 
    ds_b.index=range(len(ds_b))
    for i in ds_b.index:
    	line=ds_b.loc[i]
    	chrom=line['Chr']
    	pos=int(line['Start'])
    	genes=line['Gene.refGene'].replace(",",";").split(";")
    	
    	#sub of edge by chr
    	""" bugs maybe: the end should +1? """
    	ex_sub=ex_edge.loc[ex_init.get(chrom):ex_end.get(chrom)] #??		
    	ex_sub=ex_edge[ex_edge['symbol'].isin(genes)] #??
    	refedge='::edge::'
    	#thisminedge=11
    	thisminedge=thisminedge_set
    	confheader=''
    	for j in ex_sub.index:						  
    		ex_subline=ex_sub.loc[j]     
    		eCount=ex_subline.exonCount
    		eStart=ex_subline.exonStarts.split()
    		for k in range(len(eStart)):eStart[k]=int(eStart[k])+1
    		eEnd=ex_subline.exonEnds.split()
    		for k in range(len(eEnd)):eEnd[k]=int(eEnd[k])
    		if(ex_subline.strand=='+'):
    			for k in range(eCount):
    				if(eStart[k] > pos and eStart[k]-pos <= thisminedge_set):
    					if(eStart[k]-pos < thisminedge_set+1): thisminedge=eStart[k]-pos
    					confheader='-' # in plus strand, '-' means it is located ahead of the start of the exon
    					refedge=refedge + ex_subline['name']+":"+ex_subline.symbol+":exon"+str(k+1)+":"+str(eStart[k])+":-"+str(eStart[k]-pos)+";"
    				if(eEnd[k] < pos and pos-eEnd[k] <= thisminedge_set):
    					if(pos-eEnd[k] < thisminedge_set+1): thisminedge=pos-eEnd[k]
    					confheader='+' # in plus strand, '+' means it is located behind the end of the exon
    					refedge=refedge + ex_subline['name']+":"+ex_subline.symbol+":exon"+str(k+1)+":"+str(eEnd[k])+":+"+str(pos-eEnd[k])+";"
    		elif(ex_subline.strand=='-'):
    			for k in range(eCount):
    				if(eEnd[k] < pos and pos-eEnd[k] <= thisminedge_set):
    					if(pos-eEnd[k] < thisminedge_set+1): thisminedge=pos-eEnd[k]
    					confheader='-' # in minus strand, '-' means it is located behind the end of the exon
    					refedge=refedge + ex_subline['name'] + ":" + ex_subline.symbol + ":exon" + str(eCount-k) +":"+str(eEnd[k])+ ":-" + str(pos-eEnd[k]) + ";"
    				if(eStart[k]>pos and eStart[k]-pos <= thisminedge_set):
    					if(eStart[k]-pos < thisminedge_set+1): thisminedge=eStart[k]-pos
    					confheader='+' # in minus strand, '+' means it is located ahead of the start of the exon
    					refedge=refedge + ex_subline['name'] + ":" + ex_subline.symbol + ":exon" + str(eCount-k) + ":"+str(eStart[k])+ ":+" + str(eStart[k]-pos) + ";"
    	ds_b.loc[i,'min_edge'] = confheader+str(thisminedge)
    	ds_b.loc[i,'GeneDetail.refGene'] = str(ds_b.loc[i]['GeneDetail.refGene']) + refedge.strip(";")		  
    	if(i%1000==0 or i==ds_b.index[-1]):
    	   view_bar(i+1,len(ds_b))
    conf=[]
    for i in ds_b.index:
    	if(ds_b.loc[i]['GeneDetail.refGene'].split("::edge::")[1]==''):
    		conf.append(False)
    	else:
    		conf.append(True)
    ds_c=ds_b[conf] # it is located in the min_edge of the exons
    ds=pd.concat([ds_a,ds_c])
    ds.index=range(len(ds))
    
    
    
    #annotation:local MAF
    print "\n--adding local AF"
    ds[0] = ds['Chr'].astype(str)+'-'+ds['Start'].astype(str)+'-'+ds['End'].astype(str)+'-'+ds['Ref'].astype(str)+'-'+ds['Alt'].astype(str)
    ds = pd.merge(ds,AF,how = 'left',on=[0])
    ds = ds.drop(0,axis = 1)
    
    #annotation:OMIM_ch/OMIM_en/chpo/morbidmap/Phenotype
    print "--adding OMIM_ch/OMIM_en/CHPO/morbidmap/Phenotype"
    addOMIM_id=[]
    addOMIM_gene=[]
    addOMIM_disease_ch=[]
    addOMIM_disease_en=[]
    add_chpo_name_CH1=[]
    add_chpo_define_EN1=[]
    add_chpo_define_CH1=[]
    # update
    IP = []
    DI = []
    Inheritance = []
    Comments = []
    Gene_Locus_MIM_number = []
    Phenotype = []
    Phenotype_MIM_number = []
    add_mor_Gene = []
    add_mor_Phenotype = []
    add_mor_MIM_Number = []
    add_mor_Cyto_Location = []
    omim_phenotype = []
    ###########
    for i in ds.index:
        if(i%1000==0 or i==ds.index[-1]):
     	   view_bar(i+1,len(ds.index))
        
        targene=ds.loc[i]['Gene.refGene'].split(',')[0]
        # update 
        info_gene = DIDGene[DIDGene.Gene==targene]
        if(len(info_gene)==0):
     	   DI.append('')
        else:
     	   DI.append('YES')
        info_OMIM_Gene = OMIM_Gene_Map_Search[OMIM_Gene_Map_Search['Gene']==targene]
        if(len(info_OMIM_Gene)==0):
     	   Inheritance.append('')
     	   Comments.append('')
     	   Gene_Locus_MIM_number.append('')
     	   Phenotype.append('')
     	   Phenotype_MIM_number.append('')
        else:
     	   Inheritance.append(info_OMIM_Gene['Inheritance'].str.cat(sep = '|'))
     	   Comments.append(info_OMIM_Gene['Comments'].str.cat(sep = '|'))
     	   Gene_Locus_MIM_number.append(info_OMIM_Gene['Gene/Locus_MIM_number'].drop_duplicates().str.cat(sep = '|'))
     	   Phenotype.append(info_OMIM_Gene['Phenotype'].str.cat(sep = '||'))
     	   Phenotype_MIM_number.append(info_OMIM_Gene['Phenotype_MIM_number'].str.cat(sep='|'))		   
        info_mor = morbidmap[morbidmap.Gene == targene]
        if(len(info_mor) == 0):
     	   add_mor_Gene.append('')
     	   add_mor_Phenotype.append('')
     	   add_mor_MIM_Number.append('')
     	   add_mor_Cyto_Location.append('')
        else:
     	   add_mor_Gene.append(info_mor.iloc[0].Gene)
     	   add_mor_Phenotype.append(info_mor.iloc[0].Phenotype)
     	   add_mor_MIM_Number.append(info_mor.iloc[0].MIM_Number)
     	   add_mor_Cyto_Location.append(info_mor.iloc[0].Cyto_Location)
        ############################# 
        info_ch=OMIM_selected_ch[OMIM_selected_ch.OMIM_gene==targene]
        info_en=OMIM_selected_en[OMIM_selected_en.OMIM_gene==targene]
        if(len(info_en)==0):
     	   IP.append('')
     	   addOMIM_disease_en.append('')
        else:
     	   IP.append(info_en['IP'].str.cat(sep = '|'))
     	   addOMIM_disease_en.append(info_en['OMIM_disease'].str.cat(sep = '||'))		 
        if(len(info_ch)==0):
     	   addOMIM_id.append('')
     	   addOMIM_gene.append('')
     	   addOMIM_disease_ch.append('')
     	   add_chpo_name_CH1.append('')
     	   add_chpo_define_EN1.append('')
     	   add_chpo_define_CH1.append('')
        else:
     	   addOMIM_gene.append(info_ch.iloc[0].OMIM_gene)
     	   addOMIM_disease_ch.append(info_ch.iloc[0].OMIM_disease_ch)
     	   addOMIM_id.append(info_ch.iloc[0].OMIM_id)
     	   subphenon=phenon[phenon.infosourceid==info_ch.iloc[0].OMIM_id]
     	   if(len(subphenon)==0):
     		   add_chpo_name_CH1.append('')
     		   add_chpo_define_EN1.append('')
     		   add_chpo_define_CH1.append('')
     	   else:
     		   chpo_name=u";"
     		   chpo_define_EN=u";"
     		   chpo_define_CH=u";"
     		   for j in subphenon.index:
     			   chpo_name=chpo_name+subphenon.loc[j].chpo_name_CH1+u";"
     			   chpo_define_EN=chpo_define_EN + subphenon.loc[j].chpo_define_EN1 + u';'	
     			   chpo_define_CH=chpo_define_CH + subphenon.loc[j].chpo_define_CH1 + u";"
     		   add_chpo_name_CH1.append(chpo_name.strip(';'))
     		   add_chpo_define_EN1.append(chpo_define_EN.strip(";").strip("NA;"))
     		   add_chpo_define_CH1.append(chpo_define_CH.strip(";").strip("NA;"))
    
        #ename=''
        #for j in info_en.OMIM_disease.tolist():
     	   #ename = ename + j + "||"
        #addOMIM_disease_en.append(ename.strip('||'))
    ds['IP'] = pd.Series(IP,index = ds.index)
    ds['DI'] = pd.Series(DI,index = ds.index)
    ds['Inheritance'] = pd.Series(Inheritance,index = ds.index)
    ds['Comments'] = pd.Series(Comments,index = ds.index)
    ds['Gene_Locus_MIM_number'] = pd.Series(Gene_Locus_MIM_number,index = ds.index)
    ds['Phenotype'] = pd.Series(Phenotype,index = ds.index)
    ds['Phenotype_MIM_number'] = pd.Series(Phenotype_MIM_number,index = ds.index)
    #####		
    ds['OMIM_id']=pd.Series(addOMIM_id,index=ds.index)
    ds['OMIM_gene']=pd.Series(addOMIM_gene,index=ds.index)
    ds['OMIM_disease_ch']=pd.Series(addOMIM_disease_ch,index=ds.index)
    ds['OMIM_disease_en']=pd.Series(addOMIM_disease_en,index=ds.index)
    ds['chpo_name_CH']=pd.Series(add_chpo_name_CH1,index=ds.index)
    ds['chpo_define_EN']=pd.Series(add_chpo_define_EN1,index=ds.index)
    ds['chpo_define_CH']=pd.Series(add_chpo_define_CH1,index=ds.index)
    # update
    ds['morbidmap_Gene']=pd.Series(add_mor_Gene,index=ds.index)
    ds['morbidmap_Phenotype']=pd.Series(add_mor_Phenotype,index=ds.index)
    ds['MIM_Number']=pd.Series(add_mor_MIM_Number,index=ds.index)
    ds['Cyto_Location']=pd.Series(add_mor_Cyto_Location,index=ds.index)
    #############################
    for i in ds.index:
        taromim = ds.loc[i]['Phenotype_MIM_number']
        taromim_lst = taromim.split('|')
        tmp_lst = []
        for each_id in taromim_lst:
     	   if each_id=='':
     		   tmp_lst.append('-')
     	   else:
     		   tmp = OMIM_Phenotype[OMIM_Phenotype['OMIM_id'].str.contains(each_id)]
     		   tmp = tmp['Phenotype'].str.cat(sep = '#')
     		   tmp_lst.append(tmp)
        omim_phenotype.append('||'.join(tmp_lst))
    ds['OMIM_Phenotype'] = pd.Series(omim_phenotype,index = ds.index)
