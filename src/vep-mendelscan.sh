#!/bin/sh

sn=$1
ped=$2

##env:
source `pwd`/../../dbevn.sh


outdir="segtrio"
[ ! -d $outdir ] && mkdir $outdir


#seperate multiple alles
bcftools norm -m - $sn.gt.vcf | bcftools norm -f $ref_genome >$outdir/$sn.vcf

#annotation:vep
variant_effect_predictor.pl -i $outdir/$sn.vcf  -o $outdir/$sn.vep --cache --offline --dir /SSD750/PB3/db3/Homo/vep/ --sift b --polyphen b --symbol --canonical --force

#annotation:mendelscan
java -jar $soft_path/MendelScan/MendelScan.v1.2.1.fix.jar score $outdir/$sn.vcf --vep-file $outdir/$sn.vep --ped-file $ped --output-file $outdir/$sn.mendel.tsv --output-vcf $outdir/$sn.mendel.vcf --inheritance recessive

perl /DATA/ypliu/opt/mskcc-vcf2maf-1b16b35/vcf2maf.pl --input-vcf $outdir/$sn.vcf --output-maf $outdir/$sn.maf --vep-path /DATA/ypliu/opt/mskcc-vcf2maf-1b16b35/vep --vep-data /DATA/ypliu/opt/mskcc-vcf2maf-1b16b35/.vep --ref-fasta /DATA/sslyu/refGene/hg19.fa --filter-vcf /DATA/ypliu/opt/mskcc-vcf2maf-1b16b35/.vep/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz
sed -i '1d' $outdir/$sn.maf

/DATA/sslyu/soft/annovar/convert2annovar.pl -format vcf4 $outdir/$sn.mendel.vcf -outfile $outdir/$sn.link --includeinfo # annovar 20180416
#sed 's/^chr//' $outdir/$sn.vcf|gzip -c ->$outdir/$sn.vcf.gz
gzip -c $outdir/$sn.mendel.vcf>$outdir/$sn.mendel.vcf.gz
score.sh $outdir/$sn.mendel.vcf.gz $outdir/$sn.CADD.gz
bgzip -df $outdir/$sn.CADD.gz
sed -i '1d' $outdir/$sn.CADD

#annotation:annovar
#./mendel_to_annovar.sh $sn
/DATA/sslyu/trio_BWA-GATK_3.0/src/mendel_to_annovar.py $sn

#/DATA/sslyu/soft/annovar/table_annovar.pl $outdir/$sn.avinput $ref_annotation -buildver hg19 -out $outdir/$sn.ann -protocol refGene,avsnp147,1000g2015aug_all,1000g2015aug_eas,gnomad_genome_eas,gnomad_exomes_hom,gnomad_exomes_hom_all,gnomad_exomes_hemi,exac03_eas,esp6500_all,dbscsnv11,dbnsfp33a,revel,intervar_20180118,clinvar_20180603,ensembl,hgmd_2018_spring,bed,bed,bed,bed -operation g,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,r,r,r,r -bedfile hg19_Dead_Zone.txt,hg19_Problem_High.txt,hg19_Problem_Low.txt,hg19_HI_Predictions_Version3.bed -arg ',,,,,,,,,,,,,,,,,-colsWanted 4,-colsWanted 4,-colsWanted 4,-colsWanted 4' -remove -nastring "NA" -otherinfo # annovar 20180416 avsnp147->avsnp150 intervar_20170202->intervar_20180118 clinvar20170905->clinvar_20180603
/DATA/sslyu/soft/annovar/table_annovar.pl $outdir/$sn.avinput $ref_annotation -buildver hg19 -out $outdir/$sn.ann -protocol refGene,avsnp147,1000g2015aug_all,1000g2015aug_eas,gnomad_genome_eas,gnomad_exomes_hom,gnomad_exomes_hom_all,gnomad_exomes_hemi,exac03_eas,esp6500_all,dbscsnv11,dbnsfp33a,revel,intervar_20180118,clinvar_20181225,ensembl,hgmd_2018_spring,bed,bed,bed,bed -operation g,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,r,r,r,r -bedfile hg19_Dead_Zone.txt,hg19_Problem_High.txt,hg19_Problem_Low.txt,hg19_HI_Predictions_Version3.bed -arg ',,,,,,,,,,,,,,,,,-colsWanted 4,-colsWanted 4,-colsWanted 4,-colsWanted 4' -remove -nastring "NA" -otherinfo # annovar 20180416 avsnp147->avsnp150 intervar_20170202->intervar_20180118 clinvar_20180603->clinvar_20181225
sed -i "1s/bed\tbed2\tbed3\tbed4/Dead_Zone\tProblem_High\tProblem_Low\tHI_Predictions/" $outdir/$sn.ann\.hg19_multianno.txt

otherinfo=`cat $outdir/$sn.avinput.info`
sed -i "s/Otherinfo/`echo $otherinfo|sed "s/ /\t/g"`/" $outdir/$sn.ann.hg19_multianno.txt


