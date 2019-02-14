#!/bin/sh

:<<Note_script_version
ver:2.1
for:
update:2017/9/19
by:ykl
Note_script_version

:<<HELP
triomode="gvcf | combine"
filtmode="vqsr | hardfilt"

pre-work:
../ped: prepare *.ped for each family
../gvcf: prepare *.g.vcf for each sample
HELP


sn=$1
proc=$2
filtmode=$3

##-->Workflow Control
#usage:select cmds from following
#full functions : "combine vqsr CalculateGenotypePosteriors lowGQ denovo_anotation (merge) vep-mendelscan-annovar sep qcvcf peddy"

path=`pwd`
var=`echo $path | awk -F '/' '{print $NF}'`
if [ $var = 'peddy' ]
then	
cmd="peddy"
else

if [ "$filtmode" != "hard" -a "$filtmode" != "vqsr" ];then
    echo "ERROR: filter method not specified or wrong type('hard' / 'vqsr' only). Exiting ..."
    exit
fi

cmd="combine vqsr CalculateGenotypePosteriors lowGQ denovo_anotation vep-mendelscan-annovar sep qcvcf"
fi	

##cmd:combine
if [ "`echo $cmd|sed "s/ /\n/g"|awk '$1=="combine"{print "on"}'`" = "on" ];then
##env:
source `pwd`/../../dbevn.sh

cp ../../ped/$sn.ped ./
samples=`awk '{print $2}' $sn.ped`

echo "" > $sn.gvcf.list
for i in $samples
do
echo "../../gvcf/$i.g.vcf" >> $sn.gvcf.list
done 
sed -i '/^$/d' $sn.gvcf.list

#gvcf_argument=""  # GATK4.0
#while read -r line
#do
#gvcf_argument=$gvcf_argument" -V $line"
#done < $sn.gvcf.list

gvcf_argument="" # sentieon
while read -r line
do
gvcf_argument=$gvcf_argument" -v $line"
done < $sn.gvcf.list

outdir="triotmp"
[ ! -d $outdir ] && mkdir $outdir
#combine gvcfs:
#java -jar $soft_path/GATK/3.7/GenomeAnalysisTK.jar -T GenotypeGVCFs -R $ref_genome -L $panel --dbsnp $ref_snp \
# -stand_call_conf 10 \
# -o $outdir/$sn.raw.vcf \
# --variant $sn.gvcf.list

for sample in $samples
do
if [ ! -e ../../gvcf/$sample.g.vcf.idx ];then
/DATA/sslyu/soft/gatk-4.0.8.1/gatk IndexFeatureFile -F ../../gvcf/$sample.g.vcf
fi
done

#/DATA/sslyu/soft/gatk-4.0.8.1/gatk CombineGVCFs -R $ref_genome -L $panel --dbsnp $ref_snp -O $outdir/$sn.raw.combine.vcf $gvcf_argument # can't recognise Chinese
#/DATA/sslyu/soft/gatk-4.0.8.1/gatk GenotypeGVCFs -R $ref_genome -L $panel --dbsnp $ref_snp  -stand-call-conf 10  -O $outdir/$sn.raw.vcf  -V $outdir/$sn.raw.combine.vcf

sentieon driver -t $proc -r $ref_genome --interval $panel --algo GVCFtyper --call_conf=10 -d $ref_snp $gvcf_argument $outdir/$sn.raw.vcf
fi 
 
##cmd:vqsr
if [ "`echo $cmd|sed "s/ /\n/g"|awk '$1=="vqsr"{print "on"}'`" = "on" ];then
#sh vfilt.sh $outdir $sn $filtmode
sh /DATA/sslyu/trio_BWA-GATK_3.0/src/vfilt.sh $outdir $sn $filtmode
echo "vqsr `date`" >> finished
fi

##cmd:CalculateGenotypePosteriors
if [ "`echo $cmd|sed "s/ /\n/g"|awk '$1=="CalculateGenotypePosteriors"{print "on"}'`" = "on" ];then
/DATA/sslyu/soft/gatk-4.0.8.1/gatk CalculateGenotypePosteriors -R $ref_genome -V $outdir/$sn.raw.snp.pass.vcf -ped $sn.ped -O $outdir/$sn.snp.gt.vcf # can't recognise Chinese

/DATA/sslyu/soft/gatk-4.0.8.1/gatk CalculateGenotypePosteriors -R $ref_genome -V $outdir/$sn.raw.indel.pass.vcf -ped $sn.ped -O $outdir/$sn.indel.gt.vcf
echo "CalculateGenotypePosteriors `date`" >> finished
fi

##cmd:lowGQ
if [ "`echo $cmd|sed "s/ /\n/g"|awk '$1=="lowGQ"{print "on"}'`" = "on" ];then
/DATA/sslyu/soft/gatk-4.0.8.1/gatk VariantFiltration -R $ref_genome --genotype-filter-expression "DP<20 && GQ<20" --genotype-filter-name "lowGQ" -V $outdir/$sn.snp.gt.vcf -O $outdir/$sn.snp.gt.flt.vcf

/DATA/sslyu/soft/gatk-4.0.8.1/gatk VariantFiltration -R $ref_genome --genotype-filter-expression "DP<20 && GQ<20" --genotype-filter-name "lowGQ" -V $outdir/$sn.indel.gt.vcf -O $outdir/$sn.indel.gt.flt.vcf
echo "lowGQ `date`" >> finished
fi

##cmd:denovo_anotation
if [ "`echo $cmd|sed "s/ /\n/g"|awk '$1=="denovo_anotation"{print "on"}'`" = "on" ];then
java -jar $soft_path/GATK/3.7/GenomeAnalysisTK.jar -T VariantAnnotator -L $panel \
 -R $ref_genome \
 -A PossibleDeNovo \
 -ped $sn.ped \
 -V $outdir/$sn.snp.gt.flt.vcf \
 -o $outdir/$sn.snp.gt.novo.vcf
java -jar $soft_path/GATK/3.7/GenomeAnalysisTK.jar -T VariantAnnotator -L $panel \
 -R $ref_genome \
 -A PossibleDeNovo \
 -ped $sn.ped \
 -V $outdir/$sn.indel.gt.flt.vcf \
 -o $outdir/$sn.indel.gt.novo.vcf
#/DATA/sslyu/soft/gatk-4.0.8.1/gatk VariantAnnotator -L $panel -R $ref_genome -A PossibleDeNovo -ped ../ped/$sn.ped -V $outdir/$sn.snp.gt.flt.vcf -O $outdir/$sn.snp.gt.novo.vcf # Pedigree argument "pedigree" or "founder-id" was specified without a pedigree annotation being requested 
#/DATA/sslyu/soft/gatk-4.0.8.1/gatk VariantAnnotator -L $panel -R $ref_genome -A PossibleDeNovo -V $outdir/$sn.snp.gt.flt.vcf -O $outdir/$sn.snp.gt.novo.vcf #  VariantAnnotator is a BETA tool and is not yet ready for use in production 

#combine
#java -cp $soft_path/GATK/3.7/GenomeAnalysisTK.jar org.broadinstitute.gatk.tools.CatVariants -R $ref_genome -V $outdir/$sn.snp.gt.novo.vcf -V $outdir/$sn.indel.gt.novo.vcf -out $sn.gt.vcf -assumeSorted
#/DATA/sslyu/soft/gatk-4.0.8.1/gatk GatherVcfs -R $ref_genome -I $outdir/$sn.snp.gt.novo.vcf -I $outdir/$sn.indel.gt.novo.vcf -O $sn.gt.vcf  # indel.vcf and snp.vcf must be in genomic order, first record in indel.vcf is not after first record in snp.vcf
/DATA/sslyu/soft/gatk-4.0.8.1/gatk SortVcf -R $ref_genome -I $outdir/$sn.snp.gt.novo.vcf -I $outdir/$sn.indel.gt.novo.vcf -O $sn.gt.vcf
echo "denovo_anotation `date`" >> finished
fi

##cmd:merge
if [ "`echo $cmd|sed "s/ /\n/g"|awk '$1=="merge"{print "on"}'`" = "on" ];then
for sepsn in $samples
do
(grep ^# ../../gvcf/$sepsn.vcf; grep -v ^# ../../gvcf/$sepsn.vcf |sort -k1,1 -k2,2n) |bgzip > $sepsn.sorted.vcf.gz
tabix -p vcf $sepsn.sorted.vcf.gz
done
sample=`awk '{printf $2".sorted.vcf.gz "}' $sn.ped`
bcftools merge ${sample[$@]} > $sn.gt.vcf
echo "merge `date`" >> finished
fi

##cmd:vep-mendelscan-annovar
if [ "`echo $cmd|sed "s/ /\n/g"|awk '$1=="vep-mendelscan-annovar"{print "on"}'`" = "on" ];then
grep -v ^## $sn.gt.vcf>$sn.tmp.vcf
grep ^## $sn.gt.vcf>$sn.head.vcf

#./sort_sample.sh $sn
/DATA/sslyu/trio_BWA-GATK_3.0/src/sort_sample.py $sn
if [ -e $sn.sort.tmp.vcf ]
then
cat $sn.head.vcf $sn.sort.tmp.vcf>$sn.gt.vcf
else
cat $sn.head.vcf $sn.tmp.vcf>$sn.gt.vcf
fi

#sh vep-mendelscan.sh $sn ../../ped/$sn.mendel.ped
sh /DATA/sslyu/trio_BWA-GATK_3.0/src/vep-mendelscan.sh $sn ../../ped/$sn.mendel.ped
echo "vep-mendelscan-annovar `date`" >> finished
fi

##cmd:sep
if [ "`echo $cmd|sed "s/ /\n/g"|awk '$1=="sep"{print "on"}'`" = "on" ];then
outdir='sep'
[ ! -e $outdir ] && mkdir $outdir
for sepsn in $samples
do
#java -jar $soft_path/GenomeAnalysisTK.jar -T SelectVariants -R $ref_genome -L $panel -V $sn.gt.vcf -sn $sepsn -o $outdir/$sepsn.vcf

vcf-subset -c $sepsn $sn.gt.vcf > $outdir/$sepsn.vcf
done
echo "sep `date`" >> finished
fi

##cmd:qcvcf
if [ "`echo $cmd|sed "s/ /\n/g"|awk '$1=="qcvcf"{print "on"}'`" = "on" ];then
echo -e "sample\tindel_genomic\tsnp_genomic\tts/tv_genomic\tsnp_exonic\tts/tv_exonic" > triotmp/$sn.qcvcf
[ ! -d tmprun ] && mkdir tmprun
for sepsn in $samples
do
printf "$sepsn\t" >> triotmp/$sn.qcvcf
/DATA/sslyu/soft/annovar/convert2annovar.pl -format vcf4 sep/$sepsn.vcf -outfile tmprun/$sepsn.avinput
/DATA/sslyu/soft/annovar/table_annovar.pl tmprun/$sepsn.avinput $ref_annotation -buildver hg19 -out tmprun/$sepsn -protocol refGene -operation g -remove

awk -F "\t" 'BEGIN{indel=0}  $4=="-"||$5=="-"{indel++} $4!="-"&&$5!="-"{print $4$5} END{print indel"\t"NR-1-indel"\t"}' tmprun/$sepsn.hg19_multianno.txt | awk 'BEGIN{ts=0} $1=="AG"||$1=="GA"||$1=="CT"||$1=="TC"{ts++} END{printf $0;tv=NR-2-ts;printf (ts+0.0)/tv"\t"}' >> triotmp/$sn.qcvcf

awk -F "\t" '$6=="exonic"||$6=="exonic;splicing"{print $4$5}' tmprun/$sepsn.hg19_multianno.txt | awk '$1!~/-/{print}'| awk 'BEGIN{ts=0} $1=="AG"||$1=="GA"||$1=="CT"||$1=="TC"{ts++} END{tv=NR-1-ts;printf ts+tv"\t"(ts+0.0)/tv"\n"}' >> triotmp/$sn.qcvcf
done
rm -rf tmprun
echo "qcvcf `date`" >> finished
fi

##cmd:peddy
if [ "`echo $cmd|sed "s/ /\n/g"|awk '$1=="peddy"{print "on"}'`" = "on" ];then
##env:
source `pwd`/../dbevn.sh

if [ ! -e $sn.ped ]; then
cat ../ped/*.mendel.ped |awk '{OFS="\t"; print "0",$2,$3,$4,$5,$6}' |sort |uniq > $sn.ped
fi
samples=`awk '{print $2}' $sn.ped`

echo "" > $sn.gvcf.list
for i in $samples
do
echo "../gvcf/$i.g.vcf" >> $sn.gvcf.list
done
sed -i '/^$/d' $sn.gvcf.list

gvcf_argument="" # sentieon
while read -r line
do
gvcf_argument=$gvcf_argument" -v $line"
done < $sn.gvcf.list

outdir="pedtrio"
[ ! -d $outdir ] && mkdir $outdir
#combine gvcfs:

for i in $samples
do
if [ ! -e ../gvcf/$i.g.vcf.idx ];then
/DATA/sslyu/soft/gatk-4.0.8.1/gatk IndexFeatureFile -F ../gvcf/$i.g.vcf
fi
done

sentieon driver -t $proc -r $ref_genome --interval $panel --algo GVCFtyper --call_conf=10 -d $ref_snp $gvcf_argument $outdir/$sn.vcf

cd $outdir
vcf-sort $sn.vcf>$sn.raw.sort.vcf
sed 's/^chr//g' $sn.raw.sort.vcf>$sn.sort.vcf
bgzip -c $sn.sort.vcf>$sn.sort.vcf.gz
tabix -p vcf $sn.sort.vcf.gz
[ ! -e results ] && mkdir -p results
cd results
/home/ana005/anaconda2/bin/python -m peddy -p 4 --plot --prefix $sn ../$sn.sort.vcf.gz ../../$sn.ped
fi

