#!/bin/sh


dir=$1
sn=$2
mode=$3


##env:
source `pwd`/../../dbevn.sh


cd $dir

#mode:vqsr
if [ "$mode" = "vqsr" ];then
#SNP
resource_text="--resource $ref_vqsr_hapmap --resource_param hapmap,known=false,training=true,truth=true,prior=15.0 "
resource_text="$resource_text --resource $ref_vqsr_omni --resource_param omni,known=false,training=true,truth=false,prior=12.0 "
resource_text="$resource_text --resource $ref_vqsr_1000G --resource_param 1000G,known=false,training=true,truth=false,prior=10.0 "
resource_text="$resource_text --resource $ref_snp --resource_param dbsnp,known=true,training=false,truth=false,prior=2.0"

annotation_array="QD MQ MQRankSum ReadPosRankSum FS SOR"
for annotation in $annotation_array; do 
	annotate_text="$annotate_text --annotation $annotation"
done

sentieon driver -r $ref_genome --algo VarCal -v $sn.raw.vcf $resource_text $annotate_text --max_gaussians 4 --var_type SNP --plot_file snp.plots.R \
	--tranches_file snp.tranches snp.recal
	
sentieon driver -r $ref_genome --algo ApplyVarCal -v $sn.raw.vcf --sensitivity 99.0 --var_type SNP --tranches_file snp.tranches --recal snp.recal $sn.raw.vqsr.snp.vcf

sentieon plot vqsr -o snp.VQSR.pdf snp.plots.R 

#INDEL
resource_text="--resource $ref_1000M_indel --resource_param mills,known=false,training=true,truth=true,prior=12.0 "
resource_text="$resource_text --resource $ref_snp --resource_param dbsnp,known=true,training=false,truth=false,prior=2.0"

annotate_text=""
annotation_array="QD DP MQRankSum ReadPosRankSum FS SOR"
for annotation in $annotation_array; do 
	annotate_text="$annotate_text --annotation $annotation"
done

sentieon driver -r $ref_genome --algo VarCal -v $sn.raw.vcf $resource_text $annotate_text --max_gaussians 4 --var_type INDEL --plot_file indel.plots.R \
	--tranche 100.0 --tranche 99.9 --tranche 99.0 --tranche 90.0 \
	--tranches_file indel.tranches indel.recal
	
sentieon driver -r $ref_genome --algo ApplyVarCal -v $sn.raw.vcf --sensitivity 99.0 --var_type INDEL --tranches_file indel.tranches --recal indel.recal $sn.raw.vqsr.indel.vcf

sentieon plot vqsr -o indel.VQSR.pdf indel.plots.R 
:<<note
#snp
java -jar $soft_path/GATK/3.7/GenomeAnalysisTK.jar -T VariantRecalibrator -R $ref_genome -L $panel \
 -input $sn.raw.vcf \
 -resource:hapmap,known=false,training=true,truth=true,prior=15.0 /SSD750/PB3/db3/Homo/GATK/hapmap_3.3.hg19.sites.vcf \
 -resource:omni,known=false,training=true,truth=false,prior=12.0 /SSD750/PB3/db3/Homo/GATK/1000G_omni2.5.hg19.sites.vcf \
 -resource:1000G,known=false,training=true,truth=false,prior=10.0 /SSD750/PB3/db3/Homo/GATK/1000G_phase1.snps.high_confidence.hg19.sites.vcf \
 -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /SSD750/PB3/db3/Homo/GATK/dbsnp_138.hg19.vcf \
 -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR --maxGaussians 4 \
 -mode SNP \
 -recalFile snp.recal \
 -tranchesFile snp.tranches \
 -rscriptFile snp.plots.R

java -jar $soft_path/GATK/3.7/GenomeAnalysisTK.jar -T ApplyRecalibration -R $ref_genome -L $panel \
 -input $sn.raw.vcf \
 --ts_filter_level 99.0 -recalFile snp.recal -tranchesFile snp.tranches \
 -mode SNP \
 -o $sn.raw.vqsr.snp.vcf
 
#indel
java -jar $soft_path/GATK/3.7/GenomeAnalysisTK.jar -T VariantRecalibrator -R $ref_genome -L $panel \
 -input $sn.raw.vcf \
 -resource:mills,known=false,training=true,truth=true,prior=12.0 /SSD750/PB3/db3/Homo/GATK/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf \
 -resource:dbsnp,know=true,training=false,truth=false,prior=2.0 /SSD750/PB3/db3/Homo/GATK/dbsnp_138.hg19.vcf \
 -an QD -an DP -an FS -an SOR -an MQRankSum -an ReadPosRankSum --maxGaussians 4 \
 -mode INDEL \
 -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
 -rscriptFile indels.plots.R \
 -recalFile indel.recal \
 -tranchesFile indel.tranches 
 
java -jar $soft_path/GATK/3.7/GenomeAnalysisTK.jar -T ApplyRecalibration -R $ref_genome -L $panel \
 -input $sn.raw.vcf \
 --ts_filter_level 99.0 -recalFile indel.recal -tranchesFile indel.tranches \
 -mode INDEL \
 -o $sn.raw.vqsr.indel.vcf
note
#pass:
awk -F "\t" '$1~/^#/ ||$7=="PASS" {print}' $sn.raw.vqsr.snp.vcf > $sn.raw.snp.pass.vcf
awk -F "\t" '$1~/^#/ ||$7=="PASS" {print}' $sn.raw.vqsr.indel.vcf > $sn.raw.indel.pass.vcf

#combine:
##java -cp $soft_path/GenomeAnalysisTK.jar org.broadinstitute.gatk.tools.CatVariants -R $ref_genome -V $sn.raw.snp.pass.vcf -V $sn.raw.indel.pass.vcf -out $sn.vcf -assumeSorted
#/DATA/sslyu/soft/gatk-4.0.8.1/gatk GatherVcfs -R $ref_genome -I $sn.pass.snp.vcf -I $sn.pass.indel.vcf -O $sn.mul.vcf 
fi


##mode:hard
if [ "$mode" = "hard" ];then

/DATA/sslyu/soft/gatk-4.0.8.1/gatk IndexFeatureFile -F $sn.raw.vcf

/DATA/sslyu/soft/gatk-4.0.8.1/gatk SelectVariants -R $ref_genome -select-type SNP --variant $sn.raw.vcf -O $sn.raw.snp.vcf

/DATA/sslyu/soft/gatk-4.0.8.1/gatk VariantFiltration -R $ref_genome --variant $sn.raw.snp.vcf -O $sn.raw.snp.fil.vcf --filter-name "REJECT" --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"

/DATA/sslyu/soft/gatk-4.0.8.1/gatk SelectVariants -R $ref_genome -select-type INDEL --variant $sn.raw.vcf -O $sn.raw.indel.vcf

/DATA/sslyu/soft/gatk-4.0.8.1/gatk VariantFiltration -R $ref_genome --variant $sn.raw.indel.vcf -O $sn.raw.indel.fil.vcf --filter-name "REJECT" --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"

awk '/^#/ || $7 == "PASS"' $sn.raw.snp.fil.vcf >$sn.raw.snp.pass.vcf
awk '/^#/ || $7 == "PASS"' $sn.raw.indel.fil.vcf >$sn.raw.indel.pass.vcf

#/DATA/sslyu/soft/gatk-4.0.8.1/gatk GatherVcfs -R $ref_genome -I $sn.pass.snp.vcf -I $sn.pass.indel.vcf -O $sn.mul.vcf # snp.vcf and indel.vcf must be in genomic order
#/DATA/sslyu/soft/gatk-4.0.8.1/gatk SortVcf -R $ref_genome -I $sn.pass.snp.vcf -I $sn.pass.indel.vcf -O $sn.mul.vcf 

fi


cd ../
