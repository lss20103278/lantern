#!/bin/sh


soft_path="/SSD750/PB1/soft1"
DB="/SSD750/PB1/db1/Homo"
ref_genome="$DB/refseq/hg19.fa"
bwa_index="$DB/bwa_index/hg19.fa"
ref_1000G_indel="$DB/GATK/1000G_phase1.indels.hg19.sites.vcf"
ref_1000M_indel="$DB/GATK/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"
ref_annotation="/SSD750/PB3/db3/Homo/annovar"
ref_snp="$DB/GATK/dbsnp_138.hg19.vcf"
ref_vqsr_hapmap="/SSD750/PB3/db3/Homo/GATK/hapmap_3.3.hg19.sites.vcf"
ref_vqsr_omni="/SSD750/PB3/db3/Homo/GATK/1000G_omni2.5.hg19.sites.vcf"
ref_vqsr_1000G="/SSD750/PB3/db3/Homo/GATK/1000G_phase1.snps.high_confidence.hg19.sites.vcf"

panel="$DB/refseq/hg19_SeqCap_EZ_Exome_v3_capture_plus20bp.bed"         #WXS




sn=$1


java -jar /SSD750/PB1/soft1/GenomeAnalysisTK.jar -T VariantRecalibrator -R $ref_genome -L $panel \
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

java -jar /SSD750/PB1/soft1/GenomeAnalysisTK.jar -T ApplyRecalibration -R $ref_genome -L $panel \
 -input $sn.raw.vcf \
 --ts_filter_level 99.0 -recalFile snp.recal -tranchesFile snp.tranches \
 -mode SNP \
 -o $sn.raw.vqsr.snp.vcf
 


java -jar /SSD750/PB1/soft1/GenomeAnalysisTK.jar -T VariantRecalibrator -R $ref_genome -L $panel \
 -input $sn.raw.vcf \
 -resource:mills,known=false,training=true,truth=true,prior=12.0 /SSD750/PB3/db3/Homo/GATK/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf \
 -resource:dbsnp,know=true,training=false,truth=false,prior=2.0 /SSD750/PB3/db3/Homo/GATK/dbsnp_138.hg19.vcf \
 -an QD -an DP -an FS -an SOR -an MQRankSum -an ReadPosRankSum --maxGaussians 4 \
 -mode INDEL \
 -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
 -rscriptFile indels.plots.R \
 -recalFile indel.recal \
 -tranchesFile indel.tranches 
 
java -jar /SSD750/PB1/soft1/GenomeAnalysisTK.jar -T ApplyRecalibration -R $ref_genome -L $panel \
 -input $sn.raw.vcf \
 --ts_filter_level 99.0 -recalFile indel.recal -tranchesFile indel.tranches \
 -mode INDEL \
 -o $sn.raw.vqsr.indel.vcf




#filt pass:
awk -F "\t" '$1~/^#/ ||$7=="PASS" {print}' $sn.raw.vqsr.snp.vcf > $sn.raw.vqsr.snp.pass.vcf
awk -F "\t" '$1~/^#/ ||$7=="PASS" {print}' $sn.raw.vqsr.indel.vcf > $sn.raw.vqsr.indel.pass.vcf



#combine:
java -cp $soft_path/GenomeAnalysisTK.jar org.broadinstitute.gatk.tools.CatVariants -R $ref_genome -V $sn.raw.vqsr.snp.pass.vcf -V $sn.raw.vqsr.indel.pass.vcf -out $sn.vcf -assumeSorted





