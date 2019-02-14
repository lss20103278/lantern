#/DATA/ypliu/opt/R-3.4.3/bin/R --slave --vanilla < Run.R
export PATH="/DATA/ypliu/opt/R-3.4.3/bin/:$PATH"
.libPaths("/home/ana005/R/x86_64-pc-linux-gnu-library/3.4")
library(ExomeDepth)
#read.table("/home/ana005/data/annovar/bed/IDT/IDT_gene.bed")->bed
#read.table("/home/ana005/data/annovar/bed/wes/wes_gene.bed")->bed
#read.table("/home/ana005/data/annovar/bed/Agilent/S06588914/S06588914_Regions_cnv.bed")->bed
colnames(bed)=c("chromosome","start","end","name")
read.table("list")->bamfile
unlist(bamfile)->bamfile
bamfile=as.character(bamfile)
bamfile <- as.vector(bamfile)
bam.counts <- getBamCounts(bed.frame=bed, bam.files=bamfile, referenceFasta="/SSD750/PB1/db1/Homo/refseq/hg19.fa")
save.image("bam.counts.Rdata")
source("/DATA/sslyu/trio_BWA-GATK_3.0/src/CNV_Run.R")
q()
