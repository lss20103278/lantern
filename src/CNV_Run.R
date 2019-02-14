## set R version
#export PATH="/DATA/ypliu/opt/R-3.4.3/bin/:$PATH"

.libPaths("/home/ana005/R/x86_64-pc-linux-gnu-library/3.4")
library(ExomeDepth)
library(GenomicRanges)
load("bam.counts.Rdata")
bam.counts.dafr <- as(bam.counts[, colnames(bam.counts)], 'data.frame')
bam.counts.dafr$chromosome <- gsub(as.character(bam.counts.dafr$space), pattern = 'chr', replacement = '')

#
#read.table("GSG.bed")->GSG.bed
#colnames(GSG.bed)=c("chromosome","start","end")
#colnames(GSG.bed)=c("chromosome","start","end","name")
#

read.table("id")->id
id=as.vector(unlist(id))

#### get the annotation datasets to be used later
data(Conrad.hg19)
exons.hg19.GRanges <- GRanges(seqnames = bed$chromosome,
IRanges(start=bed$start,end=bed$end),
names = bed$name)

### prepare the main matrix of read count data
bam.counts.mat <- as.matrix(bam.counts.dafr[, grep(names(bam.counts.dafr), pattern = '*bam')])
nsamples <- ncol(bam.counts.mat)
### start looping over each sample
for (i in 1:nsamples) {
#### Create the aggregate reference set for this sample
reference_counts = as.matrix(bam.counts.mat[,-i])
colnames = dimnames(bam.counts.mat)[[2]][-i]
dimnames(reference_counts) = list(NULL, colnames) # solve the problem of just two samples
my.choice <- select.reference.set (test.counts = bam.counts.mat[,i],
reference.counts = reference_counts,
bin.length = (bam.counts.dafr$end - bam.counts.dafr$start)/1000,
n.bins.reduced = 10000)

my.reference.selected <- apply(X = bam.counts.mat[, my.choice$reference.choice, drop = FALSE],
MAR = 1,
FUN = sum)

message('Now creating the ExomeDepth object')

all.exons <- new('ExomeDepth',
test = bam.counts.mat[,i],
reference = my.reference.selected,
formula = 'cbind(test, reference) ~ 1')

################ Now call the CNVs
all.exons <- CallCNVs(x = all.exons,
transition.probability = 10^-4,
chromosome = bam.counts.dafr$space,
start = bam.counts.dafr$start,
end = bam.counts.dafr$end,
name = bam.counts.dafr$names)
########################### Now annotate the ExomeDepth object
all.exons <- AnnotateExtra(x = all.exons,
reference.annotation = Conrad.hg19.common.CNVs,
min.overlap = 0.5,
column.name = 'Conrad.hg19')

all.exons <- AnnotateExtra(x = all.exons,
reference.annotation = exons.hg19.GRanges,
min.overlap = 0.0001,
column.name = 'exons.hg19')

#output.file <- paste('Exome_', i, '.csv', sep = '')
output.file <- paste(as.character(id[i]), '_cnv.csv', sep = '')

write.csv(file = output.file, x = cbind(id[i],all.exons@CNV.calls), row.names = FALSE)

}

