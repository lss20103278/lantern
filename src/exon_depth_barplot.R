args = commandArgs(trailingOnly=TRUE)

process <- function(input_f,gene,sample){
	x_name = paste("exon_",seq(nrow(input_f)),sep="")
	max = max(input_f$V2)+50
	pdf(paste(sample,"_",gene,"_exon.depth.pdf",sep=""))
	p = barplot(input_f$V2, space=2, col="red", main=paste(gene,"exon region count of",sample), ylim=c(0,max))
	text(p, input_f$V2, labels = x_name, pos=3)
	dev.off()
}

if (length(args)!=3){
	stop("the number of arguments is wrong, first is the depth file, second is the gene name, third is the sample name", call.=FALSE)
} else {
	tab = read.table(args[1])
	gene = args[2]
	sample = args[3]
	process(tab,gene,sample)
}


