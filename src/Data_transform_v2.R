library(readxl)
library(tidyr)
library(dplyr)
args = commandArgs(trailingOnly=TRUE)


process <- function(input_f){
                #a <- read.table(input_f, header = TRUE)
                gtinfo <- dplyr::select(input_f, starts_with("GTinfo_"))
                gt <- dplyr::select(input_f, starts_with("GT_"))
                depth <- dplyr::select(input_f, starts_with("Depth_"))
                qual <- dplyr::select(input_f, starts_with("Qual_"))

                #extract column name
                pre_name_gtinfo <- colnames(gtinfo)
                sample_name_gtinfo <- gsub("GTinfo_","", pre_name_gtinfo)

                pre_name_gt <- colnames(gt)
                sample_name_gt <- gsub("GT_","",pre_name_gt)

                pre_name_depth <- colnames(depth)
                sample_name_depth <- gsub("Depth_","",pre_name_depth)

                pre_name_qual <- colnames(qual)
                sample_name_qual <- gsub("Qual_","",pre_name_qual)

                #prefix name to selected column value
                for (i in c(1:length(pre_name_gtinfo))) {
                  input_f[pre_name_gtinfo[i]] <- paste(sample_name_gtinfo[i],":",as.vector(gtinfo[[i]]), sep= "")
                  input_f[pre_name_gt[i]] <- paste(sample_name_gt[i],":",as.vector(gt[[i]]), sep= "")
                  input_f[pre_name_depth[i]] <- paste(sample_name_depth[i],":",as.vector(depth[[i]]), sep= "")
                  input_f[pre_name_qual[i]] <- paste(sample_name_qual[i],":",as.vector(qual[[i]]), sep= "")
                }

                New_a1 <- tidyr::unite(input_f,"GTinfo",pre_name_gtinfo, sep = ";") 
                New_a2 <- tidyr::unite(New_a1,"GT",pre_name_gt, sep = ";")
                New_a3 <- tidyr::unite(New_a2,"Depth",pre_name_depth, sep = ";")
                New_a4 <- tidyr::unite(New_a3,"Qual",pre_name_qual, sep = ";")
                #x_csv <- gsub(".tsv",".csv",x)
                #write CSV in R 
                write.table(New_a4, file=stdout(), sep= "\t",row.names = FALSE) #paste("Combined_",sample_name_gtinfo[1],".csv",sep="")
                #output_filename <- gsub(".txt","",i)
                #if (!file.exists(paste("Combined_",sample_name_gtinfo[1],".csv",sep=""))) {
                #  stop("Plz check your input values")
                #} else {
                #cat(paste("Combined_",sample_name_gtinfo[1],".csv",sep=""," is done. \n"))
                #}
                }

# test if there is at least one argument: if not, return an error
if (length(args)==0 ) {
  stop("At least one argument or flag must be supplied (input file).n", call.=FALSE)
} else if (length(args) >=1) {
 
  if (length(args[args %in% "--stdin"]) == 0) {
     
    for (i in args) {
      tab <- read.delim(i,header=TRUE, sep="\t")
       process(tab)
    }
    } else if (length(args[args %in% "--stdin"]) == 1 ) {
      f <- file("stdin")
      tab <- read.delim(f,header=TRUE, sep="\t")
  # default output file
 
  process(tab)
  
}}





