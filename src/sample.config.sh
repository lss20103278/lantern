#!/bin/bash

######################MTOOLBOX CONFIG FILE###################################################################
##If the default installation of MToolBox was used (install.sh), the user should specify 
##only the MANDATORY parameters. To enable other options, please remove the # before the OPTIONAL parameters
## and assign them the correct value. 
##
##
##to get support, please write to claudia.calabrese23 at gmail.com or to dome.simone at gmail.com or to roberto.preste at gmail.com
##
##
######################SET PATH TO MTOOLBOX EXECUTABLES########################################################
##
##OPTIONAL. If MToolBox default installation (install.sh) was used, samtoolsexe path is $MTOOLBOX_BIN/samtools-{samtools_version}/samtools. Otherwise please specify the FULL PATH to samtools executables.
##
#samtoolsexe=
##
##OPTIONAL. If MToolBox default installation (install.sh) was used, samtools_version is 1.3
##
#samtools_version=1.3
##
##OPTIONAL. If MToolBox default installation (install.sh) was used, musclexe path is $MTOOLBOX_BIN/muscle3.8.31_i86linux64. Otherwise, please specify the FULL PATH to muscle executables.
##
#muscleexe=
##
#OPTIONAL. If MToolBox default installation (install.sh) was used, gsnapexe path is $MTOOLBOX_BIN/gmap/bin/gsnap. Otherwise, please specify the FULL PATH to gsnap executables.
##
#gsnapexe=
##
#####################SET FILE NAMES OF REFERENCE SEQUENCE AND DATABASE TO BE USED FOR THE MAPPING STEP#######
##
#OPTIONAL. If MToolBox default installation (install.sh) was used, fasta_path is $MTOOLBOX_DIR/genome_fasta/. Otherwise, please specify the FULL PATH to fasta and fasta.fai reference sequences to be used in the mapping step.
##
#fasta_path=
##
#OPTIONAL. If MToolBox default installation (install.sh) was used, gsnapdb is $MTOOLBOX_DIR/gmapdb/. Otherwise, please specify the FULL PATH to gsnap/gmap database executables
##
#gsnapdb=
##
##OPTIONAL. If MToolBox default installation (install.sh) was used, mtdb_fasta is chrRSRS.fa (RSRS). To use rCRS sequence as reference, please set mtdb_fasta=chrM.fa Otherwise, please specify the file name of the mitochondrial reference fasta sequence you want to use.
##
mtdb_fasta=chrRSRS.fa
##
##OPTIONAL. If MToolBox default installation (install.sh) was used, hg19_fasta is hg19RSRS.fa (corresponding to the GCh37/hg19 reference multifasta including the RSRS sequence). To use the nuclear multifasta including the rCRS reference sequence, please set hg19_fasta=hg19RCRS.fa. Otherwise, please specify the file name of the nuclear reference sequence you want to use.
##
hg19_fasta=hg19RSRS.fa
##
#OPTIONAL. If MToolBox dafault installation (install.sh) was used, mtdb is chrRSRS. To use the rCRS database, please set mtdb=chrM Otherwise, please specify the name of the mitochondrial gsnap database you want to use.
##
mtdb=chrRSRS
##
#OPTIONAL. If MToolBox default installation (install.sh) was used, the default GSNAP database for nuclear DNA is hg19RSRS. To use the one with rCRS, please set humandb=hg19RCRS, which refers to hg19 nuclear genome + rCRS mitochondrial database. Otherwise, please specify the name of the nuclear (+mitochondrial) gsnap database you want to use.
##
humandb=hg19RSRS
##
######################SET PATH TO INPUT/OUTPUT and PRE/POST PROCESSING PARAMETERS############################
##
##OPTIONAL. Specify the FULL PATH of the input directory. Default is the current working directory
##
input_path=/home/ana005/data/genetic/tmp/mtDNA_20180525/MTDNA
##
##OPTIONAL. Specify the FULL PATH of the output directory. Default is the current working directory
##
output_name=/home/ana005/data/genetic/tmp/mtDNA_20180525/MTDNA/output
##
##OPTIONAL. Specify the list name of files to be analyzed. Default is use all the files with the specified file format extension 
##in the current working directory and skip this option
##
list=sample.list.txt
##
##MANDATORY. Specify the input file format extension. [fasta | bam | sam | fastq ]
##
input_type=fastq
##
##MANDATORY. Specify the mitochondrial reference to be used for the mapping step with mapExome. [RCRS | RSRS; DEFAULT is RSRS]
ref=RSRS
##
##OPTIONAL. Specify if duplicate removal by MarkDuplicates should be set on. [false | true; DEFAULT is false]
UseMarkDuplicates=true
##
##OPTIONAL. Specify if realignment around ins/dels should be set on. [false | true; DEFAULT is false]
UseIndelRealigner=false
##
##OPTIONAL: specify if to exctract only mitochondrial reads from bam file provided as input. [false | true; DEFAULT is false]
MitoExtraction=false
##

