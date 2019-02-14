#bin/bash
export PATH="/DATA/ypliu/opt/Anaconda2/bin:$PATH"
export  LD_LIBRARY_PATH="/DATA/ypliu/opt/Anaconda2/lib"

work_path="/home/ana005/data/genetic/tmp/mtDNA_20180525/MTDNA"
fastq_path="/home/ana005/data/genetic/tmp/mtDNA_20180525"
chil='18N0444'
#father= ###父亲样本名，和fastq命名一致
mother='18N0446'



cd $work_path
ln -s $fastq_path/*/1_data/*fastq .

for i in *.fastq
do
seqtk sample -s100 $i 200000 >$i.seq  # randomly select 200000 reads
done

####重命名，fastq名必须为 样本名.R1.fastq 和样本名.R2.fastq
rm *.fastq
rename '.seq' '' *
rename "_" "." *

####使用MToolBox软件分析
export PATH="/DATA/ypliu/opt/MToolBox-v.1.0/MToolBox/:$PATH"

MToolBox.sh -i sample.config.sh   ####congfig文件记得修改input_path、output_name、list三个参数，list和config都放在工作路径下,示例参考/DATA/ypliu/Analysis/analysis_mtDNA/20180330/zhangyuyang/zhangyuyang_config.sh和/DATA/ypliu/Analysis/analysis_mtDNA/20180330/zhangyuyang/zhangyuyang_list.txt 



#####注释

cd $work_path/output
cut -f 2,17-29  OUT_$chil/$chil.*.annotation.csv >annotation.temp1

awk 'BEGIN{FS=OFS="\t"}NR==FNR{a[$1]=1;b[$1]=$0;next}{if (a[$1]){print $1,$0,b[$1]}else {print $1,$0,""}}' annotation.temp1 prioritized_variants.txt |awk 'BEGIN{FS=OFS="\t"} gsub(/.$/,"",$1)' >temp1

grep -v "#" VCF_file.vcf |awk 'BEGIN{FS=OFS="\t"}{print $2,$4"-"$5}' |sed 1i'Variant Allel\tRef-Alt'>temp2 

awk 'BEGIN{FS=OFS="\t"}NR==FNR{a[$1]=1;b[$1]=$2;next}{if (a[$1]){print b[$1],$0}}' temp2 temp1 >temp3 

awk 'BEGIN{FS=OFS="\t"}{print $2,$0}' temp3|awk 'BEGIN{FS=OFS="\t"}NR==FNR{a[$1]=1;b[$1]=$2;next}{if (a[$1]){print "Yes",$0}else {print "No",$0}}' /DATA/ypliu/Analysis/analysis_mtDNA/20171027/output/polymorphic_sites -  >temp4

awk -F '\t' '{print $1"\t",$1,$2,$3,$4,$5,$6}' OUT_$chil/$chil-table.txt |sed 1's/Position\t/Variant Allel\tChild /' >temp5

#awk -F '\t' '{print $1"\t",$1,$2,$3,$4,$5,$6}' OUT_$father/$father-table.txt |sed 1's/Position\t/Variant Allel\tfather /' >temp6  ##若没有父亲则注释掉该行

awk -F '\t' '{print $1"\t",$1,$2,$3,$4,$5,$6}' OUT_$mother/$mother-table.txt |sed 1's/Position\t/Variant Allel\tmother /' >temp7

#awk 'BEGIN{FS=OFS="\t"}NR==FNR{a[$1]=1;b[$1]=$2;next}{if (a[$1]){print $0,b[$1]}}' temp6 temp5 |awk 'BEGIN{FS=OFS="\t"}NR==FNR{a[$1]=1;b[$1]=$2;next}{if (a[$1]){print $0,b[$1]}}' temp7 - >temp8   ##若没有父亲样本，则注释掉该行，选择下一行
awk 'BEGIN{FS=OFS="\t"}NR==FNR{a[$1]=1;b[$1]=$2;next}{if (a[$1]){print $0,b[$1]}}' temp7 temp5 >temp8

awk 'BEGIN{FS=OFS="\t"}{print $2,$0}' temp4 |awk 'BEGIN{FS=OFS="\t"}NR==FNR{a[$1]=1;b[$1]=$2FS$3FS$4;next}{if (a[$1]){print b[$1],$0}}' temp8 - >$chil.final

####删除不必要浪费空间的东西
rm temp*
for i in OUT_*; do cd $work_path/output/$i; rm *.fastq *.sam *.pileup; done

#####最后的$chil.final需要进行列的删除调整，参考/DATA/ypliu/Analysis/analysis_mtDNA/17N0895.final.xlsx文件进行调整，表头为No的那列表头改为Polymorphic sites








