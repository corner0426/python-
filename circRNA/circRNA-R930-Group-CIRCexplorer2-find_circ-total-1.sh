#!/bin/bash
#coding:utf-8
#此脚本为R930的circRNA分析主脚本，为山西样本分析完整版
#输入格式：$1是物种，$2是所有样本的绝对路径，每个样本存在于自己单独的目录中
#该脚本是针对R930的circRNA公司使用脚本,使用软件为CIRCexplorer2 + find circ，最终结果为circRNA的count表达矩阵；
#测试脚本,小样本
#bash circRNA-R930-group-CIRCexplorer2-find_circ-total-1.sh human /data1/yaoyh/SXJC_circ /data1/nkm/mayunlong-all-project/project-xiyan/rawdata/lncRNA #注意末尾不能再有反斜杠

:<<!
    #从上到下模块结构为：
    #一、解释说明
    #二、初始化信息和物种信息配置
    #三、样本地址整理并调用副脚本
    #四、副脚本完成监控以及后续常规分析
    #五、A/B流程判断以及差异分析
    #六、靶向预测以及富集分析
    #七、脚本末尾处理
!

#-------
#--------------------
#----------------------------------      一、解释说明
#----------------------------------------------
#----------------------------------
#--------------------
#-------
#主脚本输入格式解释：
#$1是物种，$2是任务运行目录，例如（/data1/biocloud/cloud-task/workflow/20181030_0000002），$3是原始数据存放父文件夹，在该文件夹下，每个样本单独一个文件，存放双端测序fastq
:<<!
#R930下的脚本小样本测试：



!

#使用的软件以及版本：
    #FastQC  v0.11.7
    #trim_galore   version 0.4.4_dev
    #miRdeep2   mirdeep2_0_0_8
    #miRanda    v3.3a
#使用的小脚本：
    

#文件功能备注：
    #pair.txt   样本单双端信息
    #sample.txt  样本地址目录，
    #s_count.txt  完成样本实时计数，按行计数
    #time.txt   任务以及任务时间记录，同时可用于任务进度条计算；
    #task_progress_bar.txt   任务进度条；
    #command.txt  命令记录，用于报错时检查问题所在；
    #Error1.txt  错误日志，偏流程前期的错误监控，检查依据是某一步是否生产满足足够行数的结果，如果不满足，则生成错误日志，方便排错；
    #variate.txt  变量记录；
    #Error2.txt  类似Error1.txt，偏流程后期的监控
    #project_type.txt   后期流程类型，分A和B两类，A是快速结束流程，B是进行详细的流程分析；
    
#其他：
    #---绝对路径软件：python脚本，perl脚本，mirDeep2等
    #服务器级别全局防范高内存占用软件：STAR  CIRCexplorer2  等
    #脚本缩进有三层：副脚本结束后分析、A\B流程、多分组循环；


#-------
#--------------------
#----------------------------------       二、初始化信息和物种信息配置
#----------------------------------------------
#----------------------------------
#--------------------
#-------
#_________________________________________初始化信息______________________________
#建立批量目录和文件：
mkdir -p ${2}   #生成任务运行目录
cd ${2}
mkdir workflow
mkdir star-circ-bam
mkdir ./CIRCexplorer2
mkdir result 
mkdir result/1.Fastqc  #结题质控目录top
mkdir result/1.Fastqc/Raw-fastqc 
mkdir result/1.Fastqc/Clean-fastqc   
mkdir result/2.Alignment   #结题比对目录
mkdir result/2.Alignment/Log   #结题比对目录单个文件
mkdir result/3.circRNA_characteristic
mkdir result/4.Expression  #表达矩阵
touch Error1.txt      #循环样本的报错问题
touch Error2.txt      #样本整合后是否出问题
touch variate.txt   #保存各个变量的值，用于方便错误检查
touch s_count.txt   #设计计数模式，s_count.txt保存副脚本样本完成数，行数代表完成个数
#touch task_progress_bar.txt   #任务进度百分比
touch command.txt
touch workflow/workflow.log

#变量设置
biosoft=/data1/biocloud/biosoft
export biosoft
all_genome_dir=/data1/nkm/biocloud/genome
now_dir=${2}
export  now_dir
task_No=${3}
task_date=`echo ${task_No%%_*}`

#输出日志
    echo circRNA start ${1} ${2} `date +'%Y%m%d:%H:%M'`  ${3} ${4} ${5} ${6} >> ./workflow/workflow.log
    echo circRNA start ${1} ${2} `date +'%Y%m%d:%H:%M'`  ${3} ${4} ${5} ${6} >> ${now_dir}/time.txt
    echo input:     >> variate.txt
    echo task_No   ${task_No}   >> variate.txt
    echo task_date   ${task_date}   >> variate.txt
    echo biosoft  ${biosoft}  >> variate.txt
    echo now_dir   ${now_dir}   >> variate.txt
    echo all_genome_dir   ${all_genome_dir}   >> variate.txt
    
#_________________________________________物种信息配置______________________________
species=${1}  
#人=human
if [ ${1} = human ];
then
genome_star_dir=/data1/biocloud/genome/STAR/hg38_STAR_index
#genome_bowtie2_index=/data/Public/bowtie2-hg38-index/bowtie2_index
genome_fa=/data1/biocloud/genome/human_ensemble_chr/hg38.fa
#genome_split_fa=/data/yaoyh/public_data/hg38_fasta_split
genome_gtf=/data1/biocloud/genome/human_ensemble_chr/gencode.v27.annotation.gtf
genome_dir=/data1/biocloud/genome/human_ensemble_chr
circ_ref_anno=/data/yaoyh/circRNA/CIRCexplorer-anno/hg38_ref_all.txt

species_db_name=org.Hs.eg.db
gene_name_type=ensembl
species_kegg_name=hsa

#export genome_bowtie2_index
#export genome_split_fa

fi

#小鼠=mus,没建好
if [ ${1} = mus ];
then
genome_star_dir=/data1/biocloud/genome/hg38_STAR_index
genome_fa=/data1/biocloud/genome/human_ensemble_chr/hg38.fa
genome_gtf=/data1/biocloud/genome/human_ensemble_chr/gencode.v27.annotation.gtf
genome_dir=/data1/biocloud/genome/human_ensemble_chr
circ_ref_anno=/data/yaoyh/circRNA/CIRCexplorer-anno/mm10_all.txt

fi

#水稻=rice
#if  [ $species = rice ];
#then

#仓鼠=cricetulus
if [ ${1} = cricetulus ];
then
genome_bowtie2_index=/data1/biocloud/genome/bowtie2/CriGri_bowtie2_index/CriGri
#genome_bowtie2_index=/data/Public/bowtie2-hg38-index/bowtie2_index
genome_split_fa=/data1/biocloud/genome/genome/cricetulus-griseus-ncbi
#genome_split_fa=/data/yaoyh/public_data/hg38_fasta_split
genome_gtf=/data1/biocloud/genome/CriGri_STAR_index/Cricetulus_griseus_crigri.CriGri_1.0.94.gtf
#genome_dir=/data1/biocloud/genome/CriGri_STAR_index/

mkdir ./bowtie2_circ_bam
mkdir ./find_circ

export genome_bowtie2_index
export genome_split_fa
export genome_gtf
fi

export genome_star_dir
export genome_fa
export genome_gtf
export genome_dir
export circ_ref_anno

echo  ==========================   >> variate.txt
echo genome_star_dir ${genome_star_dir} >> variate.txt
echo genome_fa ${genome_fa} >> variate.txt
echo genome_gtf ${genome_gtf} >> variate.txt
echo genome_dir ${genome_dir} >> variate.txt
echo circ_ref_anno ${circ_ref_anno} >> variate.txt

echo  ==========================   >> variate.txt
echo species configuration over `date +'%Y%m%d:%H:%M'`   >> ${now_dir}/time.txt

#建库代码
#star建库-human
#hg38_STAR_index]$ STAR --runThreadN 8 --runMode genomeGenerate --genomeDir /data1//biocloud/genome/hg38_STAR_index/ \
#--genomeFastaFiles ../human_ensemble_chr/hg38.fa --sjdbGTFfile /data1/biocloud/genome/human_ensemble_chr/gencode.v27.annotation.gtf \
#--sjdbOverhang 100
#star建库-仓鼠
#STAR --runThreadN 8 --runMode genomeGenerate --genomeDir /data/yaoyh/public_data/CriGri_STAR_index/ \
#--genomeFastaFiles Cricetulus_griseus_crigri.CriGri_1.0.dna.toplevel.fa --sjdbGTFfile Cricetulus_griseus_crigri.CriGri_1.0.94.gtf \
# --sjdbOverhang 100 --limitGenomeGenerateRAM 80920034800
#Bowtie建库-仓鼠
#bowtie2-build -f /data/yaoyh/public_data/CriGri_STAR_index/Cricetulus_griseus_crigri.CriGri_1.0.dna.toplevel.fa CriGri




#-------
#--------------------
#----------------------------------        三、样本地址整理并调用副脚本
#----------------------------------------------
#----------------------------------
#--------------------

#-------
#_________________________________________样本地址整理______________________________

#针对每个样本的双端或单端测序结果单独放在一个目录下的场景：

rm -f ./sample.txt
sleep 3
for i in $3/*
do
echo $i >>  ./sample.txt
done
unset i 
sample_count=`wc -l ./sample.txt| awk '{print $1;}'`

#_____________________________任务数量预估，以及任务进度条建立_________________________
#输入样本数量：    
sample_count=`wc -l ./sample.txt| awk '{print $1;}'`
echo   sample_count  ${sample_count}  >> variate.txt
#总time.txt预估行数，用于估算任务进度；
    #python ${biosoft}/script/deseq2-group.py $5
    #group_line=`wc -l ./group-2.txt | awk '{print $1;}'`
    task_all=$[ ${sample_count} * 7 + 16 ]
    #此步的task_all需要详细的测试和计算
    echo   task_all  ${task_all}  >> variate.txt
    export task_all


#建立任务进度统计函数



#_________________________________________循环调用副脚本______________________________
#对每个样本进行循环调用绝对路径下的miRNA-R730-cloud-bowtie-miRDeep2-analyst-tanxing-single-2.sh脚本
sample_all=`cat ./sample.txt`
i=0
for s in $sample_all
do
i=$[ ${i} + 1 ] 
#bash  ./circRNA-R930-group-CIRCexplorer2-find_circ-single-2.sh  ${1} ${s} ${now_dir} ${i} ${i} &
#Group版，$1是物种，$s是样本路径及名称，$i,$i是数字编号，用于样本索引
bash ./circRNA-R930-Group-CIRCexplorer2-find_circ-single-2.sh ${1} ${s} ${now_dir} ${i} ${i} &
#Cloud版，${4}是物种，$s是样本路径及名称，${task_No}是第二个参数，为任务编号，$i是数字编号，用于样本索引
echo ${s} circRNA-R930-Group-CIRCexplorer2-find_circ-single-2.sh start
##注：第四个参数，即${i}和${task_No}在single脚本中的用法暂时不清楚

done
unset i

echo all associate script start  `date +'%Y%m%d:%H:%M'`   >> ${now_dir}/time.txt


#_______________________
#______________________________________________
#___________________________________________________________________________ 四、副脚本完成监控以及后续常规分析
#______________________________________________
#_______________________

#______________________________________________副脚本完成监控__________________________________


#监测样本完成计数，如果样本完成数达到计数，则开始后续的各种分析

#监测完成计数，如果达到计数，则开始整理、统计和合并信息
while true
do
#s_over_num=`cat ./s_count.txt `
s_over_num=`wc -l ./s_count.txt | awk '{print $1;}'`

echo sample all is ${sample_count} , sampple done is ${s_over_num}

if [ ${s_over_num} -eq ${sample_count} ]
then
cd ${now_dir}
echo all associate script is over and main script start  `date +'%Y%m%d:%H:%M:%S'` >>  ${now_dir}/time.txt


    
    #______________________________________________小目录_________________________________
    #1.统计FASTQC和比对信息
    #2.整合表达矩阵
    #3.对Error1结果进行判断，执行A流程（提前结束）或者B流程（后续分析）

#_____________________________1.统计FASTQC和比对信息_________________________________
cd ${now_dir}/result/1.Fastqc/Raw-fastqc
for i in *.zip
do
unzip ${i}
done
unset i
cp 
perl ${biosoft}/script/tiqu_fastqc_rawdata.pl
mv rawdata_fastqc_result.txt  ../
cd ${now_dir}
cd result/1.Fastqc/Clean-fastqc
for i in *.zip
do
unzip ${i}
done
unset i
perl ${biosoft}/script/tiqu_fastqc_cleandata.pl
mv cleandata_fastqc_result.txt ../
cd ${now_dir}/result/1.Fastqc/
mkdir data
mv Raw-fastqc  ./data/
mv Clean-fastqc  ./data/
cd ${now_dir}
echo 1. fastqc-stat merge `date +'%Y%m%d:%H:%M:%S'` >> ${now_dir}/time.txt


    #_______________________________#2.整合表达矩阵________________________________
## CIRCexplorer2
if [ ${1} = human ] || [ ${1} = mus ];
then
#2 统计比对信息（STAR）
python ${biosoft}/script/circRNA_alignment_statistic.py ${now_dir}/result/2.Alignment/Log ${now_dir}/result/2.Alignment/CIRCexplorer2_ST@AR_alignment_sum.csv

#3 环状RNA特征分析
##3.1环状RNA基因组注释
#mkdir result/3.circRNA_characteristic
mkdir ${now_dir}/result/3.circRNA_characteristic/annotation
##I.确定gtf/gff文件已经转成bed格式(chr	start	end	name/feature	score/gene	strand)
anno_file=/data1/biocloud/genome/genome/annotation/gencode.v27.annotation.bed
##II.执行bed6_CIRCexplorer2.py,制作项目环状RNA的bed文件
python ${biosoft}/script/bed6_CIRCexplorer2.py ${now_dir}/CIRCexplorer2 ${now_dir}/result/3.circRNA_characteristic/annotation/CIRCexplorer2_6item.bed
###III.对gtf_bed和CIRCexplorer2进行intersection
#该步骤耗时，因此移到RPM filter后只注释高表达的环状RNA
#cd ${now_dir}/result/3.circRNA_characteristic/annotation
#nohup python ${biosoft}/script/circ_feature_anno.py CIRCexplorer2_6item.bed ${anno_file} CIRCexplorer2_anno.txt &
cd ${now_dir}

##3.2获取长度及外显子个数
mkdir result/3.circRNA_characteristic/len_exon_num
cd result/3.circRNA_characteristic/len_exon_num

python ${biosoft}/script/circRNA-len_exon_num.py ${now_dir}/CIRCexplorer2 circRNA-len_exon_num.txt

###作图

cd ${now_dir}

#4 提取环状RNA count && RPM矩阵
##count矩阵
#python ${biosoft}/script/CIRCexplorer2_count_extract.py ${now_dir}/CIRCexplorer2 ${now_dir}/result/3.Expression/circRNA_count/circRNA_matrix.csv
rm -rf CIRCexplorer2_sorted_name.txt
cd ${now_dir}/CIRCexplorer2
ls ./ | xargs stat -c "%n" | sort -n  >>  ../CIRCexplorer2_sorted_name.txt
cd ${now_dir}
CIRCexplorer2_name=`cat CIRCexplorer2_sorted_name.txt`

rm -rf CIRCexplorer2_output_files
for i in ${CIRCexplorer2_name}
do
list_csv=""
i_name3=${i#*@}
list_csv=CIRCexplorer2/${i}/${i_name3}_circularRNA_known.txt
echo ${list_csv} >> CIRCexplorer2_output_files
done

python ${biosoft}/script/circRNA-count-CIRCexplorer2.py CIRCexplorer2_output_files result/4.Expression/circRNA_count_matrix.txt

unset list_csv
## RPM矩阵
rm -rf STAR_output_files
for i in ${CIRCexplorer2_name}
do
list_csv=""
i_name3=${i#*@}
list_csv=star-circ-bam/${i}/${i_name3}Log.final.out
echo ${list_csv} >> STAR_output_files
done

python ${biosoft}/script/circRNA-mapped-reads.py STAR_output_files result/4.Expression/total_mapped_reads.csv

cd result/4.Expression/
python ${biosoft}/script/circRNA-RPM-CIRCexplorer2.py circRNA_count_matrix.txt total_mapped_reads.csv circRNA_rpm_matrix.txt
cd ${now_dir}
unset list_csv

fi #CIRCexplorer2 统计结束

## find circ 
if [ ${1} = cricetulus ];
then

#2 统计比对信息（bowtie2）


#3 环状RNA特征分析
##3.1环状RNA基因组注释
#mkdir result/3.circRNA_characteristic
mkdir result/3.circRNA_characteristic/annotation

anno_file=/data1/biocloud/genome/genome/annotation/GCF_CriGri_1.0.bed
python ${biosoft}/script/bed6_CIRI.py ${now_dir}/CIRI_output ${now_dir}/result/3.circRNA_characteristic/annotation/CIRI_6item.bed 
cd ${now_dir}/result/3.circRNA_characteristic/annotation
python ${biosoft}/script/circ_feature_anno.py CIRI_6item.bed ${anno_file} CIRI_anno.txt

cd ${now_dir}
#4 提取环状RNA count矩阵
rm -rf CIRI_sorted_name.txt
cd ${now_dir}/CIRI_output
ls ./ | xargs stat -c "%n" | sort -n  >>  ../CIRI_sorted_name.txt
cd ${now_dir}
CIRI_name=` cat ./CIRI_sorted_name.txt  `

rm -rf CIRI_output_files
for i in ${CIRI_name}
do
list_csv=""
i_name3=${i#*@}
list_csv=CIRI_output/${i}/${i_name3}_CIRI_out.txt
echo ${list_csv} >> CIRI_output_files
#list_name=`echo -e ${list_name}"\t"${i_name3}`
done

#echo $list_name > result/3.Expression/circRNA_count/circRNA_count-matrix.txt
python ${biosoft}/script/circRNA-count-extract.py CIRI_output_files result/4.Expression/circRNA_count_matrix.txt

fi #find circ 统计结束

#_______________________
#______________________________________________
#___________________________________________________________________________ 五、A/B流程判断以及到差异分析
#______________________________________________
    #________________________________A/B流程判断________________________________
    #对Error1结果进行判断，执行A流程（提前结束）或者B流程（后续分析）
doc_line_number=`wc -l ${now_dir}/result/4.Expression/circRNA_rpm_matrix.txt | awk '{print $1;}'`
if [ ${doc_line_number} -le 3 ] || [ ! -f ${now_dir}/result/4.Expression/circRNA_rpm_matrix.txt ];
then
echo circRNA-rpm-matrix >> ${now_dir}/Error1.txt
fi
unset doc_line_number

cd ${now_dir}
Error1_line_number=`wc -l ./Error1.txt | awk '{print $1;}'`
echo ${Error1_line_number} >> ${now_dir}/variate.txt
if [ ${Error1_line_number} -ge 1 ] ;
then
#A流程（提前结束）
echo  task is over
echo A > project_type.txt

fi
#此处fi结束A流程

    #________________________________B流程分析________________________________
    #B流程分析目录：
    # 五：
    #3  对分组信息进行预处理，生成两个文档  ,同时对表达矩阵进行count筛选
    #4  样本热图
    #5  循环做deseq2的差异分析以及差异火山图

if [ ${Error1_line_number} -lt 1 ]
then
#________________________________分组预处理，相关性热图以及count筛选

cd ${now_dir}/result/4.Expression/
#开始用R处理
#4 去除低表达环状RNA（RPM<0.1），并绘制样本correlation热图
cat <<EOF> filter_cor.R 
#!${biosoft}/R-3.5.1/bin/R
#！(1)筛选出最大RPM大于0.1的环状RNA
circ_RPM_matrix <- read.csv('circRNA_rpm_matrix.txt', row.names = 1, header = TRUE)#read.csv自动设置分隔符为逗号，rownames = 1会将第一列设为列名，同时去掉circ_id
circ_filter_rpm_logic<-c()
for (i in 1:nrow(circ_RPM_matrix)) {
  circ_filter_rpm_logic<-c(circ_filter_rpm_logic, max(as.numeric(circ_RPM_matrix[i,1:ncol(circ_RPM_matrix)]))>=0.1)
}
circ_RPM_filter_matirx <- circ_RPM_matrix[circ_filter_rpm_logic,]
#circ_RPM_filter_matirx <- circ_RPM_filter_matirx[!duplicated(circ_RPM_filter_matirx[,1]),]
write.csv(circ_RPM_filter_matirx, 'circ_RPM_filter_matrix.csv')

#！(2)样本相关性热图,基于RPM矩阵
circ_cor <- cor(circ_RPM_filter_matirx)
library("pheatmap")
pdf("circ_RPM_correctoin.pdf")
pheatmap(circ_cor, cluster_rows = FALSE, cluster_cols = FALSE, fontsize_row=10, fontsize_col=10,fontsize=7,display_numbers = T,fontsize_number = 10)
dev.off()

EOF

chmod 755 filter_cor.R
sleep 1
${biosoft}/R-3.5.1/bin/Rscript filter_cor.R
sleep 1

##环状RNA注释

nohup python ${biosoft}/script/circ_feature_anno.py ${now_dir}/result/3.circRNA_characteristic/annotation/CIRCexplorer2_6item.bed circ_RPM_filter_matrix.csv ${anno_file} ${now_dir}/result/3.circRNA_characteristic/annotation/CIRCexplorer2_anno.txt &

#差异分析_pre_circRNA筛选
cat <<EOF> count_pre_matrix.R 
rpm_file <- read.csv('circ_RPM_filter_matrix.csv')
count_file <- read.csv('circRNA_count_matrix.txt')
count_file <- count_file[match(rpm_file[,1], count_file[,1]),] #只保留满足rpm基因名和count基因名相互匹配的数据
write.csv(count_file, 'circRNA_count_filter_matrix.csv', row.names = FALSE)
EOF
chmod 755 count_pre_matrix.R
sleep 1
${biosoft}/R-3.5.1/bin/Rscript count_pre_matrix.R
sleep 1



#B.6从差异分析filter_results中获取候选circRNA,进行miRanda
miRNA_seq=/data1/biocloud/genome/miRNA-genome/human/ref_miRDeep2/mature.hsa.fa
##从差异环状RNA获取bed文件
cd ${now_dir}/result/5.DE/${group_double_name}
python ${biosoft}/script/diff_bed12_circexplorer2.py ${group_double_name}-circ-deseq2_result-filter.csv ${now_dir}/CIRCexplorer2 ${group_double_name}-diff_circ.bed

bedtools getfasta -fi ${genome_fa} -bed ${group_double_name}-diff_circ.bed -s -name -fo ${group_double_name}-diff_circ_seq.fa

miranda ${miRNA_seq} ${group_double_name}-diff_circ_seq.fa -en 20 -strict -out ${group_double_name}-diff_miranda.txt 

python ${biosoft}/script/extract_miranda_target.py ${group_double_name}-diff_miranda.txt ${group_double_name}-diff_target.txt 2

cd ${now_dir}





fi
#此处fi针对Error1的分流，此处B流程结束

break
#针对于if判断单个样本脚本是否执行完（245行），若满足if条件，则执行到该位置，用break结束while循环
fi
#此处fi针对满足样本的数量后的主流程的结束

cd ${now_dir}
unset s_over_num
echo  waiting sample done
sleep 30
done

unset species

#输出结束时间
echo total-end `date +'%Y%m%d:%H:%M:%S'` >> ${now_dir}/time.txt