#!/bin/bash
#此脚本为circRNA分析的副脚本 CIRCexplorer2为主，跑完后即记数，find_circ也跑，但不是重点
#输入设计：$1是物种，$2是样本信息，$3是脚本工作初始执行路径
#此脚本只针对circRNA的分析

#从上到下模块结构为：
#一、初始化信息
#二、物种信息
#三、清洗和质控
#四、联配、组装和定量
#五、脚本末尾处理

#-------
#--------------------
#----------------------------------    一
#----------------------------------------------
#----------------------------------
#--------------------
#-------

#一、初始化信息-----------------------------------

#剩余内存监控函数：

function R930_mem_check()
{
#两个输入，第一个剩余内存大小，第二个是waiting的对象名称,每隔10s判断一次；现有剩余内存大于$1,即可往下执行；部分复杂的判断调用比较少，就不封装函数了
while true
do
free_mem=` free -g | grep 'buffers/cache' | grep -v 'grep' | awk '{print $4;}' `
if [ ${free_mem} -ge ${1} ]
then
break
fi
sleep 10
echo waiting 10 s  ${s_name}  $2
done
}

#建立任务进度函数
function task_progress_bar()
{
    time_count=`wc -l time.txt | awk '{print $1;}' `
    task_over_now=$[ ${time_count} * 100 / ${task_all} ] 
    echo ${task_over_now} >> task_progress_bar.txt
}

#监控内存大小设置
fastqc_mem=5
trim_galore_mem=5
star_mem=70
star_have=3
bowtie2_mem=50
bowtie2_have=3
#rsem_thread=5
#stringtie_thread=6

#-------
#--------------------
#----------------------------------    二
#----------------------------------------------
#----------------------------------
#--------------------
#-------

#二、物种信息-----------------------------------------------------
#暂时省略，可用主脚本的export输出的物种相关变量

#-------
#--------------------
#----------------------------------    三
#----------------------------------------------
#----------------------------------
#--------------------
#-------

#三、清洗和质控------------------------------------------------------

#是否是gz压缩文件，以及解压缩；
for i in ${2}/*gz
do
if [ -f ${i} ]
then
gzip -d ${i}
fi
done
unset i

#每个样本建立子文件夹,以及名称简化
cd $3
s=$2
s_name=`echo ${s##*/}`
echo ${s_name}
cd ${s}
mkdir ./trim
cd $3

#数据是fastq或fq结尾判断
end=0
for i in  ${s}/*fq
do
if [ -f $i ];then
end="fq"
fi
done
unset i
for i in  ${s}/*fastq
do
if [ -f $i ];then
end="fastq"
fi
done
unset i

#pair判断，双端=2，单端=1
pair=0
for i in  ${s}/*fq;do
if [ -f $i ];then
pair=$[ ${pair} + 1 ]
fi
done
unset i
for i in  ${s}/*fastq;do
if [ -f $i ];then
pair=$[ ${pair} + 1 ]
fi
done
unset i
echo ${s_name} is ${pair}
echo ${s_name}  ${pair} ${end} >> ${3}/pair.txt

#记录单个样本开始时间
echo ${s_name}-start `date +'%Y%m%d:%H:%M:%S'`  >> ${3}/time.txt

#剩余可用内存大于fastqc_mem时通过，并执行下一步进行FASTQC质控，每隔10s执行一次
R930_mem_check 50 Fastqc
#raw_data fastqc执行
for i in  ${s}/*${end};do
if [ -f $i ];then
fastqc   -o  ./result/1.Fastqc/Raw-fastqc $i  
fi
done
echo ${s_name}-fastqc-over `date +'%Y%m%d:%H:%M:%S'`  >> ${3}/time.txt
task_progress_bar
cd ${3}

#trim_galore单双端判断
if [ ${pair} -eq 2 ];then
trim_galore_pair="--paired"
fi
if [ ${pair} -eq 1 ];then
trim_galore_pair=""
fi

#剩余可用内存大于trim_galore_mem时通过，并执行trim_galore，每隔10s执行一次
R930_mem_check 50 trim_galore
#执行trim_galore
trim_galore ${trim_galore_pair} -q 20 --length 14 --fastqc_args "-o ./result/1.Fastqc/Clean-fastqc" -o ${s}/trim  ${s}/*.${end} 
echo ${s_name}-trim_galore-over `date +'%Y%m%d:%H:%M:%S'`  >> ${3}/time.txt
task_progress_bar
#错误check，如果文件行数达不到要求或者文件不存在，则想错误日志文件存入记录
doc_line_number=`wc -l ${s}/trim/*.fq | head -n 1 | awk '{print $1;}'`
#echo ${doc_line_number}
if [ ${doc_line_number} -le 3 ] || [ ! -f ${s}/trim/*1.fq ];
then
echo ${s_name} sample is wrong: the error is in trim_galore step, please check the fastq input  >> ${3}/Error1.txt
fi
unset doc_line_number

#gzip压缩文件
#gzip  ${s}/*.${end}

cd ${3}
  
#-------
#--------------------
#----------------------------------      四
#----------------------------------------------
#----------------------------------
#--------------------
#-------

#四、联配、组装和定量


## CIRCexplorer2
if [ ${1} = human ] || [ ${1} = mus ];
then
#测试用的if


#设置睡眠，睡眠时间随机确定，取值范围在1—101秒之内;此处的意义为避免多个STAR程序（单个内存>30G）同时运行造成服务器内存耗尽。
time_sleep=$[$RANDOM%100+1]
echo $s_name first_sleep ${time_sleep} >> ${3}/STAR_waiting.txt
echo waiting $s_name first_sleep ${time_sleep} s STAR
sleep ${time_sleep}

#剩余可用内存大于star_mem,并且已经运行的STAR程序数量小于star_have时，进入STAR运行程序；
while true
do
free_7=` free -g | grep 'buffers/cache' | grep -v 'grep' | awk '{print $4;}' `
if [ ${free_7} -ge ${star_mem} ]
then
star_count=`ps -ef | grep STAR | grep -v "grep" | awk '{print $2;}' | wc -l`
if [ ${star_count} -lt ${star_have} ]
then
break 
fi
unset star_count
fi
sleep 10
done

#第二次睡眠，睡眠时间随机确定，取值范围在1—101秒之内;
time_sleep2=$[$RANDOM%200+50]
echo $s_name second_sleep ${time_sleep2} >>  ${3}/STAR_waiting.txt
echo waiting $s_name second_sleep  ${time_sleep2} s STAR
sleep ${time_sleep2}
while true
do
free_7=` free -g | grep 'buffers/cache' | grep -v 'grep' | awk '{print $4;}' `
if [ ${free_7} -ge ${star_mem} ]
then
star_count=`ps -ef | grep STAR | grep -v "grep" | awk '{print $2;}' | wc -l`
if [ ${star_count} -lt ${star_have} ]
then
echo ${s_name}  star-start
echo ${s_name}  star-start  `date +'%Y%m%d:%H:%M:%S'`  >>  ${3}/time.txt 
break 
fi
unset star_count
fi
sleep ${time_sleep2}
done


#记录单个样本开始时间
echo ${s_name}-CIRCexplorer2 start `date +'%Y%m%d:%H:%M:%S'`  >> ${3}/time.txt


#STAR运行,以及单双端判断：
if [ ${pair} -eq 2 ];then
mkdir ./star-circ-bam/${5}@${s_name}
STAR --runThreadN 8 --chimSegmentMin 10 --genomeDir ${genome_star_dir} --readFilesIn ${s}/trim/*val_1.fq ${s}/trim/*val_2.fq --outFileNamePrefix  ./star-circ-bam/${5}@${s_name}/${s_name} > star-circ-bam/${5}@${s_name}/${s_name}".star.log"
echo $s_name  STAR-end  `date +'%Y%m%d:%H:%M:%S'`  >> ${3}/time.txt
fi

if [ ${pair} -eq 1 ];then
mkdir ./star-circ-bam/${5}@${s_name}  
STAR --runThreadN 8 --chimSegmentMin 10 --genomeDir ${genome_star_dir} --readFilesIn ${s}/trim/*.fq --outFileNamePrefix  ./star-circ-bam/${5}@${s_name}/${s_name} > star-circ-bam/${5}@${s_name}/${s_name}".star.log"
echo $s_name  STAR-end  `date +'%Y%m%d:%H:%M:%S'`  >> ${3}/time.txt
fi

#错误check，如果没有生成联配主要文件，则记录到错误日志
doc_line_number=`wc -l ./star-circ-bam/${5}@${s_name}/${s_name}Chimeric.out.junction | awk '{print $1;}'`
if [ ${doc_line_number} -le 3 ];
then
echo ${s_name} sample STAR is wrong: the error is in STAR step, please check the fastq input  >> ${3}/Error1.txt
fi
unset doc_line_number

mkdir ./result/2.Alignment/Log/${s_name}
cp ./star-circ-bam/${5}@${s_name}/*Log.final.out  ./result/2.Alignment/Log/${s_name}/ 

#CIRCexplorer2 parse commands
mkdir ./CIRCexplorer2/${5}@${s_name}
cd ./CIRCexplorer2/${5}@${s_name}
echo $s_name  CIRCexplorer2 parse start  `date +'%Y%m%d:%H:%M:%S'`  >> ${3}/time.txt
CIRCexplorer2 parse -t STAR ../../star-circ-bam/${5}@${s_name}/${s_name}Chimeric.out.junction > ${s_name}"_CIRCexplorer2_parse.log"
echo $s_name  CIRCexplorer2 parse end  `date +'%Y%m%d:%H:%M:%S'`  >> ${3}/time.txt

header=`sed -n -r '1s/(.).*/\1/p' back_spliced_junction.bed` ###v2.3.2输出方式
if [ ${header} != 'c' ]
then
awk '{print "chr"$1, $2, $3, $4, $5, $6}' back_spliced_junction.bed > ${s_name}"_back_spliced_junction.chr.bed"
else
mv back_spliced_junction.bed  ${s_name}"_back_spliced_junction.chr.bed"
fi

#等待一段时间并监控内存是否大于10GB
time_sleep3=$[$RANDOM%60+1]
sleep ${time_sleep3}
while true
do
free_7=` free -g | grep 'buffers/cache' | grep -v 'grep' | awk '{print $4;}' `
if [ ${free_7} -ge 10 ]
then
break
fi
sleep 10
echo waiting 10 s  start CIRCexplorer2 anno
done

#CIRCexplorer2 anno commands
echo $s_name  CIRCexplorer2 anno start  `date +'%Y%m%d:%H:%M:%S'`  >> ${3}/time.txt
CIRCexplorer2 annotate -r ${circ_ref_anno} -g ${genome_fa} -b ${s_name}"_back_spliced_junction.chr.bed" -o ${s_name}"_circularRNA_known.txt" > ${s_name}"_CIRCexplorer2_annotate.log"
echo $s_name  CIRCexplorer2 anno end  `date +'%Y%m%d:%H:%M:%S'`  >> ${3}/time.txt
cd ${3}

#任务第一次完成计数
echo $s_name >> ${3}/s_count.txt

time_sleep4=$[$RANDOM%10+1]
sleep ${time_sleep4}

fi
#测试用  CIRCexplorer2只提供human & mouse相关注释文件, 若非以上物种，直接使用下面代码-find circ.


if [ ${1} = cricetulus ];
then
## find circ 


#设置睡眠，睡眠时间随机确定，取值范围在1—101秒之内;此处的意义为避免多个STAR程序（单个内存>30G）同时运行造成服务器内存耗尽。
time_sleep=$[$RANDOM%100+1]
echo $s_name first_sleep ${time_sleep} >> ${3}/BOWTIE2_waiting.txt
echo waiting $s_name first_sleep ${time_sleep} s BOWTIE2
sleep ${time_sleep}

#剩余可用内存大于star_mem,并且已经运行的STAR程序数量小于star_have时，进入STAR运行程序；
while true
do
free_7=` free -g | grep 'buffers/cache' | grep -v 'grep' | awk '{print $4;}' `
if [ ${free_7} -ge ${bowtie2_mem} ]
then
bowtie2_count=`ps -ef | grep bowtie2 | grep -v "grep" | awk '{print $2;}' | wc -l`
if [ ${bowtie2_count} -lt ${bowtie2_have} ]
then
break 
fi
unset bowtie2_count
fi
sleep 10
done

#第二次睡眠，睡眠时间随机确定，取值范围在1—101秒之内;
time_sleep2=$[$RANDOM%200+50]
echo $s_name second_sleep ${time_sleep2} >>  ${3}/BOWTIE2_waiting.txt
echo waiting $s_name second_sleep  ${time_sleep2} s BOWTIE2
sleep ${time_sleep2}
while true
do
free_7=` free -g | grep 'buffers/cache' | grep -v 'grep' | awk '{print $4;}' `
if [ ${free_7} -ge ${bowtie2_mem} ]
then
bowtie2_count=`ps -ef | grep bowtie2 | grep -v "grep" | awk '{print $2;}' | wc -l`
if [ ${bowtie2_count} -lt ${bowtie2_have} ]
then
echo ${s_name}  bowtie2-start
echo ${s_name}  bowtie2-start  `date +'%Y%m%d:%H:%M:%S'`  >>  ${3}/time.txt 
break 
fi
unset bowtie2_count
fi
sleep ${time_sleep2}
done


#记录单个样本开始时间
echo ${s_name}-find_circ start `date +'%Y%m%d:%H:%M:%S'`  >> ${3}/time.txt

#Bowtie2运行,以及单双端判断：

if [ ${pair} -eq 2 ];then
mkdir ./bowtie2_circ_bam/${5}@${s_name} #第一次比对bam文件和未比对fa文件
#STAR --runThreadN 8 --chimSegmentMin 10 --genomeDir ${genome_star_dir} --readFilesIn ${s}/trim/*val_1.fq ${s}/trim/*val_2.fq --outFileNamePrefix  ./star_circ/${s_name}/${s_name}"_" > ./star_circ/${s_name}/${s_name}".star.log"
bowtie2 -p 8 --very-sensitive --mm -M20 --score-min=C,-15,0 -x ${genome_bowtie2_index} -q -1 ${s}/trim/*val_1.fq -2 ${s}/trim/*val_2.fq 2>./bowtie2_circ_bam/${5}@${s_name}/${s_name}.step1_bowtie2.log | samtools view -hbuS - | samtools sort -o ./bowtie2_circ_bam/${5}@${s_name}/${s_name}".bam"

samtools view -hf 4 ./bowtie2_circ_bam/${5}@${s_name}/${s_name}".bam" | samtools view -Sb - > ./bowtie2_circ_bam/${5}@${s_name}/${s_name}".unmapped.bam"

python2 /home/yaoyh/sofeware/find_circ/unmapped2anchors.py ./bowtie2_circ_bam/${5}@${s_name}/${s_name}".unmapped.bam" | gzip > ./bowtie2_circ_bam/${5}@${s_name}/${s_name}".anchors.qfa.gz"

echo $s_name  bowtie2-step1-end  `date +'%Y%m%d:%H:%M:%S'`  >> ${3}/time.txt

mkdir ./result/2.Alignment/Log/${s_name}
cp ./bowtie2_circ_bam/${5}@${s_name}/*step1_bowtie2.log ./result/2.Alignment/Log/${s_name}/ 
fi

if [ ${pair} -eq 1 ];then
mkdir ./bowtie2_circ_bam/${5}@${s_name}
#STAR --runThreadN 8 --chimSegmentMin 10 --genomeDir ${genome_star_dir} --readFilesIn ${s}/trim/*trimmed.fq  --outFileNamePrefix   ./star_circ/${s_name}/${s_name}"_"
bowtie2 -p 8 --very-sensitive --mm -M20 --score-min=C,-15,0 -x ${genome_bowtie2_index} -U ${s}/trim/*val.fq  | samtools view -hbuS - | samtools sort - output
samtools view -hf 4 output.bam | samtools view -Sb - > unmapped.bam
python2 /home/yaoyh/sofeware/find_circ/unmapped2anchors.py unmapped.bam | gzip > anchors.qfa.gz
echo $s_name  bowtie2-step1-end  `date +'%Y%m%d:%H:%M:%S'`  >> ${3}time.txt
fi

#cp ./bowtie2_circ_bam/${5}@${s_name}/${s_name}.step1_bowtie2.log ./result/2.Alignment/Log 

#剩余可用内存大于10G时通过，每隔10s执行一次
sleep 10
while true
do
free_7=` free -g | grep 'buffers/cache' | grep -v 'grep' | awk '{print $4;}' `
if [ ${free_7} -ge 10 ]
then
break
fi
sleep 10
done

###circular detection based on re-alignment
mkdir ./find_circ/${5}@${s_name}
#cd ./CIRCexplorer2/${s_name}
echo $s_name  re-alignment start  `date +'%Y%m%d:%H:%M:%S'`  >> ${3}/time.txt
#CIRCexplorer2 parse -t STAR ../../star_circ/${s_name}/${s_name}"_Chimeric.out.junction" > ${s_name}"_CIRCexplorer2_parse.log"
bowtie2 -p 8 --reorder --mm -M20 --score-min=C,-15,0 -q -x ${genome_bowtie2_index} -U ./bowtie2_circ_bam/${5}@${s_name}/${s_name}".anchors.qfa.gz" | python2 /home/yaoyh/sofeware/find_circ/find_circ.py -G ${genome_split_fa} -s ./find_circ/${5}@${s_name}/${s_name}".find_circ.sites.log" \
> ./find_circ/${5}@${s_name}/${s_name}".find_circ.sites.bed" 2> ./find_circ/${5}@${s_name}/${s_name}".find_circ.sites.reads"

echo $s_name  re-alignment end  `date +'%Y%m%d:%H:%M:%S'`  >> ${3}/time.txt

###Filter the outputs
echo $s_name  Final Filter start  `date +'%Y%m%d:%H:%M:%S'`  >>${3}/time.txt
#CIRCexplorer2 annotate -r ${circ_ref_anno} -g ${genome_fa} -b ${s_name}"_back_spliced_junction.chr.bed" -o ${s_name}"_circularRNA_known.txt" > ${s_name}"_CIRCexplorer2_annotate.log"
grep circ ./find_circ/${5}@${s_name}/${s_name}".find_circ.sites.bed" | grep -v chrM | python2 /home/yaoyh/sofeware/find_circ/sum.py -2,3 | python2 /home/yaoyh/sofeware/find_circ/scorethresh.py -16 1 | \
python2 /home/yaoyh/sofeware/find_circ/scorethresh.py -15 2 | python2 /home/yaoyh/sofeware/find_circ/scorethresh.py -14 2 | python2 /home/yaoyh/sofeware/find_circ/scorethresh.py 7 2 | \
python2 /home/yaoyh/sofeware/find_circ/scorethresh.py 8 40 | python2 /home/yaoyh/sofeware/find_circ/scorethresh.py 9 40 | python2 /home/yaoyh/sofeware/find_circ/scorethresh.py -17 100000 > ./find_circ/${5}@${s_name}/${s_name}".findcirc_40x40.bed"

grep circ ./find_circ/${5}@${s_name}/${s_name}".find_circ.sites.bed" | grep -v chrM | python2 /home/yaoyh/sofeware/find_circ/sum.py -2,3 | python2 /home/yaoyh/sofeware/find_circ/scorethresh.py -16 1 | \
python2 /home/yaoyh/sofeware/find_circ/scorethresh.py -15 2 | python2 /home/yaoyh/sofeware/find_circ/scorethresh.py -14 2 | python2 /home/yaoyh/sofeware/find_circ/scorethresh.py 7 2 | \
python2 /home/yaoyh/sofeware/find_circ/scorethresh.py 8,9 35 | python2 /home/yaoyh/sofeware/find_circ/scorethresh.py -17 100000 > ./find_circ/${5}@${s_name}/${s_name}".findcirc_filter.bed"
echo $s_name  Final Filter end  `date +'%Y%m%d:%H:%M:%S'`  >>${3}/time.txt


#任务第一次完成计数
echo $s_name >> ${3}/s_count.txt

fi


#-------
#--------------------
#----------------------------------      五
#----------------------------------------------
#----------------------------------
#--------------------
#-------

#删除占用空间大的文件
#释放变量
unset a
unset s
unset pair
unset species
unset genome_dir

