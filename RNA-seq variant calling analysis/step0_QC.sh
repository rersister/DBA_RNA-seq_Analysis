#!/bin/sh 
#PBS -q q512G
#PBS -l nodes=1:ppn=20,mem=80gb,walltime=40:00:00 
#HSCHED -s hschedd 

# author zhenggm  QC为 fastQC  的安装地址
QC="/xtdisk/fangxd_group/liangyue/miniconda/bin/fastqc"
wkdir="/xtdisk/fangxd_group/zhenggm/DBA_BulkRNA_10"
dataDir="/xtdisk/fangxd_group/zhenggm/DBA_BulkRNA_10/rawdata"
cd ${wkdir}


# -o 表示输出路径   后面是输入fq的文件路径
cat sample.txt |while read id 
do 
	$QC -o ${wkdir}/step_0_QC ${dataDir}/${id}_1.fq.gz;
	$QC -o ${wkdir}/step_0_QC ${dataDir}/${id}_2.fq.gz
	echo "**** ${id} has been done ****"
done