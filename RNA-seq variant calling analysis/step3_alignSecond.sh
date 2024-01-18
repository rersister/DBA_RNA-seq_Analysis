#!/bin/sh 
#PBS -q q512G
#PBS -l nodes=1:ppn=5,mem=100gb,walltime=40:00:00 
#HSCHED -s hschedd 


##利用star 构建参考基因组 index 脚本
STAR_path="~/miniconda/bin/STAR"
STAR2genome='~/DBA_BulkRNA_10/GATK-RNA-seq2call_variants/2.star2genome'
wk_path='~/DBA_BulkRNA_10/GATK-RNA-seq2call_variants'
dataDir='~/DBA_BulkRNA_10'
cd ${wk_path}

#gzip 解压
#gzip -c -d Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz > Homo_sapiens.GRCh38.fa
#gzip -c -d Homo_sapiens.GRCh38.99.gtf.gz > Homo_sapiens.GRCh38.99.gtf 

cat sample.txt |while read id
do
	dir=${STAR2genome}/${id}_index
	#output diretory
	if [ ! -d ${dir} ]
		then mkdir -p ${dir}
	fi
	# #上一步生成的index路径  #mapping结果的路径及前缀 #唯一比对的质量值设为255 #直接输出排序好的BAM文件
	arr=${id}
	fq1=${dataDir}/step1_trim/${id}_1_val_1.fq.gz
	fq2=${dataDir}/step1_trim/${id}_2_val_2.fq.gz
	$STAR_path --runMode alignReads \
		--runThreadN 20 \
		--genomeDir ${dir} \
		--outFileNamePrefix ${wk_path}/3.mapping-2/${arr} \
		--outSAMmapqUnique 255 \
		--outSAMtype BAM SortedByCoordinate \
		--outBAMsortingThreadN 10 \
		--readFilesCommand zcat \
		--readFilesIn $fq1 $fq2
done
