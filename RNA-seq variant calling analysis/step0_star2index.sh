#!/bin/sh 
#PBS -q q512G
#PBS -l nodes=1:ppn=5,mem=100gb,walltime=40:00:00 
#HSCHED -s hschedd 


##利用star 构建参考基因组 index 脚本
STAR_path="/xtdisk/fangxd_group/liangyue/miniconda/bin/STAR"
file_path="/xtdisk/fangxd_group/zhenggm/red_cell_pub_up_Anal/STAR/GRCh37"
index_path='/xtdisk/fangxd_group/zhenggm/red_cell_pub_up_Anal/STAR/GRCh37'
## gtf 注释文件下载  https://www.gencodegenes.org/human/release_43lift37.html
## fasta 文件下载  nohup wget  --no-check-certificate https://193.62.193.139/pub/grch37/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz &
cd ${index_path}
#gzip 解压
#gzip -c -d Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz > Homo_sapiens.GRCh38.fa
#gzip -c -d Homo_sapiens.GRCh38.99.gtf.gz > Homo_sapiens.GRCh38.99.gtf 

$STAR_path --runMode genomeGenerate \
	--runThreadN 10 \
	--genomeDir ./star-index \
	--genomeFastaFiles ${index_path}/GRCh37.primary_assembly.genome.fa \
	--sjdbGTFfile ${file_path}/Homo_sapiens.GRCh37.gtf \
	--sjdbOverhang 100