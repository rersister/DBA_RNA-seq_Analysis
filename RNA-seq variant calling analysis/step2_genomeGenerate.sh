#!/bin/sh 
#PBS -q q512G
#PBS -l nodes=1:ppn=5,mem=100gb,walltime=40:00:00 
#HSCHED -s hschedd 

STAR_path="~/miniconda/bin/STAR"
STAR_index="~/red_cell_pub_up_Anal/STAR/GRCh37/star-index"
refer_fast='~/red_cell_pub_up_Anal/STAR/GRCh37'
wk_path='~/DBA_BulkRNA_10/GATK-RNA-seq2call_variants'
dataDir='~/DBA_BulkRNA_10'
STAR2genome='~/DBA_BulkRNA_10/GATK-RNA-seq2call_variants/2.star2genome'
cd ${wk_path}



###star 进行map
cat sample.txt |while read id
do
	dir=${STAR2genome}/${id}_index
	#output diretory
	if [ ! -d ${dir} ]
		then mkdir -p ${dir}
	fi
	$STAR_path --runThreadN 40 \
	--runMode genomeGenerate \
	--genomeDir ${dir} \
	--genomeFastaFiles ${refer_fast}/GRCh37.primary_assembly.genome.fa \
	--sjdbGTFfile ${refer_fast}/Homo_sapiens.GRCh37.gtf \
	--sjdbFileChrStartEnd ${wk_path}/1.mapping/${id}SJ.out.tab \
	--sjdbOverhang 140 \
	--limitSjdbInsertNsj 1100000
done
