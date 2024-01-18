#!/bin/sh 
#PBS -q core40
#PBS -l nodes=1:ppn=5,mem=80gb,walltime=400:00:00 
#HSCHED -s hschedd 


STAR_path="/xtdisk/fangxd_group/liangyue/miniconda/bin/STAR"
STAR_index="/xtdisk//fangxd_group/zhenggm/red_cell_pub_up_Anal/STAR/GRCh37/star-index"
wk_path='/xtdisk/fangxd_group/zhenggm/DBA_BulkRNA_10/GATK-RNA-seq2call_variants'
dataDir='/xtdisk/fangxd_group/zhenggm/DBA_BulkRNA_10'
cd ${wk_path}

###star 进行map
cat sample.txt |while read id
do
	arr=${id}
	fq1=${dataDir}/step1_trim/${id}_1_val_1.fq.gz
	fq2=${dataDir}/step1_trim/${id}_2_val_2.fq.gz
	$STAR_path --runMode alignReads \
		--runThreadN 20 \
		--genomeDir ${STAR_index} \
		--outFileNamePrefix ${wk_path}/1.mapping/${id} \
		--outSAMmapqUnique 255 \
		--outSAMtype BAM SortedByCoordinate \
		--outBAMsortingThreadN 10 \
		--readFilesCommand zcat \
		--readFilesIn $fq1 $fq2
done

