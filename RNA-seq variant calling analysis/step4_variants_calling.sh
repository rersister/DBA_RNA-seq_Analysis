#!/bin/sh 
#PBS -q q512G
#PBS -l nodes=1:ppn=5,mem=100gb,walltime=40:00:00 
#HSCHED -s hschedd 

gatk='/pnas/fangxd_group/zhenggm/software/yes/bin/gatk'
samtools='/pnas/fangxd_group/zhenggm/software/yes/envs/samtools-1.16/bin/samtools'
ref_genome='/xtdisk/fangxd_group/zhenggm/red_cell_pub_up_Anal/STAR/GRCh37/GRCh37.primary_assembly.genome.fa'
dbsnp='/xtdisk/fangxd_group/zhenggm/DBA_BulkRNA_10/GATK-RNA-seq2call_variants/dbsnp/all_coChr.vcf'
pre_BAM='/xtdisk/fangxd_group/zhenggm/DBA_BulkRNA_10/GATK-RNA-seq2call_variants/3.mapping-2'
outdir_pre='/xtdisk/fangxd_group/zhenggm/DBA_BulkRNA_10/GATK-RNA-seq2call_variants/4.variants-calling' #输出文件路径outdir=$5 #输出文件路径

cat sample.txt |while read id
do
	BAM=${pre_BAM}/${id}Aligned.sortedByCoord.out.bam
	RGID=${id}
	library=${id}
	sample=${id}
	outdir=${outdir_pre}/${id}
	#output diretory
	if [ ! -d ${outdir}/variants ]
		then mkdir -p ${outdir}/variants
	fi
	if [ ! -d ${outdir}/processBAM ]
		then mkdir -p ${outdir}/processBAM
	fi

	#add RG & sort
	time $gatk AddOrReplaceReadGroups \
		--INPUT $BAM \
		--OUTPUT $outdir/processBAM/${sample}.RG.bam \
		--RGID $RGID \
		--RGLB $library \
		--RGPL illumina \
		--RGPU snpcall \
		--RGSM ${sample} && echo "** ADD RG done **"

	#mark duplicates & index
	time $gatk MarkDuplicates \
		-I $outdir/processBAM/${sample}.RG.bam \
		-M $outdir/processBAM/${sample}.markdup_metrics.txt \
		-O $outdir/processBAM/${sample}.RG.markdup.bam &&\
	time $samtools index $outdir/processBAM/${sample}.RG.markdup.bam && echo "** ${sample}.RG.markdup.bam index done **"

	#split N Trim
	time $gatk SplitNCigarReads \
		--R $ref_genome \
		--I $outdir/processBAM/${sample}.RG.markdup.bam \
		--O $outdir/processBAM/${sample}.RG.markdup.split.bam && echo "** ${sample}.RG.markdup.bam split N done **"

	#HaplotypeCaller
	time $gatk HaplotypeCaller \
		-R $ref_genome \
		-I $outdir/processBAM/${sample}.RG.markdup.split.bam \
		-O $outdir/variants/${sample}.vcf.gz \
		-dont-use-soft-clipped-bases \
		--standard-min-confidence-threshold-for-calling 20 \
		--dbsnp $dbsnp && echo "** ${sample}.vcf done **"


	#选择SNP，并filter  #过滤的阈值可以根据实验进行设置
	time $gatk SelectVariants \
		-select-type SNP \
		-V $outdir/variants/${sample}.vcf.gz \
		-O $outdir/variants/${sample}.snp.vcf.gz && \

	time $gatk VariantFiltration \
		-V $outdir/variants/${sample}.vcf.gz \
		-O $outdir/variants/${sample}.snp.filter.vcf.gz \
		-R $ref_genome \
		--filter-expression 'QD < 2.0 || FS > 60.0 || MQ <25.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0' \
		--filter-name "my_filter" && echo "** ${sample}.snp filter done **"
done