
#!/bin/sh 
#PBS -q core24
#PBS -l nodes=1:ppn=5,mem=50gb,walltime=40:00:00 
#HSCHED -s hschedd 

wk_path='~/DBA_BulkRNA_10/GATK-RNA-seq2call_variants/Anovar/annovar'
outdir='~/DBA_BulkRNA_10/GATK-RNA-seq2call_variants/4.variants-calling' #输出文件路径outdir=$5 #输出文件路径


cat sample.txt |while read id
do
	cd ${wk_path}
	sample=${id}
	vcf_file=${outdir}/${sample}/variants/${sample}.snp.filter.vcf.gz
	perl convert2annovar.pl -format vcf4 ${vcf_file} -out ${outdir}/${sample}/variants/${sample}.filter.vcf.avinput
	avinput=${outdir}/${sample}/variants/${sample}.filter.vcf.avinput
	perl annotate_variation.pl -geneanno -dbtype refGene \
		-out ${outdir}/${sample}/variants/${sample}.geneanno \
		--buildver hg19 ${avinput} humandb/
done
