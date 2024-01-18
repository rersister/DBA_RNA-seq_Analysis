
#!/bin/sh 
#PBS -q core24
#PBS -l nodes=1:ppn=5,mem=50gb,walltime=40:00:00 
#HSCHED -s hschedd 

out_path='~/DBA_BulkRNA_10/GATK-RNA-seq2call_variants/5.findRP_variants'
file_path='~/DBA_BulkRNA_10/GATK-RNA-seq2call_variants/4.variants-calling' #输出文件路径outdir=$5 #输出文件路径
wk_path='~/DBA_BulkRNA_10/GATK-RNA-seq2call_variants'
cd ${wk_path}

cat sample.txt |while read id
do
	sample=${id}
	grep '\<RPS' ${file_path}/${sample}/variants/${sample}.geneanno.variant_function > ${out_path}/${id}.geneanno.variant_function_findRPS
	grep '\<RPL' ${file_path}/${sample}/variants/${sample}.geneanno.variant_function > ${out_path}/${id}.geneanno.variant_function_findRPL
	grep '\<RPS' ${file_path}/${sample}/variants/${sample}.geneanno.exonic_variant_function > ${out_path}/${id}.geneanno.exonic_variant_function_findRPS
	grep '\<RPL' ${file_path}/${sample}/variants/${sample}.geneanno.exonic_variant_function > ${out_path}/${id}.geneanno.exonic_variant_function_findRPL
done
