#!/bin/bash
set -o pipefail
## step1:filter_all
SOAPnuke1.5.6 filter -n 0.1 -q 0.5 -i -l 12 -Q 2 -5 1 -E 50 -G -A 0.3 -1 Cleandata/SAMPLE_R1.fq.gz -2 Cleandata/SAMPLE_R2.fq.gz -f AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -r AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA -M 2 -o process/SiHa -C SAMPLE_1.fq -D SAMPLE_2.fq && \
fqcheck -r process/SiHa/SAMPLE_1.fq.gz -c process/SiHa/SAMPLE_1.fqcheck && \
fqcheck -r process/SiHa/SAMPLE_2.fq.gz -c process/SiHa/SAMPLE_2.fqcheck && \
perl fqcheck_distribute.pl process/SiHa/SAMPLE_1.fqcheck process/SiHa/SAMPLE_2.fqcheck -o process/SiHa/SAMPLE. && \
mv process/SiHa/SAMPLE_1.fq.gz process/SiHa/SAMPLE_2.fq.gz process/SiHa/*.png process/clean_data/ && \
if [-e process/SiHa/Raw_SAMPLE_1.fq.gz ];then rm -rf process/SiHa/Raw_SAMPLE_1.fq.gz;fi && \
if [-e process/SiHa/Raw_SAMPLE_2.fq.gz ];then rm -rf process/SiHa/Raw_SAMPLE_2.fq.gz;fi && \

## Check Filter Stat
perl soapnuke_stat.pl process/SiHa/Basic_Statistics_of_Sequencing_Quality.txt process/SiHa/Statistics_of_Filtered_Reads.txt > process/SiHa/Siha.stat && \
perl catRS.pl process/SiHa_DNA.stat_all.list process/clean_data/SiHa_DNA.xls Siha_DNA && \

## Mapping using bwa
bwa mem -t 30 -M -R "@RG\tID:SiHa\tPL:ILLUMINA\tPU:SiHa_DNA\tLB:SiHa_DNA\tSM:SiHa_DNA\tCN:BGI" hg19.fasta process/clean_data/SAMPLE_1.fq.gz process/clean_data/SAMPLE_2.fq.gz | samtools view -Sb -o process/SiHa/SAMPLE.bam - && \
java -Xmx3g -Djava.io.tmpdir=process/java_tmp -XX:MaxPermSize=512m -XX:-UseGCOverheadLimit -jar picard.jar SortSam I=process/SiHa/SAMPLE.bam O=process/SiHa/SAMPLE.sort.bam SO=coordinate MAX_RECORDS_IN_RAM=1000000 && \
rm -rf process/SiHa/SAMPLE.bam && \

## rmdup_SiHa_DNA
java -Xmx2G -Djava.io.tmpdir=process/java_tmp -XX:MaxPermSize=512m -XX:-UseGCOverheadLimit -jar picard.jar MarkDuplicates REMOVE_DUPLICATES=true I=process/SiHa/SAMPLE.sort.bam O=result/result_alignment/SAMPLE.bam METRICS_FILE=process/SAMPLE.bam.mat TMP_DIR=process/tmp && \
java -Xmx2G -jar picard.jar BuildBamIndex I=result/result_alignment/SAMPLE.bam O=result/result_alignment/SAMPLE.bam.bai && \

## clean_laneBam
rm -rf process/SiHa/SAMPLE.sort.bam && \

##step6_bam_Stat
perl bamstats.pl -i /WES/outdir_hg19/result/result_alignment/SiHa_DNA.bam -r shell/bed_ex_region.sort.bed -o process/statistics -plot && \
cp process/statistics/cumuPlot.png result/result_alignment/SiHa_DNA.Cumulative.png && \
cp process/statistics/histPlot.png result/result_alignment/SiHa_DNA.Depth.png && \
cp process/statistics/information.xls result/result_alignment/SiHa_DNA.SummaryTable.xls && \
cp process/statistics/insert.png result/result_alignment/SiHa_DNA.Insert.png && \

##step7_GATKreallnRecal
java -Xmx4g -Djava.io.tmpdir=process/java_tmp -jar GenomeAnalysisTK.jar -T RealignerTargetCreator -R hg19.fasta -I result/result_alignment/Siha_DNA.bam -L shell/Siha_DNA/cvr_ex_region.sort.bed -known Database/hg19/gatk/Mills_and_1000G_gold_standard.indels.hg19.vcf -known Database/hg19/gatk/1000G_phase1.indels.hg19.vcf -o process/Siha_DNA/Siha_DNA.intervals && \
java -Xmx4g -Djava.io.tmpdir=process/java_tmp -jar GenomeAnalysisTK.jar -T IndelRealigner -R hg19.fasta -I result/result_alignment/Siha_DNA.bam -known Database/hg19/gatk/Mills_and_1000G_gold_standard.indels.hg19.vcf -known Database/hg19/gatk/1000G_phase1.indels.hg19.vcf -targetIntervals process/Siha_DNA/Siha_DNA.intervals -o process/Siha_DNA/Siha_DNA.realign.bam && \
java -Xmx4g -Djava.io.tmpdir=process/java_tmp -jar GenomeAnalysisTK.jar -nct 5 -T BaseRecalibrator -R hg19.fasta -I process/Siha_DNA.realign.bam -knownSites Database/hg19/gatk/dbsnp_138.hg19.vcf -knownSites Database/hg19/gatk/Mills_and_1000G_gold_standard.indels.hg19.vcf -knownSites Database/hg19/gatk/1000G_phase1.indels.hg19.vcf -o process/Siha_DNA/Siha_DNA.realign.recal.table && \
java -Xmx4g -Djava.io.tmpdir=process/java_tmp -jar GenomeAnalysisTK.jar -nct 5 -T PrintReads -R hg19.fasta -I process/Siha_DNA.realign.bam -BQSR process/Siha_DNA/Siha_DNA.realign.recal.table -o process/Siha_DNA/Siha_DNA.realign.recal.bam && \

##step8_clean_Siha_DNA
rm -f process/Siha_DNA/Siha_DNA.intervals process/Siha_DNA/Siha_DNA.realign.bam process/Siha_DNA/Siha_DNA.realign.bai && \

##step9_callGVCF_GATK
java -Xmx5G -Djava.io.tmpdir=process/java_tmp -jar GenomeAnalysisTK.jar -T HaplotypeCaller -R hg19.fasta -I process/Siha_DNA/Siha_DNA.realign.recal.bam -L shell/Siha_DNA/cvr_ex_region.sort.bed --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 -o process/Siha_DNA/callGVCF_GATK/Siha_DNA.g.vcf.gz && \
java -Xmx2G -Djava.io.tmpdir=process/java_tmp -jar GenomeAnalysisTK.jar -T GenotypeGVCFs -R hg19.fasta --variant process/Siha_DNA/callGVCF_GATK/Siha_DNA.g.vcf.gz -o process/Siha_DNA/callGVCF_GATK/Siha_DNA.vcf.gz -stand_call_conf 10 -allSites && \
rm -rf process/Siha_DNA/Siha_DNA.realign.recal.bam && \

##step10_GATK_snp
export PERL5LIB="/ifs4/BC_PUB/biosoft/pipeline/DNA/DNA_Human_WES/DNA_Human_WES_2016b/lib:$PERL5LIB"
java -Xmx3G -Djava.io.tmpdir=process/java_tmp -jar GenomeAnalysisTK.jar  -T SelectVariants -R hg19.fasta -V /WES/outdir_hg19/process/Siha_DNA/callGVCF_GATK/Siha_DNA.vcf.gz -selectType SNP --excludeNonVariants -o /WES/outdir_hg19/process/Siha_DNA/snp_GATK/Siha_DNA.raw.snp.vcf.gz && \
java -Xmx3G -Djava.io.tmpdir=process/java_tmp -jar GenomeAnalysisTK.jar -T VariantFiltration -R hg19.fasta -V /WES/outdir_hg19/process/Siha_DNA/snp_GATK/Siha_DNA.raw.snp.vcf.gz \
--filterExpression "QD < 2.0 || FS > 60.0 || MQ <40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filterName "filter"  -o /WES/outdir_hg19/process/Siha_DNA/snp_GATK/Siha_DNA.filtered_snp.vcf.gz && \
/software/bcftools-1.2/bcftools view -e 'ALT=="*"' -f "PASS" -o /WES/outdir_hg19/process/Siha_DNA/snp_GATK/anno/Siha_DNA.filtered_snp.vcf.gz -O z /WES/outdir_hg19/process/Siha_DNA/snp_GATK/Siha_DNA.filtered_snp.vcf.gz && \
perl /software/annodb-3.3.0/annodb.pl --mode snp --if vcf --optdb 'Conservative,Cancer,Disease,Functional,ENCODE,Population' --genedb 'proteinatlas,omim,pharmgkb,dbnsfp,cgc' --remove --verdbsnp v149_hg19 --buildver hg19 --outfile /WES/outdir_hg19/process/Siha_DNA/snp_GATK/anno/Siha_DNA.filtered_snp.vcf /WES/outdir_hg19/process/Siha_DNA/snp_GATK/anno/Siha_DNA.filtered_snp.vcf.gz && \
perl /software/annodb-3.3.0/Annodb_stat_for_all.pl -i /WES/outdir_hg19/process/Siha_DNA/snp_GATK/anno/Siha_DNA.filtered_snp.vcf.genome_summary.csv -o /WES/outdir_hg19/process/Siha_DNA/snp_GATK/anno/Siha_DNA.snp.stat -s Siha_DNA -v snp && \
mkdir -p /WES/outdir_hg19/result/Siha_DNA/result_variation/snp && \
cp /WES/outdir_hg19/process/Siha_DNA/snp_GATK/anno/Siha_DNA.filtered_snp.vcf.gz /WES/outdir_hg19/result/Siha_DNA/result_variation/snp/Siha_DNA.snp.vcf.gz && \
cp /WES/outdir_hg19/process/Siha_DNA/snp_GATK/anno/Siha_DNA.filtered_snp.vcf.genome_summary.csv /WES/outdir_hg19/result/Siha_DNA/result_variation/snp/Siha_DNA.snp.annot.csv && \
cp /WES/outdir_hg19/process/Siha_DNA/snp_GATK/anno/Siha_DNA.snp.stat.all_stat /WES//outdir_hg19/result/Siha_DNA/result_variation/snp/Siha_DNA.snp.AnnotationTable.xls && \
cp /WES/outdir_hg19/process/Siha_DNA/snp_GATK/anno/Siha_DNA.snp.stat.coding_stat /WES/outdir_hg19/result/Siha_DNA/result_variation/snp/Siha_DNA.coding.snp.AnnotationTable.xls && \
cp /WES/outdir_hg19/process/Siha_DNA/snp_GATK/anno/Siha_DNA.filtered_snp.vcf.exome_summary.csv /WES/outdir_hg19/result/Siha_DNA/result_variation/snp/Siha_DNA.snp.cds_annot.csv && \
echo ==========end at : `date` ========== && \
echo -|awk -v S=$SECONDS '{printf "task run time\t%02d:%02d:%02d\n",S/(60*60),S%(60*60)/60,S%60}' && \
echo Still_waters_run_deep 1>&2 && \
echo Still_waters_run_deep > /WES//outdir_hg19/shell/Siha_DNA/snp_GATK/GATK_snp.sh.sign

echo ==========start at : `date` ==========
rm -f /WES/outdir_hg19/process/Siha_DNA/snp_GATK/Siha_DNA.raw.snp.vcf.gz && \
rm -f /WES/outdir_hg19/process/Siha_DNA/snp_GATK/Siha_DNA.filtered_snp.vcf.gz && \
echo ==========end at : `date` ========== && \
echo -|awk -v S=$SECONDS '{printf "task run time\t%02d:%02d:%02d\n",S/(60*60),S%(60*60)/60,S%60}' && \
echo Still_waters_run_deep 1>&2 && \
echo Still_waters_run_deep > /WES/outdir_hg19/shell/Siha_DNA/snp_GATK/rm.sh.sign

echo ==========start at : `date` ==========
export PERL5LIB=/ifs4/BC_PUB/biosoft/pipeline/DNA/DNA_Human_WES/DNA_Human_WES_2016b/lib\
:/ifs4/BC_PUB/biosoft/pipeline/DNA/DNA_HiSeqExome/DNA_HiSeqExome_2015a/bin/annodb/lib/lib/perl5/x86_64-linux-thread-multi/auto/Bio/DB\
:/hwfssz1/ST_MCHRI/CLINIC/SOFTWARES/lib/perl5/5.24.1/\
:/hwfssz1/ST_MCHRI/CLINIC/SOFTWARES/lib/perl5/5.24.1/x86_64-linux/\
:$PERL5LIB
export PATH="/ifs4/BC_PUB/biosoft/pipeline/DNA/DNA_Human_WES/DNA_Human_WES_2016b/bin:$PATH"
/software/java/jre1.8.0_101/bin/java -Xmx3G -Djava.io.tmpdir=/WES/outdir_hg19/java_tmp -jar /software/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar -T VariantFiltration -R hg19.fasta -V /WES/outdir_hg19/process/Siha_DNA/indel_GATK/Siha_DNA.raw.indel.vcf.gz \
--filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" --filterName "filter"  -o /WES/outdir_hg19/process/Siha_DNA/indel_GATK/Siha_DNA.filtered_indel.vcf.gz && \
/software/bcftools-1.2/bcftools view -f "PASS" -o /WES/outdir_hg19/process/Siha_DNA/indel_GATK/anno/Siha_DNA.filtered_indel.vcf.gz -O z /WES/outdir_hg19/process/Siha_DNA/indel_GATK/Siha_DNA.filtered_indel.vcf.gz && \
perl /software/indel_stat.pl /WES/outdir_hg19/process/Siha_DNA/indel_GATK/anno/Siha_DNA.filtered_indel.vcf.genome_summary.csv /WES/outdir_hg19/process/Siha_DNA/indel_GATK/anno/Siha_DNA.indel_len && \
perl /software/indel_stat.pl /WES/outdir_hg19/process/Siha_DNA/indel_GATK/anno/Siha_DNA.filtered_indel.vcf.exome_summary.csv /WES/outdir_hg19/process/Siha_DNA/indel_GATK/anno/Siha_DNA.indel_cds_len && \
perl /software/indel_lenght_R.pl Siha_DNA /WES/outdir_hg19/process/Siha_DNA/indel_GATK/anno/Siha_DNA.indel_len.xls /WES/outdir_hg19/process/Siha_DNA/indel_GATK/anno/Siha_DNA.indel_cds_len.xls /WES/outdir_hg19/process/Siha_DNA/indel_GATK/anno && \
mkdir -p /WES/outdir_hg19/result/Siha_DNA/result_variation/indel && \
cp /WES/outdir_hg19/process/Siha_DNA/indel_GATK/anno/Siha_DNA.filtered_indel.vcf.gz /WES/outdir_hg19/result/Siha_DNA/result_variation/indel/Siha_DNA.indel.vcf.gz && \
cp /WES/outdir_hg19/process/Siha_DNA/indel_GATK/anno/Siha_DNA.filtered_indel.vcf.genome_summary.csv /WES/outdir_hg19/result/Siha_DNA/result_variation/indel/Siha_DNA.indel.annot.csv && \
cp /WES/outdir_hg19/process/Siha_DNA/indel_GATK/anno/Siha_DNA.filtered_indel.vcf.exome_summary.csv /WES/outdir_hg19/result/Siha_DNA/result_variation/indel/Siha_DNA.indel.cds_annot.csv && \
cp /WES/outdir_hg19/process/Siha_DNA/indel_GATK/anno/Siha_DNA.indel.stat.all_stat /WES/outdir_hg19/result/Siha_DNA/result_variation/indel/Siha_DNA.indel.AnnotationTable.xls && \
cp /WES/outdir_hg19/process/Siha_DNA/indel_GATK/anno/Siha_DNA.indel.stat.coding_stat /WES/outdir_hg19/result/Siha_DNA/result_variation/indel/Siha_DNA.coding.indel.AnnotationTable.xls && \
cp /WES/outdir_hg19/process/Siha_DNA/indel_GATK/anno/Siha_DNA.indel_cds_len.png /WES/outdir_hg19/result/Siha_DNA/result_variation/indel/Siha_DNA.indel.cds.len.png && \
echo ==========end at : `date` ========== && \
echo -|awk -v S=$SECONDS '{printf "task run time\t%02d:%02d:%02d\n",S/(60*60),S%(60*60)/60,S%60}' && \
echo Still_waters_run_deep 1>&2 && \
echo Still_waters_run_deep > /WES/outdir_hg19/shell/Siha_DNA/indel_GATK/GATK_indel.sh.sign

echo ==========start at : `date` ==========
rm -f /WES/outdir_hg19/process/Siha_DNA/indel_GATK/Siha_DNA.raw.indel.vcf.gz && \
rm -f /WES/outdir_hg19/process/Siha_DNA/indel_GATK/Siha_DNA.filtered_indel.vcf.gz /WES/outdir_hg19/process/Siha_DNA/indel_GATK/anno/Siha_DNA.indel*.xls && \
echo ==========end at : `date` ========== && \
echo -|awk -v S=$SECONDS '{printf "task run time\t%02d:%02d:%02d\n",S/(60*60),S%(60*60)/60,S%60}' && \
echo Still_waters_run_deep 1>&2 && \
echo Still_waters_run_deep > /WES/outdir_hg19/shell/Siha_DNA/indel_GATK/rm.sh.sign
