#!/bin/bash
echo ==========start at : `date` ==========
# HISAT alignment
export LD_LIBRARY_PATH=/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_RNAseq/RNA_RNAseq_2017a/Alignment/../software/RNA_lib:$LD_LIBRARY_PATH && \
cd /hwfssz5/ST_MCHRI/COHORT/USER/zhongjixing/DPW/RNA-seq/result_2/process/GenomeMapping_HISAT/Siha1 && \
/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_RNAseq/RNA_RNAseq_2017a/Alignment/../software/hisat2-2.0.4/hisat2 --phred64 --sensitive --no-discordant --no-mixed -I 1 -X 1000 -x /ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_RNAref/RNA_RNAref_version5.0_beta/Database/hg19/GenomeHisat2Index/chrALL -1 /hwfssz5/ST_MCHRI/COHORT/USER/zhongjixing/DPW/RNA-seq/result_2/upload/CleanData/Siha1_1.fq.gz -2 /hwfssz5/ST_MCHRI/COHORT/USER/zhongjixing/DPW/RNA-seq/result_2/upload/CleanData/Siha1_2.fq.gz 2>Siha1.Map2GenomeStat.xls | /ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_RNAseq/RNA_RNAseq_2017a/Alignment/../software/samtools view -b -S -o Siha1.bam - && \

# AddRG & Reorder & Sort BAM for downstream analysis
if [ ! -d java_tmp ];then mkdir -p java_tmp;fi && \
/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_RNAseq/RNA_RNAseq_2017a/Alignment/../software/java -Xmx4G -Djava.io.tmpdir=java_tmp -jar /hwfssz1/ST_MCHRI/STEMCELL/USER/wuliang2/biosoftware/picard.jar AddOrReplaceReadGroups I=Siha1.bam O=Siha1.AddRG.bam RGID=Siha1 RGLB=Siha1_library RGPL=illumina RGPU=machine RGSM=Siha1 VALIDATION_STRINGENCY=SILENT && \
/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_RNAseq/RNA_RNAseq_2017a/Alignment/../software/java -Xmx4G -Djava.io.tmpdir=java_tmp -jar /hwfssz1/ST_MCHRI/STEMCELL/USER/wuliang2/biosoftware/picard.jar ReorderSam I=Siha1.AddRG.bam O=Siha1.AddRG.Reorder.bam R=/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_RNAref/RNA_RNAref_version5.0_beta/Database/hg19/GenomeGatkIndex/chrALL.sort.fa VALIDATION_STRINGENCY=SILENT && \
/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_RNAseq/RNA_RNAseq_2017a/Alignment/../software/java -Xmx4G -Djava.io.tmpdir=java_tmp -jar /hwfssz1/ST_MCHRI/STEMCELL/USER/wuliang2/biosoftware/picard.jar SortSam I=Siha1.AddRG.Reorder.bam O=Siha1.AddRG.Reorder.Sort.bam SO=coordinate VALIDATION_STRINGENCY=SILENT && \
/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_RNAseq/RNA_RNAseq_2017a/Alignment/../software/samtools index Siha1.AddRG.Reorder.Sort.bam && \
mv Siha1.AddRG.Reorder.Sort.bam Siha1.AddRG.Reorder.Sort.bam.bai /hwfssz5/ST_MCHRI/COHORT/USER/zhongjixing/DPW/RNA-seq/result_2/upload/IGV/bam && \
rm -rf java_tmp && \
rm Siha1.bam Siha1.AddRG.bam Siha1.AddRG.Reorder.bam && \
cp Siha1.Map2GenomeStat.xls /hwfssz5/ST_MCHRI/COHORT/USER/zhongjixing/DPW/RNA-seq/result_2/Analysis_Report/BGI_result/MapStat/GenomeMapping  && \
echo ==========end at : `date` ========== && \
echo Still_waters_run_deep 1>&2 && \
echo Still_waters_run_deep > /hwfssz5/ST_MCHRI/COHORT/USER/zhongjixing/DPW/RNA-seq/result_2/shell/GenomeMapping_HISAT/Siha1/Alignment_Siha1.sh.sign
