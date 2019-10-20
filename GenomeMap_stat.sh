#!/bin/bash
echo ==========start at : `date` ==========
perl /ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_RNAseq/RNA_RNAseq_2017a/Alignment/MapStat.pl -indir /hwfssz5/ST_MCHRI/COHORT/USER/zhongjixing/DPW/RNA-seq/result_2/process/GenomeMapping_HISAT -output /hwfssz5/ST_MCHRI/COHORT/USER/zhongjixing/DPW/RNA-seq/result_2/process/GenomeMapping_HISAT/GenomeMappingSummary.xls -type PE && \
cp /hwfssz5/ST_MCHRI/COHORT/USER/zhongjixing/DPW/RNA-seq/result_2/process/GenomeMapping_HISAT/GenomeMappingSummary.xls /hwfssz5/ST_MCHRI/COHORT/USER/zhongjixing/DPW/RNA-seq/result_2/Analysis_Report/BGI_result/MapStat/GenomeMapping && \
cd /hwfssz5/ST_MCHRI/COHORT/USER/zhongjixing/DPW/RNA-seq/result_2/upload/IGV/bam/../../ && \
md5sum IGV/bam/*bam >> md5.txt && \
echo ==========end at : `date` ========== && \
echo Still_waters_run_deep 1>&2 && \
echo Still_waters_run_deep > /hwfssz5/ST_MCHRI/COHORT/USER/zhongjixing/DPW/RNA-seq/result_2/shell/GenomeMapping_HISAT/GenomeMap_stat.sh.sign
