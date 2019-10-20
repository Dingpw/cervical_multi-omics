#!/bin/bash
echo ==========start at : `date` ==========
perl /ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_RNAseq/RNA_RNAseq_2017a/Filter/filter_stat.pl -indir /hwfssz5/ST_MCHRI/COHORT/USER/zhongjixing/DPW/RNA-seq/result_2/process/Filter_SOAPnuke -output /hwfssz5/ST_MCHRI/COHORT/USER/zhongjixing/DPW/RNA-seq/result_2/process/Filter_SOAPnuke/FilterSummary.xls && \
cp /hwfssz5/ST_MCHRI/COHORT/USER/zhongjixing/DPW/RNA-seq/result_2/process/Filter_SOAPnuke/FilterSummary.xls /hwfssz5/ST_MCHRI/COHORT/USER/zhongjixing/DPW/RNA-seq/result_2/Analysis_Report/BGI_result/CleanData && \
cd /hwfssz5/ST_MCHRI/COHORT/USER/zhongjixing/DPW/RNA-seq/result_2/upload/CleanData/../ && \
md5sum CleanData/*fq.gz > md5.txt && \
echo ==========end at : `date` ========== && \
echo Still_waters_run_deep 1>&2 && \
echo Still_waters_run_deep > /hwfssz5/ST_MCHRI/COHORT/USER/zhongjixing/DPW/RNA-seq/result_2/shell/Filter_SOAPnuke/Filter_stat.sh.sign
