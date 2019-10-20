#!/bin/bash
echo ==========start at : `date` ==========
/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_RNAseq/RNA_RNAseq_2017a/Rm_rRNA/bin/Soap/2bwt-builder /hwfssz5/ST_MCHRI/COHORT/USER/zhongjixing/DPW/RNA-seq/result_2/process/Rm_rRNA/index/Human_rRNA_NCBI.fa && \
echo ==========end at : `date` ========== && \
echo Still_waters_run_deep 1>&2 && \
echo Still_waters_run_deep > /hwfssz5/ST_MCHRI/COHORT/USER/zhongjixing/DPW/RNA-seq/result_2/shell/Rm_rRNA/builder.sh.sign
