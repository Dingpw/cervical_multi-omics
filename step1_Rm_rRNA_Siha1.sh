#!/bin/bash
echo ==========start at : `date` ==========
/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_RNAseq/RNA_RNAseq_2017a/Rm_rRNA/bin/Soap/soap_mm_gz -a /hwfssz5/ST_MCHRI/COHORT/USER/zhongjixing/DPW/cervical/HIC-RNA/ANYWSH170024_PM-YWSH170024-01_2019-04-11/Rawdata/Siha1/Siha1_R1.fq.gz -b /hwfssz5/ST_MCHRI/COHORT/USER/zhongjixing/DPW/cervical/HIC-RNA/ANYWSH170024_PM-YWSH170024-01_2019-04-11/Rawdata/Siha1/Siha1_R2.fq.gz -D /hwfssz5/ST_MCHRI/COHORT/USER/zhongjixing/DPW/RNA-seq/result_2/process/Rm_rRNA/index/Human_rRNA_NCBI.fa.index -m 0 -x 1000 -v 5 -r 2 -p 3 -o /hwfssz5/ST_MCHRI/COHORT/USER/zhongjixing/DPW/RNA-seq/result_2/process/Rm_rRNA/Siha1/Siha1/Siha1_rRNA.PESoap.gz -2 /hwfssz5/ST_MCHRI/COHORT/USER/zhongjixing/DPW/RNA-seq/result_2/process/Rm_rRNA/Siha1/Siha1/Siha1_rRNA.PESoapSingle.gz && \
/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_RNAseq/RNA_RNAseq_2017a/Rm_rRNA/bin/rRNAFilter.pl -fq /hwfssz5/ST_MCHRI/COHORT/USER/zhongjixing/DPW/cervical/HIC-RNA/ANYWSH170024_PM-YWSH170024-01_2019-04-11/Rawdata/Siha1/Siha1_R1.fq.gz,/hwfssz5/ST_MCHRI/COHORT/USER/zhongjixing/DPW/cervical/HIC-RNA/ANYWSH170024_PM-YWSH170024-01_2019-04-11/Rawdata/Siha1/Siha1_R2.fq.gz -soap /hwfssz5/ST_MCHRI/COHORT/USER/zhongjixing/DPW/RNA-seq/result_2/process/Rm_rRNA/Siha1/Siha1/Siha1_rRNA.PESoap.gz,/hwfssz5/ST_MCHRI/COHORT/USER/zhongjixing/DPW/RNA-seq/result_2/process/Rm_rRNA/Siha1/Siha1/Siha1_rRNA.PESoapSingle.gz -output /hwfssz5/ST_MCHRI/COHORT/USER/zhongjixing/DPW/RNA-seq/result_2/process/Rm_rRNA/Siha1/Siha1/Siha1_rRNAremoved && \
/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_RNAseq/RNA_RNAseq_2017a/Rm_rRNA/bin/fqcheck -r /hwfssz5/ST_MCHRI/COHORT/USER/zhongjixing/DPW/RNA-seq/result_2/process/Rm_rRNA/Siha1/Siha1/Siha1_rRNAremoved_1.fq.gz -c /hwfssz5/ST_MCHRI/COHORT/USER/zhongjixing/DPW/RNA-seq/result_2/process/Rm_rRNA/Siha1/Siha1/1.fqcheck && \
/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_RNAseq/RNA_RNAseq_2017a/Rm_rRNA/bin/fqcheck -r /hwfssz5/ST_MCHRI/COHORT/USER/zhongjixing/DPW/RNA-seq/result_2/process/Rm_rRNA/Siha1/Siha1/Siha1_rRNAremoved_2.fq.gz -c /hwfssz5/ST_MCHRI/COHORT/USER/zhongjixing/DPW/RNA-seq/result_2/process/Rm_rRNA/Siha1/Siha1/2.fqcheck && \
/ifs4/BC_PUB/biosoft/pipeline/RNA/RNA_RNAseq/RNA_RNAseq_2017a/Rm_rRNA/bin/fqcheck_distribute.pl /hwfssz5/ST_MCHRI/COHORT/USER/zhongjixing/DPW/RNA-seq/result_2/process/Rm_rRNA/Siha1/Siha1/1.fqcheck /hwfssz5/ST_MCHRI/COHORT/USER/zhongjixing/DPW/RNA-seq/result_2/process/Rm_rRNA/Siha1/Siha1/2.fqcheck -o /hwfssz5/ST_MCHRI/COHORT/USER/zhongjixing/DPW/RNA-seq/result_2/process/Rm_rRNA/Siha1/Siha1/Siha1_rRNAremoved. && \
rm /hwfssz5/ST_MCHRI/COHORT/USER/zhongjixing/DPW/RNA-seq/result_2/process/Rm_rRNA/Siha1/Siha1/Siha1_rRNA.PESoap.gz && \
rm /hwfssz5/ST_MCHRI/COHORT/USER/zhongjixing/DPW/RNA-seq/result_2/process/Rm_rRNA/Siha1/Siha1/Siha1_rRNA.PESoapSingle.gz  && \
echo ==========end at : `date` ========== && \
echo Still_waters_run_deep 1>&2 && \
echo Still_waters_run_deep > /hwfssz5/ST_MCHRI/COHORT/USER/zhongjixing/DPW/RNA-seq/result_2/shell/Rm_rRNA/Siha1/Rm_rRNA_Siha1.sh.sign