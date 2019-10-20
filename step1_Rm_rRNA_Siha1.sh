#!/bin/bash

## Build index
2bwt-builder Human_rRNA_NCBI.fa && \

## Rm rRNA
soap_mm_gz -a SAMPLE1_R1.fq.gz -b SAMPLE_R2.fq.gz -D Human_rRNA_NCBI.fa.index -m 0 -x 1000 -v 5 -r 2 -p 3 -o SAMPLE1_rRNA.PESoap.gz -2 SAMPLE1_rRNA.PESoapSingle.gz && \
rRNAFilter.pl -fq SAMPLE1_R1.fq.gz,SAMPLE_R2.fq.gz -soap SAMPLE1_rRNA.PESoap.gz,SAMPLE1_rRNA.PESoapSingle.gz -output SAMPLE1_rRNAremoved && \
fqcheck -r SAMPLE_rRNAremoved_1.fq.gz -c 1.fqcheck && \
fqcheck -r SAMPLE_rRNAremoved_2.fq.gz -c 2.fqcheck && \
fqcheck_distribute.pl 1.fqcheck 2.fqcheck -o SAMPLE_rRNAremoved. && \
rm SAMPLE1_rRNA.PESoap.gz && \
rm SAMPLE1_rRNA.PESoapSingle.gz  && \

## Check Filter Stat
perl filter_stat.pl -indir /Filter_SOAPnuke -output /Filter_SOAPnuke/FilterSummary.xls && \
cp /Filter_SOAPnuke/FilterSummary.xls /Analysis_Report/BGI_result/CleanData && \
cd /upload/CleanData/../ && \
md5sum CleanData/*fq.gz > md5.txt && \

## Filter
tile=`perl findNtile2.pl -fq1 SAMPLE_rRNAremoved_1.fq.gz -fq2 SAMPLE_rRNAremoved_2.fq.gz  -seqType 0 ` && \
SOAPnuke filter -l 5 -q 0.5 -n 0.05 -Q 1 -5 0  -c 1 -1 SAMPLE_rRNAremoved_1.fq.gz -2 SAMPLE_rRNAremoved_2.fq.gz -f AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -r AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA $tile -o SAMPLE -C SAMPLE_1.fq.gz -D SAMPLE_2.fq.gz -R SAMPLE_1.rawdata.fq.gz -W SAMPLE_2.rawdata.fq.gz && \
mv SAMPLE/SAMPLE_1.fq.gz SAMPLE/SAMPLE_2.fq.gz CleanData && ln -fs CleanData/SAMPLE_1.fq.gz SAMPLE && ln -fs CleanData/SAMPLE_2.fq.gz SAMPLE && \
