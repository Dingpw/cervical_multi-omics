#!/bin/bash

soap_mm_gz -a SAMPLE1_R1.fq.gz -b SAMPLE1_R2.fq.gz -D Human_rRNA_NCBI.fa.index -m 0 -x 1000 -v 5 -r 2 -p 3 -o SAMPLE1_rRNA.PESoap.gz -2 SAMPLE1_rRNA.PESoapSingle.gz && \
rRNAFilter.pl -fq SAMPLE1_R1.fq.gz,SAMPLE1_R2.fq.gz -soap SAMPLE1_rRNA.PESoap.gz,SAMPLE1_rRNA.PESoapSingle.gz -output SAMPLE1_rRNAremoved && \
fqcheck -r SAMPLE1_rRNAremoved_1.fq.gz -c 1.fqcheck && \
fqcheck -r SAMPLE1_rRNAremoved_2.fq.gz -c 2.fqcheck && \
fqcheck_distribute.pl 1.fqcheck 2.fqcheck -o SAMPLE1_rRNAremoved. && \
rm SAMPLE1_rRNA.PESoap.gz && \
rm SAMPLE1_rRNA.PESoapSingle.gz  && \
