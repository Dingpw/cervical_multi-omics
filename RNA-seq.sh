#!/bin/bash

## step1_Build_index
bwt-builder process/Rm_rRNA/index/Human_rRNA_NCBI.fa && \

## step2_Rm_rRNA
soap_mm_gz -a SAMPLE1_R1.fq.gz -b SAMPLE_R2.fq.gz -D process/Rm_rRNA/index/Human_rRNA_NCBI.fa.index -m 0 -x 1000 -v 5 -r 2 -p 3 -o process/Rm_rRNA/SAMPLE_rRNA.PESoap.gz -2 process/Rm_rRNA/SAMPLE_rRNA.PESoapSingle.gz && \
rRNAFilter.pl -fq SAMPLE1_R1.fq.gz,SAMPLE_R2.fq.gz -soap process/Rm_rRNA/SAMPLE_rRNA.PESoap.gz,process/Rm_rRNA/SAMPLE_rRNA.PESoapSingle.gz -output process/Rm_rRNA/SAMPLE_rRNAremoved && \
fqcheck -r process/Rm_rRNA/SAMPLE_rRNAremoved_1.fq.gz -c process/Rm_rRNA/1.fqcheck && \
fqcheck -r process/Rm_rRNA/SAMPLE_rRNAremoved_2.fq.gz -c process/Rm_rRNA/2.fqcheck && \
fqcheck_distribute.pl process/Rm_rRNA/1.fqcheck process/Rm_rRNA/2.fqcheck -o process/Rm_rRNA/SAMPLE_rRNAremoved. && \
rm process/Rm_rRNA/SAMPLE_rRNA.PESoap.gz && \
rm process/Rm_rRNA/SAMPLE_rRNA.PESoapSingle.gz  && \

## step3_Filter
tile=`perl findNtile2.pl -fq1 process/Rm_rRNA/SAMPLE_rRNAremoved_1.fq.gz -fq2 process/Rm_rRNA/SAMPLE_rRNAremoved_2.fq.gz  -seqType 0 ` && \
SOAPnuke filter -l 5 -q 0.5 -n 0.05 -Q 1 -5 0  -c 1 -1 process/Rm_rRNA/SAMPLE_rRNAremoved_1.fq.gz -2 process/Rm_rRNA/SAMPLE_rRNAremoved_2.fq.gz -f AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -r AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA $tile -o process/Filter_SOAPnuke -C SAMPLE_1.fq.gz -D SAMPLE_2.fq.gz -R SAMPLE_1.rawdata.fq.gz -W SAMPLE_2.rawdata.fq.gz && \
mv process/Filter_SOAPnuke/SAMPLE_1.fq.gz process/Filter_SOAPnuke/SAMPLE_2.fq.gz CleanData && ln -fs CleanData/SAMPLE_1.fq.gz process/Filter_SOAPnuke && ln -fs CleanData/SAMPLE_2.fq.gz process/Filter_SOAPnuke && \

## step4_Check_Filter_stat
perl filter_stat.pl -indir process/Filter_SOAPnuke -output process/Filter_SOAPnuke/FilterSummary.xls && \
cp process/Filter_SOAPnuke/FilterSummary.xls /Analysis_Report/CleanData && \
cd CleanData/../ && \
md5sum CleanData/*fq.gz > md5.txt && \

## step5_Genomapping
# HISAT alignment
cd process/GenomeMapping_HISAT/SAMPLE && \
hisat2 --phred64 --sensitive --no-discordant --no-mixed -I 1 -X 1000 -x Database/hg19/GenomeHisat2Index/chrALL -1 CleanData/SAMPLE_1.fq.gz -2 CleanData/SAMPLE_2.fq.gz 2>SAMPLE.Map2GenomeStat.xls | samtools view -b -S -o SAMPLE.bam - && \

# AddRG & Reorder & Sort BAM for downstream analysis
if [ ! -d java_tmp ];then mkdir -p java_tmp;fi && \
java -Xmx4G -Djava.io.tmpdir=java_tmp -jar picard.jar AddOrReplaceReadGroups I=SAMPLE.bam O=SAMPLE.AddRG.bam RGID=SAMPLE RGLB=SAMPLE_library RGPL=illumina RGPU=machine RGSM=SAMPLE VALIDATION_STRINGENCY=SILENT && \
java -Xmx4G -Djava.io.tmpdir=java_tmp -jar picard.jar ReorderSam I=SAMPLE.AddRG.bam O=SAMPLE.AddRG.Reorder.bam R=Database/hg19/GenomeGatkIndex/chrALL.sort.fa VALIDATION_STRINGENCY=SILENT && \
java -Xmx4G -Djava.io.tmpdir=java_tmp -jar picard.jar SortSam I=SAMPLE.AddRG.Reorder.bam O=SAMPLE.AddRG.Reorder.Sort.bam SO=coordinate VALIDATION_STRINGENCY=SILENT && \
samtools index SAMPLE.AddRG.Reorder.Sort.bam && \
mv SAMPLE.AddRG.Reorder.Sort.bam SAMPLE.AddRG.Reorder.Sort.bam.bai IGV/bam && \
rm -rf java_tmp && \
rm SAMPLE.bam SAMPLE.AddRG.bam SAMPLE.AddRG.Reorder.bam && \
cp SAMPLE.Map2GenomeStat.xls Analysis_Report/MapStat/GenomeMapping  && \

## Check alignment stat
perl MapStat.pl -indir process/GenomeMapping_HISAT -output process/GenomeMapping_HISAT/GenomeMappingSummary.xls -type PE && \
cp process/GenomeMapping_HISAT/GenomeMappingSummary.xls Analysis_Report/MapStat/GenomeMapping && \
cd IGV/bam/../../ && \
md5sum IGV/bam/*bam >> md5.txt && \

## step6_GeneExpression
## build Index
if [ ! -d process/GeneExp_RSEM/rsem-build ];then mkdir -p process/GeneExp_RSEM/rsem-build;fi && \
cat Database/hg19/GeneBowtie2Index/refMrna.fa >process/GeneExp_RSEM/rsem-build/refMrna.fa && \
cat Database/hg19/Annotation_kegg76/refMrna.fa.gene2mark >process/GeneExp_RSEM/gene2tr.txt && \
export LD_LIBRARY_PATH=/GeneExp/../software/RNA_lib:$LD_LIBRARY_PATH && \
rsem-prepare-reference process/GeneExp_RSEM/rsem-build/refMrna.fa process/GeneExp_RSEM/rsem-build/refMrna.fa --bowtie2 --bowtie2-path BOWTIE_PATH --transcript-to-gene-map process/GeneExp_RSEM/gene2tr.txt && \
perl fastaDeal.pl -attr id:len process/GeneExp_RSEM/rsem-build/refMrna.fa >process/GeneExp_RSEM/TranscriptLength.txt && \
awk '{print $1"\t1\t"$2}' process/GeneExp_RSEM/TranscriptLength.txt >process/GeneExp_RSEM/TranscriptLength.bed  && \

## Gene Expression
export LD_LIBRARY_PATH=/GeneExp/../software/RNA_lib:$LD_LIBRARY_PATH && \
bowtie2 -q --phred64 --sensitive --dpad 0 --gbar 99999999 --mp 1,1 --np 1 --score-min L,0,-0.1 -I 1 -X 1000 --no-mixed --no-discordant  -p 1 -k 200 -x process/GeneExp_RSEM/rsem-build/refMrna.fa -1 CleanData/SAMPLE_1.fq.gz -2 CleanData/SAMPLE_2.fq.gz 2>process/GeneExp_RSEM/SAMPLE.Map2GeneStat.xls | samtools view -S -b -o process/GeneExp_RSEM/SAMPLE.bam - && \
rsem-calculate-expression --paired-end -p 8 --bam process/GeneExp_RSEM/SAMPLE.bam process/GeneExp_RSEM/rsem-build/refMrna.fa process/GeneExp_RSEM/SAMPLE && \
awk '{if($7!=0.00)print $1"\t"$2"\t"$3"\t"$5"\t"$7}' process/GeneExp_RSEM/SAMPLE.genes.results |grep -v '^ERCC' >process/GeneExp_RSEM/SAMPLE.gene.fpkm.xls && \
awk '{if($7!=0.00)print $1"\t"$2"\t"$3"\t"$5"\t"$7}' process/GeneExp_RSEM/SAMPLE.isoforms.results |grep -v '^ERCC' >process/GeneExp_RSEM/SAMPLE.transcript.fpkm.xls && \
cp process/GeneExp_RSEM/SAMPLE.gene.fpkm.xls process/GeneExp_RSEM/SAMPLE.transcript.fpkm.xls Analysis_Report/Quantify/GeneExpression/GeneExpression  && \

## Check mapping stat
BowtieMapStat.pl -bam process/GeneExp_RSEM/SAMPLE.bam -key process/GeneExp_RSEM/SAMPLE.Bowtie2Gene -seqType PE -samtools samtools -gene2tr process/GeneExp_RSEM/gene2tr.txt && \
cp process/GeneExp_RSEM/SAMPLE.Bowtie2Gene.MapReadsStat.xls Analysis_Report/MapStat/GeneMapping  && \

## Reads Coverage
samtools depth -b process/GeneExp_RSEM/TranscriptLength.bed process/GeneExp_RSEM/SAMPLE.transcript.sorted.bam >process/GeneExp_RSEM/SAMPLE.depth.txt && \
perl cover_stat.pl process/GeneExp_RSEM/SAMPLE.depth.txt process/GeneExp_RSEM/TranscriptLength.txt SAMPLE process/GeneExp_RSEM/SAMPLE && \
cp process/GeneExp_RSEM/SAMPLE/SAMPLE.ReadsCoverage.pdf process/GeneExp_RSEM/SAMPLE/SAMPLE.ReadsCoverage.png Analysis_Report/Quantify/GeneExpression/ReadsCoverage  && \

## Reads Random
perl ReadsRandomInGene_forBowtie_ggplot2.pl -len process/GeneExp_RSEM/TranscriptLength.txt -bam process/GeneExp_RSEM/SAMPLE.bam -seqType PE -prefix process/GeneExp_RSEM/ && \
cp process/GeneExp_RSEM/SAMPLE.ReadsRandom.pdf process/GeneExp_RSEM/SAMPLE.ReadsRandom.png Analysis_Report/Quantify/GeneExpression/ReadsRandom  && \

## SeqSaturation
SeqSaturation_Bowtie.pl -bam process/GeneExp_RSEM/SAMPLE.transcript.bam -seqType PE -gene2tr process/GeneExp_RSEM/gene2tr.txt -prefix process/GeneExp_RSEM -cutoff 1 && \
mv process/GeneExp_RSEM/SAMPLE.SeqSaturation* Analysis_Report/Quantify/GeneExpression/SeqSaturation && \

## statistic
# all samples fpkm statistic
perl AllGeneStat.pl process/GeneExp_RSEM process/GeneExp_RSEM/AllSamples.GeneExpression.FPKM.xls && \
perl AllTranscriptStat.pl process/GeneExp_RSEM process/GeneExp_RSEM/AllSamples.TranscriptExpression.FPKM.xls && \
cp process/GeneExp_RSEM/AllSamples.GeneExpression.FPKM.xls Analysis_Report/Quantify/GeneExpression/GeneExpression && \
perl ExpStat.pl -indir process/GeneExp_RSEM -output process/GeneExp_RSEM/GeneExpressionSummary.xls && \
cp process/GeneExp_RSEM/GeneExpressionSummary.xls Analysis_Report/Quantify/GeneExpression && \
# correlation heatmap
perl drawCorrelationHeatmap.pl -indir Analysis_Report/Quantify/GeneExpression/GeneExpression -outdir Analysis_Report/Quantify/GeneExpression/CorrelationHeatmap && \
mv Analysis_Report/Quantify/GeneExpression/CorrelationHeatmap/correlation-heatmap.R process/GeneExp_RSEM && \
# hcluster tree
perl drawHclustTree.pl -indir Analysis_Report/Quantify/GeneExpression/GeneExpression -outdir Analysis_Report/Quantify/GeneExpression/HclusterTree && \
perl drawbox_density-plot.pl -i Analysis_Report/Quantify/GeneExpression/GeneExpression -o Analysis_Report/Quantify/GeneExpression/GeneExpression && \
perl drawstacked_bar.pl -indir Analysis_Report/Quantify/GeneExpression/GeneExpression -outdir Analysis_Report/Quantify/GeneExpression/GeneExpression && \
