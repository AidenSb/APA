#!/bin/bash
for f in  *_S1_L001_I1_001.fastq.gz;
do
sampname=${f:0:11};
id=${sampname}_Cellranger_out;
cellranger count --localcores 32 --localmem 200 --id=$id --fastqs=/data/APAproject/post_qual/data/Alexandra_Grubman/with_SRA/ --sample=$sampname --transcriptome=/home/aiden/data/refgenome/refdata-gex-GRCh38-2020-A/;
done


