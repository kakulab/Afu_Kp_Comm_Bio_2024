#!/bin/bash

WKDIR=/data/Michi/RNAseq_intro


# Mapping of all .bam files with NGM 

for i in $WKDIR/*.bam
do
ngm -q $i -r $WKDIR/*.fasta -o $i.sam -t 6  # -t 6 is optional - means 6 threads of the processor are used, if you don't know what to do, remove it
done

# convert .sam to .bam using samtools (only required if mapping output is .sam)

for i in $WKDIR/*.sam
do
samtools view -bS $i > $i.bam
done

# sort -bam files using samtools

for i in $WKDIR/*.sam.bam
do
samtools sort $i $i.sorted
done

# flagstat analysis

for i in $WKDIR/*.sorted.bam
do
samtools flagstat $i > $i.flagstat.txt
done


# read count with htseq-count

for i in $WKDIR/*.sorted.bam
do
htseq-count -f bam -s no -t gene -i ID $i $WKDIR/*.gff > $i.count.txt
done

# clear count files for flags (removing rRNA reads is already done for during counting by selecting -t gene)

for i in $WKDIR/*.count.txt
do
head -n -5 $i > $i.crop.txt
done

