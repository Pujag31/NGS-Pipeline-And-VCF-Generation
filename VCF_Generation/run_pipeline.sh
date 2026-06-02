#!/bin/bash

FASTQ_FILE=$1

mkdir -p output

echo "Running FastQC..."

fastqc $FASTQ_FILE -o output/

echo "Running Alignment..."

bwa mem reference/ref.fa $FASTQ_FILE > output/aligned.sam

echo "Converting SAM to BAM..."

samtools view -S -b output/aligned.sam > output/aligned.bam

echo "Sorting BAM..."

samtools sort output/aligned.bam -o output/sorted.bam

echo "Indexing BAM..."

samtools index output/sorted.bam

echo "Calling Variants..."

bcftools mpileup -f reference/ref.fa output/sorted.bam | \
bcftools call -mv -Ov -o output/variants.vcf

echo "Pipeline Completed Successfully!"