# scripts/01_align.sh
REF="data/reference/hg38.fa"
SAMPLE="sample1"

bwa mem -t 8 -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA" \
  $REF \
  data/raw/${SAMPLE}_R1.fastq.gz \
  data/raw/${SAMPLE}_R2.fastq.gz \
| samtools sort -o results/aligned/${SAMPLE}.bam

samtools index results/aligned/${SAMPLE}.bam