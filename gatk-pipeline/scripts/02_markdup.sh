# scripts/02_markdup.sh
SAMPLE="sample1"

gatk MarkDuplicates \
  -I results/aligned/${SAMPLE}.bam \
  -O results/marked_dups/${SAMPLE}.markdup.bam \
  -M results/marked_dups/${SAMPLE}.metrics.txt

samtools index results/marked_dups/${SAMPLE}.markdup.bam