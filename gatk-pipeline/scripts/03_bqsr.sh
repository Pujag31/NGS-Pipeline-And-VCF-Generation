# scripts/03_bqsr.sh
REF="data/reference/hg38.fa"
KNOWN_SITES="data/reference/dbsnp_146.hg38.vcf.gz"
SAMPLE="sample1"

gatk BaseRecalibrator \
  -I results/marked_dups/${SAMPLE}.markdup.bam \
  -R $REF \
  --known-sites $KNOWN_SITES \
  -O results/bqsr/${SAMPLE}.recal.table

gatk ApplyBQSR \
  -I results/marked_dups/${SAMPLE}.markdup.bam \
  -R $REF \
  --bqsr-recal-file results/bqsr/${SAMPLE}.recal.table \
  -O results/bqsr/${SAMPLE}.bqsr.bam