# scripts/04_haplotypecaller.sh
REF="data/reference/hg38.fa"
SAMPLE="sample1"

gatk HaplotypeCaller \
  -R $REF \
  -I results/bqsr/${SAMPLE}.bqsr.bam \
  -O results/variants/${SAMPLE}.g.vcf.gz \
  -ERC GVCF \
  --native-pair-hmm-threads 8