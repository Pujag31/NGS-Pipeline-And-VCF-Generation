# scripts/05_genotypegvcf.sh
REF="data/reference/hg38.fa"

# Combine GVCFs (if multiple samples)
gatk CombineGVCFs \
  -R $REF \
  -V results/variants/sample1.g.vcf.gz \
  -V results/variants/sample2.g.vcf.gz \
  -O results/variants/cohort.g.vcf.gz

# Genotype
gatk GenotypeGVCFs \
  -R $REF \
  -V results/variants/cohort.g.vcf.gz \
  -O results/variants/cohort.vcf.gzss