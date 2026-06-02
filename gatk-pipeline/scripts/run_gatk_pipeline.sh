#!/bin/bash
# ============================================================
#  GATK VARIANT CALLING PIPELINE — MASTER BATCH SCRIPT
#  Supports multiple samples automatically
#  Usage: bash run_gatk_pipeline.sh
# ============================================================

set -euo pipefail  # Exit on error, undefined vars, pipe failures

# ─────────────────────────────────────────
#  CONFIGURATION — Edit these paths
# ─────────────────────────────────────────
REF="data/reference/hg38.fa"
KNOWN_SITES="data/reference/dbsnp_146.hg38.vcf.gz"
RAW_DIR="data/raw"
THREADS=8

# Output directories
ALIGNED="results/aligned"
MARKDUP="results/marked_dups"
BQSR="results/bqsr"
VARIANTS="results/variants"
LOGS="results/logs"

# ─────────────────────────────────────────
#  SETUP
# ─────────────────────────────────────────
mkdir -p $ALIGNED $MARKDUP $BQSR $VARIANTS $LOGS

# Detect all samples from paired FASTQ files (e.g., sample1_R1.fastq.gz)
SAMPLES=($(ls ${RAW_DIR}/*_R1.fastq.gz 2>/dev/null | xargs -n1 basename | sed 's/_R1.fastq.gz//'))

if [ ${#SAMPLES[@]} -eq 0 ]; then
  echo "❌ No FASTQ files found in ${RAW_DIR}/"
  echo "   Expected format: sampleName_R1.fastq.gz / sampleName_R2.fastq.gz"
  exit 1
fi

echo "============================================"
echo "  GATK BATCH PIPELINE"
echo "  Samples detected: ${#SAMPLES[@]}"
echo "  $(printf '  - %s\n' "${SAMPLES[@]}")"
echo "============================================"

# ─────────────────────────────────────────
#  STEP 1: Index Reference (once)
# ─────────────────────────────────────────
log_step() { echo -e "\n[$(date '+%H:%M:%S')] ▶ $1"; }

log_step "Checking reference index..."
if [ ! -f "${REF}.bwt" ]; then
  echo "  Indexing reference genome..."
  bwa index $REF
  samtools faidx $REF
  gatk CreateSequenceDictionary -R $REF
else
  echo "  Reference already indexed. Skipping."
fi

# ─────────────────────────────────────────
#  PER-SAMPLE PROCESSING
# ─────────────────────────────────────────
GVCF_INPUTS=""  # For joint genotyping later

for SAMPLE in "${SAMPLES[@]}"; do
  echo ""
  echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
  echo "  Processing sample: ${SAMPLE}"
  echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

  R1="${RAW_DIR}/${SAMPLE}_R1.fastq.gz"
  R2="${RAW_DIR}/${SAMPLE}_R2.fastq.gz"
  LOG="${LOGS}/${SAMPLE}.log"

  # ── STEP 2: Align with BWA-MEM ──
  log_step "[${SAMPLE}] Step 1/4 — Aligning reads..."
  if [ ! -f "${ALIGNED}/${SAMPLE}.bam" ]; then
    bwa mem -t $THREADS \
      -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tLB:lib1\tPL:ILLUMINA" \
      $REF $R1 $R2 2>>"$LOG" \
    | samtools sort -@ $THREADS -o "${ALIGNED}/${SAMPLE}.bam"
    samtools index "${ALIGNED}/${SAMPLE}.bam"
    echo "  ✔ Alignment done → ${ALIGNED}/${SAMPLE}.bam"
  else
    echo "  ⏭ Already aligned. Skipping."
  fi

  # ── STEP 3: Mark Duplicates ──
  log_step "[${SAMPLE}] Step 2/4 — Marking duplicates..."
  if [ ! -f "${MARKDUP}/${SAMPLE}.markdup.bam" ]; then
    gatk MarkDuplicates \
      -I "${ALIGNED}/${SAMPLE}.bam" \
      -O "${MARKDUP}/${SAMPLE}.markdup.bam" \
      -M "${MARKDUP}/${SAMPLE}.metrics.txt" \
      --TMP_DIR /tmp 2>>"$LOG"
    samtools index "${MARKDUP}/${SAMPLE}.markdup.bam"
    echo "  ✔ Duplicates marked → ${MARKDUP}/${SAMPLE}.markdup.bam"
  else
    echo "  ⏭ Already done. Skipping."
  fi

  # ── STEP 4: BQSR ──
  log_step "[${SAMPLE}] Step 3/4 — Base Quality Score Recalibration..."
  if [ ! -f "${BQSR}/${SAMPLE}.bqsr.bam" ]; then
    gatk BaseRecalibrator \
      -I "${MARKDUP}/${SAMPLE}.markdup.bam" \
      -R $REF \
      --known-sites $KNOWN_SITES \
      -O "${BQSR}/${SAMPLE}.recal.table" 2>>"$LOG"

    gatk ApplyBQSR \
      -I "${MARKDUP}/${SAMPLE}.markdup.bam" \
      -R $REF \
      --bqsr-recal-file "${BQSR}/${SAMPLE}.recal.table" \
      -O "${BQSR}/${SAMPLE}.bqsr.bam" 2>>"$LOG"

    samtools index "${BQSR}/${SAMPLE}.bqsr.bam"
    echo "  ✔ BQSR done → ${BQSR}/${SAMPLE}.bqsr.bam"
  else
    echo "  ⏭ Already done. Skipping."
  fi

  # ── STEP 5: HaplotypeCaller ──
  log_step "[${SAMPLE}] Step 4/4 — HaplotypeCaller (gVCF)..."
  if [ ! -f "${VARIANTS}/${SAMPLE}.g.vcf.gz" ]; then
    gatk HaplotypeCaller \
      -R $REF \
      -I "${BQSR}/${SAMPLE}.bqsr.bam" \
      -O "${VARIANTS}/${SAMPLE}.g.vcf.gz" \
      -ERC GVCF \
      --native-pair-hmm-threads $THREADS 2>>"$LOG"
    echo "  ✔ gVCF done → ${VARIANTS}/${SAMPLE}.g.vcf.gz"
  else
    echo "  ⏭ Already done. Skipping."
  fi

  GVCF_INPUTS+=" -V ${VARIANTS}/${SAMPLE}.g.vcf.gz"
done

# ─────────────────────────────────────────
#  STEP 6: Joint Genotyping (all samples)
# ─────────────────────────────────────────
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  Joint Genotyping (all samples)"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

log_step "Combining GVCFs..."
if [ ! -f "${VARIANTS}/cohort.g.vcf.gz" ]; then
  gatk CombineGVCFs \
    -R $REF \
    $GVCF_INPUTS \
    -O "${VARIANTS}/cohort.g.vcf.gz" 2>>"${LOGS}/joint.log"
  echo "  ✔ Combined gVCF → ${VARIANTS}/cohort.g.vcf.gz"
fi

log_step "Running GenotypeGVCFs..."
if [ ! -f "${VARIANTS}/cohort.vcf.gz" ]; then
  gatk GenotypeGVCFs \
    -R $REF \
    -V "${VARIANTS}/cohort.g.vcf.gz" \
    -O "${VARIANTS}/cohort.vcf.gz" 2>>"${LOGS}/joint.log"
  echo "  ✔ Final VCF → ${VARIANTS}/cohort.vcf.gz"
fi

# ─────────────────────────────────────────
#  STEP 7: Variant Filtering (VQSR / hard filter)
# ─────────────────────────────────────────
log_step "Applying hard filters to SNPs and INDELs..."

# SNPs
gatk SelectVariants -R $REF -V "${VARIANTS}/cohort.vcf.gz" \
  --select-type-to-include SNP \
  -O "${VARIANTS}/cohort.snps.vcf.gz" 2>>"${LOGS}/joint.log"

gatk VariantFiltration -R $REF \
  -V "${VARIANTS}/cohort.snps.vcf.gz" \
  --filter-expression "QD < 2.0"    --filter-name "QD2" \
  --filter-expression "FS > 60.0"   --filter-name "FS60" \
  --filter-expression "MQ < 40.0"   --filter-name "MQ40" \
  --filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
  --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
  -O "${VARIANTS}/cohort.snps.filtered.vcf.gz" 2>>"${LOGS}/joint.log"

# INDELs
gatk SelectVariants -R $REF -V "${VARIANTS}/cohort.vcf.gz" \
  --select-type-to-include INDEL \
  -O "${VARIANTS}/cohort.indels.vcf.gz" 2>>"${LOGS}/joint.log"

gatk VariantFiltration -R $REF \
  -V "${VARIANTS}/cohort.indels.vcf.gz" \
  --filter-expression "QD < 2.0"   --filter-name "QD2" \
  --filter-expression "FS > 200.0" --filter-name "FS200" \
  --filter-expression "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
  -O "${VARIANTS}/cohort.indels.filtered.vcf.gz" 2>>"${LOGS}/joint.log"

echo ""
echo "============================================"
echo "  ✅ PIPELINE COMPLETE"
echo "  SNPs    → ${VARIANTS}/cohort.snps.filtered.vcf.gz"
echo "  INDELs  → ${VARIANTS}/cohort.indels.filtered.vcf.gz"
echo "  Logs    → ${LOGS}/"
echo "============================================"
