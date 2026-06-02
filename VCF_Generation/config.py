import os

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

REFERENCE = os.path.join(
    BASE_DIR,
    "reference/genome/GRCh38/genome.fa"
)

RAW_DATA = os.path.join(BASE_DIR, "data/raw")
TRIMMED = os.path.join(BASE_DIR, "data/trimmed")

OUTPUT = os.path.join(BASE_DIR, "output")

BAM_DIR = os.path.join(OUTPUT, "bam")
VCF_DIR = os.path.join(OUTPUT, "vcf")
REPORT_DIR = os.path.join(OUTPUT, "reports")

KNOWN_SNPS = os.path.join(
    BASE_DIR,
    "reference/known_sites/snps"
)

KNOWN_INDELS = os.path.join(
    BASE_DIR,
    "reference/known_sites/indels"
)