# Final NGS WES/WGS Dashboard

This is the final project version. It uses Flask, not Streamlit.

## Features

- Select FASTQ files already on server, or upload through browser
- Runs pipeline in background
- Saves all output under `jobs`
- Shows job status and pipeline log
- Downloads VCF, BAM, PDF report, logs, QC files

## Pipeline

1. FastQC
2. fastp trimming
3. BWA-MEM alignment to GRCh38
4. samtools sort and index
5. GATK MarkDuplicates
6. GATK BaseRecalibrator and ApplyBQSR
7. GATK HaplotypeCaller
8. PDF report generation

## Setup

```bash
cd final_dashboard
python3 -m pip install -r requirements.txt
cp config.example.json config.json
python3 app.py
```

Open:

```text
http://localhost:5000
```

For large FASTQ files, copy them into `input_fastq` first and use the server selection option.
