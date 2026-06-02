import argparse
import json
import subprocess
import os
from pathlib import Path

from report_generator import write_pdf_report


def write_status(job_dir, value):
    (job_dir / "status.txt").write_text(value + "\n", encoding="utf-8")
def run(command, log):
    log.write("\n$ " + " ".join(str(part) for part in command) + "\n")
    log.flush()

    env = os.environ.copy()
    env["PATH"] = "/home/puja_g/miniconda3/bin:" + env.get("PATH", "")

    subprocess.run(
        [str(part) for part in command],
        stdout=log,
        stderr=subprocess.STDOUT,
        check=True,
        env=env
    )

def pipe_bwa_to_sort(sample, r1, r2, sorted_bam, config, log):
    tools = config["tools"]
    threads = str(config.get("threads", 4))
    sort_memory = config.get("samtools_sort_memory", "1G")
    read_group = f"@RG\\tID:{sample}\\tSM:{sample}\\tPL:ILLUMINA\\tLB:{sample}\\tPU:{sample}"

    bwa_cmd = [
        tools["bwa"],
        "mem",
        "-t",
        threads,
        "-R",
        read_group,
        config["reference_fasta"],
        r1,
        r2,
    ]
    sort_cmd = [
        tools["samtools"],
        "sort",
        "-@",
        threads,
        "-m",
        sort_memory,
        "-o",
        sorted_bam,
        "-",
    ]

    log.write("\n$ " + " ".join(map(str, bwa_cmd)) + " | " + " ".join(map(str, sort_cmd)) + "\n")
    log.flush()
    env = os.environ.copy()
    env["PATH"] = "/home/puja_g/miniconda3/bin:" + env.get("PATH", "")
    bwa = subprocess.Popen(
        [str(x) for x in bwa_cmd],
        stdout=subprocess.PIPE,
        stderr=log,
        env=env)
    sort = subprocess.Popen(
        [str(x) for x in sort_cmd],
        stdin=bwa.stdout,
        stdout=log,
        stderr=log,
        env=env
        )

    if bwa.stdout:
        bwa.stdout.close()

    sort_code = sort.wait()
    bwa_code = bwa.wait()
    if bwa_code != 0:
        raise subprocess.CalledProcessError(bwa_code, bwa_cmd)
    if sort_code != 0:
        raise subprocess.CalledProcessError(sort_code, sort_cmd)


def gatk(config, tool, args):
    command = [config["tools"]["gatk"]]
    java_options = config.get("gatk_java_options")
    if java_options:
        command.extend(["--java-options", java_options])
    command.append(tool)
    command.extend(args)
    return command


def run_pipeline(sample, r1, r2, job_dir, config):
    tools = config["tools"]
    reference = config["reference_fasta"]
    known_sites = config.get("known_sites", [])

    qc_dir = job_dir / "qc"
    trimmed_dir = job_dir / "trimmed"
    bam_dir = job_dir / "bam"
    vcf_dir = job_dir / "vcf"
    report_dir = job_dir / "report"
    logs_dir = job_dir / "logs"
    for directory in (qc_dir, trimmed_dir, bam_dir, vcf_dir, report_dir, logs_dir):
        directory.mkdir(parents=True, exist_ok=True)

    log_path = logs_dir / "pipeline.log"
    fastp_json = qc_dir / f"{sample}.fastp.json"
    fastp_html = qc_dir / f"{sample}.fastp.html"
    trimmed_r1 = trimmed_dir / f"{sample}_R1.trimmed.fastq.gz"
    trimmed_r2 = trimmed_dir / f"{sample}_R2.trimmed.fastq.gz"
    sorted_bam = bam_dir / f"{sample}.sorted.bam"
    marked_bam = bam_dir / f"{sample}.marked.bam"
    metrics = bam_dir / f"{sample}.markdup.metrics.txt"
    recal_table = bam_dir / f"{sample}.recal.table"
    recal_bam = bam_dir / f"{sample}.recal.bam"
    final_bam = recal_bam if known_sites else marked_bam
    vcf = vcf_dir / f"{sample}.raw.vcf.gz"
    pdf = report_dir / f"{sample}_report.pdf"

    with log_path.open("w", encoding="utf-8") as log:
        write_status(job_dir, "RUNNING: QC and trimming")
        if "fastqc" in tools:
            run([tools["fastqc"], "-t", config.get("threads", 4), "-o", qc_dir, r1, r2], log)

        run(
            [
                tools["fastp"],
                "-i",
                r1,
                "-I",
                r2,
                "-o",
                trimmed_r1,
                "-O",
                trimmed_r2,
                "--thread",
                config.get("threads", 4),
                "--json",
                fastp_json,
                "--html",
                fastp_html,
            ],
            log,
        )

        write_status(job_dir, "RUNNING: Alignment and BAM sorting")
        pipe_bwa_to_sort(sample, trimmed_r1, trimmed_r2, sorted_bam, config, log)
        run([tools["samtools"], "index", sorted_bam], log)

        write_status(job_dir, "RUNNING: Mark duplicates")
        run(
            gatk(
                config,
                "MarkDuplicates",
                ["-I", sorted_bam, "-O", marked_bam, "-M", metrics, "--CREATE_INDEX", "true"],
            ),
            log,
        )

        if known_sites:
            write_status(job_dir, "RUNNING: BQSR")
            bqsr_args = ["-R", reference, "-I", marked_bam, "-O", recal_table]
            for known_site in known_sites:
                bqsr_args.extend(["--known-sites", known_site])
            run(gatk(config, "BaseRecalibrator", bqsr_args), log)
            run(
                gatk(
                    config,
                    "ApplyBQSR",
                    ["-R", reference, "-I", marked_bam, "--bqsr-recal-file", recal_table, "-O", recal_bam],
                ),
                log,
            )
        else:
            log.write("\nNo known_sites configured. Skipping BQSR.\n")

        write_status(job_dir, "RUNNING: Variant calling")
        run(gatk(config, "HaplotypeCaller", ["-R", reference, "-I", final_bam, "-O", vcf]), log)

        write_status(job_dir, "RUNNING: Report")
        write_pdf_report(sample, job_dir, vcf, final_bam, metrics, log_path, pdf)
        write_status(job_dir, "DONE")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    parser.add_argument("--r1", required=True)
    parser.add_argument("--r2", required=True)
    parser.add_argument("--job-dir", required=True)
    parser.add_argument("--config", required=True)
    args = parser.parse_args()

    job_dir = Path(args.job_dir)
    job_dir.mkdir(parents=True, exist_ok=True)

    try:
        with Path(args.config).open("r", encoding="utf-8") as handle:
            config = json.load(handle)
        run_pipeline(args.sample, Path(args.r1), Path(args.r2), job_dir, config)
    except Exception as error:
        with (job_dir / "error.txt").open("w", encoding="utf-8") as handle:
            handle.write(str(error) + "\n")
        write_status(job_dir, "FAILED")
        raise


if __name__ == "__main__":
    main()
