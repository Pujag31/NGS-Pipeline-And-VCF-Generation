import argparse
import json
import subprocess
from pathlib import Path


def run(command, log):
    log.write("\n$ " + " ".join(map(str, command)) + "\n")
    log.flush()
    subprocess.run(list(map(str, command)), stdout=log, stderr=subprocess.STDOUT, check=True)


def run_bwa_sort(sample, r1, r2, out_bam, config, log):
    tools = config["tools"]
    ref = config["reference_fasta"]
    threads = str(config.get("threads", 4))
    sort_mem = config.get("samtools_sort_memory", "1G")

    rg = f"@RG\\tID:{sample}\\tSM:{sample}\\tPL:ILLUMINA\\tLB:{sample}\\tPU:{sample}"

    bwa_cmd = [tools["bwa"], "mem", "-t", threads, "-R", rg, ref, r1, r2]
    sort_cmd = [tools["samtools"], "sort", "-@", threads, "-m", sort_mem, "-o", out_bam, "-"]

    bwa = subprocess.Popen(list(map(str, bwa_cmd)), stdout=subprocess.PIPE, stderr=log)
    sort = subprocess.Popen(list(map(str, sort_cmd)), stdin=bwa.stdout, stdout=log, stderr=log)

    if bwa.stdout:
        bwa.stdout.close()

    sort_code = sort.wait()
    bwa_code = bwa.wait()

    if bwa_code != 0:
        raise subprocess.CalledProcessError(bwa_code, bwa_cmd)
    if sort_code != 0:
        raise subprocess.CalledProcessError(sort_code, sort_cmd)


def gatk(config, tool, args):
    cmd = [config["tools"]["gatk"]]
    java_options = config.get("gatk_java_options")
    if java_options:
        cmd += ["--java-options", java_options]
    cmd += [tool] + args
    return cmd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    parser.add_argument("--r1", required=True)
    parser.add_argument("--r2", required=True)
    parser.add_argument("--config", default="config.json")
    parser.add_argument("--outdir", default="results")
    args = parser.parse_args()

    sample = args.sample
    outdir = Path(args.outdir) / sample
    outdir.mkdir(parents=True, exist_ok=True)

    with open(args.config, "r", encoding="utf-8") as f:
        config = json.load(f)

    ref = config["reference_fasta"]
    known_sites = config.get("known_sites", [])
    tools = config["tools"]

    sorted_bam = outdir / f"{sample}.sorted.bam"
    marked_bam = outdir / f"{sample}.marked.bam"
    metrics = outdir / f"{sample}.markdup.metrics.txt"
    recal_table = outdir / f"{sample}.recal.table"
    recal_bam = outdir / f"{sample}.recal.bam"
    vcf = outdir / f"{sample}.raw.vcf.gz"
    log_path = outdir / "pipeline.log"

    with open(log_path, "w", encoding="utf-8") as log:
        run_bwa_sort(sample, args.r1, args.r2, sorted_bam, config, log)
        run([tools["samtools"], "index", sorted_bam], log)

        run(gatk(config, "MarkDuplicates", [
            "-I", sorted_bam,
            "-O", marked_bam,
            "-M", metrics,
            "--CREATE_INDEX", "true"
        ]), log)

        if known_sites:
            bqsr_cmd = [
                "-R", ref,
                "-I", marked_bam,
                "-O", recal_table
            ]
            for site in known_sites:
                bqsr_cmd += ["--known-sites", site]

            run(gatk(config, "BaseRecalibrator", bqsr_cmd), log)

            run(gatk(config, "ApplyBQSR", [
                "-R", ref,
                "-I", marked_bam,
                "--bqsr-recal-file", recal_table,
                "-O", recal_bam
            ]), log)

            caller_input = recal_bam
        else:
            log.write("\nNo known_sites found. Skipping BQSR.\n")
            caller_input = marked_bam

        run(gatk(config, "HaplotypeCaller", [
            "-R", ref,
            "-I", caller_input,
            "-O", vcf
        ]), log)

    print("DONE")
    print("VCF:", vcf)
    print("LOG:", log_path)


if __name__ == "__main__":
    main()
