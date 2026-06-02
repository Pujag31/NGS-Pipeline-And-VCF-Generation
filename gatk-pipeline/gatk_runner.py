import argparse
import json
import subprocess
from pathlib import Path


def run_command(command, log_file):
    log_file.write("\n$ " + " ".join(str(part) for part in command) + "\n")
    log_file.flush()
    subprocess.run(
        [str(part) for part in command],
        stdout=log_file,
        stderr=subprocess.STDOUT,
        check=True,
    )


def run_bwa_and_sort(sample_name, read1, read2, sorted_bam, config, log_file):
    tools = config["tools"]
    reference = config["reference_fasta"]
    threads = str(config.get("threads", 2))
    sort_memory = config.get("samtools_sort_memory", "1G")
    read_group = (
        f"@RG\\tID:{sample_name}\\tSM:{sample_name}\\tPL:ILLUMINA"
        f"\\tLB:{sample_name}\\tPU:{sample_name}"
    )

    bwa_command = [
        tools["bwa"],
        "mem",
        "-t",
        threads,
        "-R",
        read_group,
        reference,
        read1,
        read2,
    ]
    sort_command = [
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

    log_file.write(
        "\n$ "
        + " ".join(str(part) for part in bwa_command)
        + " | "
        + " ".join(str(part) for part in sort_command)
        + "\n"
    )
    log_file.flush()

    bwa_process = subprocess.Popen(
        [str(part) for part in bwa_command],
        stdout=subprocess.PIPE,
        stderr=log_file,
    )
    sort_process = subprocess.Popen(
        [str(part) for part in sort_command],
        stdin=bwa_process.stdout,
        stdout=log_file,
        stderr=log_file,
    )

    if bwa_process.stdout:
        bwa_process.stdout.close()

    sort_return_code = sort_process.wait()
    bwa_return_code = bwa_process.wait()

    if bwa_return_code != 0:
        raise subprocess.CalledProcessError(bwa_return_code, bwa_command)
    if sort_return_code != 0:
        raise subprocess.CalledProcessError(sort_return_code, sort_command)


def gatk_command(config, tool_name, args):
    command = [config["tools"]["gatk"]]
    java_options = config.get("gatk_java_options")
    if java_options:
        command.extend(["--java-options", java_options])
    command.append(tool_name)
    command.extend(args)
    return command


def run_pipeline(sample_name, read1, read2, output_dir, config):
    reference = config["reference_fasta"]
    known_sites = config.get("known_sites", [])
    tools = config["tools"]

    output_dir.mkdir(parents=True, exist_ok=True)
    sorted_bam = output_dir / f"{sample_name}.sorted.bam"
    marked_bam = output_dir / f"{sample_name}.marked.bam"
    metrics = output_dir / f"{sample_name}.markdup.metrics.txt"
    recal_table = output_dir / f"{sample_name}.recal.table"
    recal_bam = output_dir / f"{sample_name}.recal.bam"
    raw_vcf = output_dir / f"{sample_name}.raw.vcf.gz"
    log_path = output_dir / "pipeline.log"
    status_path = output_dir / "status.txt"

    with log_path.open("w", encoding="utf-8") as log:
        try:
            status_path.write_text("RUNNING\n", encoding="utf-8")
            run_bwa_and_sort(sample_name, read1, read2, sorted_bam, config, log)
            run_command([tools["samtools"], "index", sorted_bam], log)

            run_command(
                gatk_command(
                    config,
                    "MarkDuplicates",
                    [
                        "-I",
                        sorted_bam,
                        "-O",
                        marked_bam,
                        "-M",
                        metrics,
                        "--CREATE_INDEX",
                        "true",
                    ],
                ),
                log,
            )

            if known_sites:
                bqsr_args = [
                    "-R",
                    reference,
                    "-I",
                    marked_bam,
                    "-O",
                    recal_table,
                ]
                for known_site in known_sites:
                    bqsr_args.extend(["--known-sites", known_site])

                run_command(gatk_command(config, "BaseRecalibrator", bqsr_args), log)
                run_command(
                    gatk_command(
                        config,
                        "ApplyBQSR",
                        [
                            "-R",
                            reference,
                            "-I",
                            marked_bam,
                            "--bqsr-recal-file",
                            recal_table,
                            "-O",
                            recal_bam,
                        ],
                    ),
                    log,
                )
                caller_input = recal_bam
            else:
                log.write("\nNo known_sites configured. Skipping BQSR.\n")
                caller_input = marked_bam

            run_command(
                gatk_command(
                    config,
                    "HaplotypeCaller",
                    [
                        "-R",
                        reference,
                        "-I",
                        caller_input,
                        "-O",
                        raw_vcf,
                    ],
                ),
                log,
            )
            status_path.write_text("DONE\n", encoding="utf-8")
        except Exception as error:
            log.write(f"\nFAILED: {error}\n")
            status_path.write_text("FAILED\n", encoding="utf-8")
            raise


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    parser.add_argument("--read1", required=True)
    parser.add_argument("--read2", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--config", default="config.json")
    args = parser.parse_args()

    with Path(args.config).open("r", encoding="utf-8") as handle:
        config = json.load(handle)

    run_pipeline(
        args.sample,
        Path(args.read1),
        Path(args.read2),
        Path(args.output_dir),
        config,
    )


if __name__ == "__main__":
    main()
