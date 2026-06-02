import subprocess
import uuid
from pathlib import Path

import streamlit as st


BASE_DIR = Path(__file__).resolve().parent
INPUT_DIR = BASE_DIR / "input_fastq"
RUN_DIR = BASE_DIR / "runs"
CONFIG_PATH = BASE_DIR / "config.json"
RUNNER_PATH = BASE_DIR / "gatk_runner.py"


def list_fastqs():
    INPUT_DIR.mkdir(exist_ok=True)
    files = []
    for pattern in ("*.fastq", "*.fq", "*.fastq.gz", "*.fq.gz"):
        files.extend(INPUT_DIR.glob(pattern))
    return sorted(files, key=lambda path: path.name.lower())


def safe_sample_name(value):
    cleaned = "".join(char if char.isalnum() or char in "._-" else "_" for char in value)
    return cleaned.strip("._-") or "sample"


def get_status(run_dir):
    status_path = run_dir / "status.txt"
    if status_path.exists():
        return status_path.read_text(encoding="utf-8").strip()
    return "STARTED"


st.set_page_config(page_title="GATK Pipeline", layout="centered")
st.title("GATK Pipeline")

if not CONFIG_PATH.exists():
    st.error("config.json is missing.")
    st.stop()

fastqs = list_fastqs()
if not fastqs:
    st.warning("Put FASTQ files inside input_fastq, then refresh.")
    st.stop()

sample_name = safe_sample_name(st.text_input("Sample name", "sample_001"))
read1 = st.selectbox("FASTQ R1", fastqs, format_func=lambda path: path.name)
read2 = st.selectbox("FASTQ R2", fastqs, format_func=lambda path: path.name)

if st.button("Start Pipeline", type="primary"):
    if read1 == read2:
        st.error("R1 and R2 must be different files.")
        st.stop()

    run_dir = RUN_DIR / f"{sample_name}_{uuid.uuid4().hex[:12]}"
    run_dir.mkdir(parents=True, exist_ok=True)

    stdout_path = run_dir / "launcher.log"
    command = [
        "python3",
        str(RUNNER_PATH),
        "--sample",
        sample_name,
        "--read1",
        str(read1),
        "--read2",
        str(read2),
        "--output-dir",
        str(run_dir),
        "--config",
        str(CONFIG_PATH),
    ]

    with stdout_path.open("w", encoding="utf-8") as stdout:
        subprocess.Popen(
            command,
            cwd=BASE_DIR,
            stdout=stdout,
            stderr=subprocess.STDOUT,
            start_new_session=True,
        )

    st.success("Pipeline started. You can refresh this page to check progress.")
    st.code(str(run_dir), language="text")

st.subheader("Runs")
run_dirs = sorted([path for path in RUN_DIR.glob("*") if path.is_dir()], reverse=True)

for run_dir in run_dirs[:10]:
    status = get_status(run_dir)
    st.write(f"{run_dir.name}: {status}")
    log_path = run_dir / "pipeline.log"
    vcf_files = list(run_dir.glob("*.vcf.gz"))

    if vcf_files:
        with vcf_files[0].open("rb") as handle:
            st.download_button(
                f"Download {vcf_files[0].name}",
                data=handle,
                file_name=vcf_files[0].name,
                key=f"vcf-{run_dir.name}",
            )

    if log_path.exists():
        with st.expander(f"Log: {run_dir.name}"):
            st.text_area(
                "pipeline.log",
                log_path.read_text(encoding="utf-8", errors="replace")[-8000:],
                height=240,
                key=f"log-{run_dir.name}",
            )
