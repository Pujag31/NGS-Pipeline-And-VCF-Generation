import streamlit as st
import subprocess, os, shutil, zipfile, platform, multiprocessing, threading, io
import pandas as pd
from datetime import datetime
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib import colors

# ─────────────────────────────────────────────────────────────────────
# STORAGE PHILOSOPHY — ZERO PERSISTENT FILES
#
#   ALL files (FASTQ upload, BAM, VCF, PDF) live only in temp/
#   As soon as a file is needed for download it is read into memory
#   (bytes) and handed to st.download_button, then temp/ is wiped.
#   Nothing is ever written to output/, uploads/, or any other folder.
#   Disk usage at any moment = one sample's working files only.
# ─────────────────────────────────────────────────────────────────────

st.set_page_config(page_title="NGS WES/WGS Pipeline", layout="wide")
st.title("🧬 NGS WES/WGS Pipeline")
st.success("Frontend Running  |  Zero-storage mode: all files auto-deleted after download")

BASE_DIR  = os.path.dirname(os.path.abspath(__file__))
REFERENCE = os.path.join(BASE_DIR, "reference", "genome", "GrCh38", "genome.fa")

# ── CPU detection ─────────────────────────────────────────────────────
TOTAL_CPUS = multiprocessing.cpu_count()
THREADS = min(4, TOTAL_CPUS)
SORT_MEM   = "4G"

# ── Only temp/ and logs/ are created — nothing else ──────────────────
os.makedirs("temp", exist_ok=True)
os.makedirs("logs", exist_ok=True)
LOG_FILE   = "logs/pipeline.log"
IS_WINDOWS = platform.system() == "Windows"

# ── Tool check ────────────────────────────────────────────────────────
def tool_available(name):
    try:
        cmd = "where" if IS_WINDOWS else "which"
        return subprocess.run([cmd, name], capture_output=True).returncode == 0
    except:
        return False

TOOLS = {t: tool_available(t) for t in ["fastqc", "bwa", "samtools", "bcftools"]}

with st.expander("🔧 Tool Status", expanded=False):
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("FastQC",   "✅" if TOOLS["fastqc"]   else "❌")
    c2.metric("BWA",      "✅" if TOOLS["bwa"]      else "❌")
    c3.metric("Samtools", "✅" if TOOLS["samtools"] else "❌")
    c4.metric("Bcftools", "✅" if TOOLS["bcftools"] else "❌")
    missing = [t for t, ok in TOOLS.items() if not ok]
    if missing:
        st.warning(f"Missing: {', '.join(missing)}")
        st.code("conda install -c bioconda fastqc bwa samtools bcftools")

if not os.path.exists(REFERENCE):
    st.error(f"Reference genome not found: `{REFERENCE}`")
else:
    st.info(f"Reference: `{REFERENCE}`")

# ── Sidebar ───────────────────────────────────────────────────────────
st.sidebar.header("⚙️ Pipeline Options")
st.sidebar.caption(f"CPUs: {TOTAL_CPUS}  |  Threads: {THREADS}")
st.sidebar.info(
    "**Zero-storage mode**\n\n"
    "No files are saved to disk permanently.\n"
    "Everything runs in temp/ and is wiped after each sample.\n"
    "Download your VCF + PDF before moving to the next sample."
)

# ── Logger ────────────────────────────────────────────────────────────
def log(msg):
    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    with open(LOG_FILE, "a") as f:
        f.write(f"[{ts}] {msg}\n")

def run_cmd(cmd, step):
    log(f"START: {step}")
    r = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if r.returncode != 0:
        log(f"FAILED: {step} -- {r.stderr[-300:]}")
        raise Exception(r.stderr[-500:])
    log(f"DONE: {step}")
    return r



# ── FastQC parser ─────────────────────────────────────────────────────
def parse_fastqc(txt_path, sample):
    data = {"Sample": sample, "Reads": 0, "GC": 0, "Length": "0"}
    with open(txt_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if "Total Sequences" in line:
                data["Reads"]  = int(line.split("\t")[-1].strip())
            elif "%GC" in line:
                data["GC"]     = int(line.split("\t")[-1].strip())
            elif "Sequence length" in line:
                data["Length"] = line.split("\t")[-1].strip()
    return data

# ── PDF built into memory (no disk write) ────────────────────────────
def generate_pdf_bytes(sample, qc_data, total_variants, snps, indels,
                       flagstat_text, elapsed_align):
    buf    = io.BytesIO()
    doc    = SimpleDocTemplate(buf)
    styles = getSampleStyleSheet()
    story  = [
        Paragraph(f"NGS Pipeline Report -- {sample}", styles["Title"]),
        Spacer(1, 10),
        Paragraph(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}", styles["Normal"]),
        Paragraph(f"Threads used: {THREADS} of {TOTAL_CPUS} available", styles["Normal"]),
        Spacer(1, 10),
        Paragraph("QC Summary (FastQC)", styles["Heading2"]),
        Paragraph(
            f"<b>Total Reads:</b> {qc_data['Reads']:,}<br/>"
            f"<b>GC Content:</b> {qc_data['GC']}%<br/>"
            f"<b>Read Length:</b> {qc_data['Length']}",
            styles["BodyText"]
        ),
        Spacer(1, 10),
        Paragraph("Alignment Summary (BWA MEM)", styles["Heading2"]),
        Paragraph(
            f"<b>Alignment time:</b> {elapsed_align}s<br/>"
            f"<b>Flagstat:</b><br/>{flagstat_text.replace(chr(10), '<br/>')}",
            styles["BodyText"]
        ),
        Spacer(1, 10),
        Paragraph("Variant Calling Summary (bcftools)", styles["Heading2"]),
        Paragraph(
            f"<b>Total Variants:</b> {total_variants}<br/>"
            f"<b>SNPs:</b> {snps}<br/>"
            f"<b>Indels:</b> {indels}",
            styles["BodyText"]
        ),
        Spacer(1, 10),
        Paragraph("Storage", styles["Heading2"]),
        Paragraph(
            "Zero-storage mode: all intermediate files were deleted from disk "
            "after pipeline completion. Only downloaded files were retained.",
            styles["BodyText"]
        ),
    ]
    doc.build(story)
    buf.seek(0)
    return buf.read()

# ── ALIGNMENT — BWA | samtools view | samtools sort, fully piped ──────
def run_alignment_fast(r1_path, r2_path, sample):
    """
    Pipe chain: bwa mem -> samtools view -> samtools sort -> temp BAM.
    No SAM is ever written to disk.
    stderr of all three processes is drained by background threads
    to prevent OS pipe buffer deadlocks.
    """
    bam_out = os.path.join("temp", f"{sample}.bam")
    log(f"START: bwa mem piped (threads={THREADS})")

    bwa_err_buf  = []
    view_err_buf = []
    sort_err_buf = []

    def _drain(stream, buf):
        try:
            buf.append(stream.read().decode(errors="replace"))
        except Exception:
            pass

    # 1. BWA MEM
    bwa = subprocess.Popen(
    [
        "bwa", "mem",
        "-t", str(THREADS),
        "-R", f"@RG\\tID:{sample}\\tSM:{sample}\\tPL:ILLUMINA",
        REFERENCE,
        r1_path,
        r2_path
    ],
    stdout=subprocess.PIPE,
    stderr=subprocess.PIPE
)

    # 2. samtools view (SAM -> unsorted BAM stream)
    view = subprocess.Popen(
        ["samtools", "view", "-@", str(THREADS), "-Sb", "-"],
        stdin=bwa.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    bwa.stdout.close()
    threading.Thread(target=_drain, args=(view.stderr, view_err_buf), daemon=True).start()

    # 3. samtools sort (unsorted BAM stream -> sorted BAM on disk)
    sort = subprocess.Popen(
        ["samtools", "sort", "-@", str(THREADS), "-m", SORT_MEM, "-o", bam_out, "-"],
        stdin=view.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    view.stdout.close()
    t_sort_err = threading.Thread(target=_drain, args=(sort.stderr, sort_err_buf), daemon=True)
    t_sort_err.start()

    # Wait tail-to-head
    sort.wait(); view.wait(); bwa.wait()
    t_sort_err.join(timeout=5)

    if bwa.returncode  != 0:
        raise Exception(f"BWA MEM failed (rc={bwa.returncode}):\n{''.join(bwa_err_buf)[-600:]}")
    if view.returncode != 0:
        raise Exception(f"samtools view failed (rc={view.returncode}):\n{''.join(view_err_buf)[-400:]}")
    if sort.returncode != 0:
        raise Exception(f"samtools sort failed (rc={sort.returncode}):\n{''.join(sort_err_buf)[-400:]}")

    log(f"DONE: bwa+sort -> {bam_out}")

    # Index
    run_cmd(["samtools", "index", bam_out], "samtools index")
    return bam_out

# ─────────────────────────────────────────────────────────────────────
# UPLOAD
# ─────────────────────────────────────────────────────────────────────
r1_file = st.file_uploader(
    "📂 Upload R1 FASTQ",
    type=["fastq", "fq", "gz"],
    key="r1"
)

r2_file = st.file_uploader(
    "📂 Upload R2 FASTQ",
    type=["fastq", "fq", "gz"],
    key="r2"
)

# ─────────────────────────────────────────────────────────────────────
# MAIN PIPELINE
# ─────────────────────────────────────────────────────────────────────
if files:
    qc_dashboard = []

    for file in files:
        sample = file.name.split(".")[0]

        st.markdown("---")
        st.header(f"🧬 {sample}")



        # ── Save uploaded FASTQ to temp (only copy on disk) ───────────
        r1_path = os.path.join("temp", r1_file.name)
        r2_path = os.path.join("temp", r2_file.name)
        with open(r1_path, "wb") as f:
            f.write(r1_file.getbuffer())
            with open(r2_path, "wb") as f:
                f.write(r2_file.getbuffer())

        # ─────────────────────────────────────────────────────────────
        # STEP 1 — FastQC
        # ─────────────────────────────────────────────────────────────
        st.subheader("🧪 Step 1 — FastQC")
        qc_data = {"Sample": sample, "Reads": 0, "GC": 0, "Length": "0"}

        if TOOLS["fastqc"]:
            try:
                run_cmd(["fastqc", input_path, "-o", "temp", "-t", str(THREADS)], "FastQC")

                # Find and extract fastqc_data.txt from the zip
                qc_zip = next(
                    (os.path.join("temp", f) for f in os.listdir("temp") if f.endswith("_fastqc.zip")),
                    None
                )
                if qc_zip:
                    with zipfile.ZipFile(qc_zip) as z:
                        txt_entry = next(x for x in z.namelist() if "fastqc_data.txt" in x)
                        z.extract(txt_entry, "temp")
                        qc_data = parse_fastqc(os.path.join("temp", txt_entry), sample)

                    # Delete FastQC zip + html — we only needed the parsed numbers
                    for f in os.listdir("temp"):
                        if f.endswith("_fastqc.zip") or f.endswith("_fastqc.html"):
                            os.remove(os.path.join("temp", f))

                qc_dashboard.append(qc_data)
                r1, r2, r3 = st.columns(3)
                r1.metric("Total Reads",  f"{qc_data['Reads']:,}")
                r2.metric("GC Content",   f"{qc_data['GC']}%")
                r3.metric("Read Length",  qc_data["Length"])
                st.success("FastQC done — zip/html deleted from disk")

            except Exception as e:
                st.error(f"FastQC failed: {e}")
        else:
            st.warning("FastQC not installed — skipping")

        # ─────────────────────────────────────────────────────────────
        # STEP 2 — Alignment
        # ─────────────────────────────────────────────────────────────
        st.subheader(f"🔗 Step 2 — Alignment (BWA MEM, {THREADS} threads)")

        if not os.path.exists(REFERENCE):
            st.error(f"Reference not found: `{REFERENCE}` — cannot align")
            wipe_temp()
            continue

        if not (TOOLS["bwa"] and TOOLS["samtools"]):
            st.error(f"Missing: {[t for t in ['bwa','samtools'] if not TOOLS[t]]}")
            wipe_temp()
            continue

        flagstat_text = ""
        elapsed_align = 0
        try:
            t0 = datetime.now()
            with st.spinner(f"BWA MEM → samtools view → samtools sort (piped, {THREADS} threads)…"):
                bam = run_alignment_fast(input_path, sample)
            elapsed_align = int((datetime.now() - t0).total_seconds())

            # Delete the uploaded FASTQ now — no longer needed
            try:
                os.remove(input_path)
                log(f"Deleted FASTQ from temp: {input_path}")
            except Exception:
                pass

            # Flagstat
            flagstat_text = subprocess.run(
                ["samtools", "flagstat", "-@", str(THREADS), bam],
                capture_output=True, text=True
            ).stdout

            mapped_pct = "N/A"
            for line in flagstat_text.split("\n"):
                if "mapped (" in line:
                    mapped_pct = line.split("(")[1].split("%")[0].strip() + "%"
                    break

            a1, a2 = st.columns(2)
            a1.metric("Alignment time", f"{elapsed_align}s")
            a2.metric("Reads mapped",   mapped_pct)
            with st.expander("Flagstat details"):
                st.code(flagstat_text, language="bash")
            st.success(f"BAM sorted + indexed in temp ({os.path.getsize(bam)/1e6:.1f} MB)")

        except Exception as e:
            st.error(f"Alignment failed: {e}")
            wipe_temp()
            continue

        # ─────────────────────────────────────────────────────────────
        # STEP 3 — Variant Calling
        # ─────────────────────────────────────────────────────────────
        st.subheader("🧬 Step 3 — Variant Calling (bcftools)")

        if not TOOLS["bcftools"]:
            st.error("bcftools not installed")
            wipe_temp()
            continue

        vcf_path      = os.path.join("temp", f"{sample}.vcf")
        snp_vcf_path  = os.path.join("temp", f"{sample}_snps.vcf")
        indel_vcf_path= os.path.join("temp", f"{sample}_indels.vcf")
        total_variants = snps = indels = 0

        try:
            with st.spinner("bcftools mpileup | bcftools call…"):
                mpileup = subprocess.Popen(
                    ["bcftools", "mpileup", "--threads", str(THREADS), "-f", REFERENCE, bam],
                    stdout=subprocess.PIPE, stderr=subprocess.PIPE
                )
                with open(vcf_path, "w") as vcf_out:
                    call = subprocess.run(
                        ["bcftools", "call", "-mv", "-Ov", "--threads", str(THREADS)],
                        stdin=mpileup.stdout, stdout=vcf_out, stderr=subprocess.PIPE
                    )
                mpileup.wait()
                if call.returncode != 0:
                    raise Exception(call.stderr.decode()[-400:])

            # Delete BAM + BAI now — variant calling is done
            for ext in [bam, bam + ".bai"]:
                try:
                    os.remove(ext)
                except Exception:
                    pass
            log("BAM + BAI deleted from temp after variant calling")

            # Split SNPs / Indels
            run_cmd(["bcftools", "view", "--threads", str(THREADS),
                     "-v", "snps",   vcf_path, "-o", snp_vcf_path],  "SNP split")
            run_cmd(["bcftools", "view", "--threads", str(THREADS),
                     "-v", "indels", vcf_path, "-o", indel_vcf_path], "Indel split")

            # Count variants
            with open(vcf_path) as vf:
                for line in vf:
                    if not line.startswith("#"):
                        total_variants += 1
                        p = line.split("\t")
                        if len(p) >= 5:
                            if len(p[3]) == len(p[4]): snps   += 1
                            else:                       indels += 1

            v1, v2, v3 = st.columns(3)
            v1.metric("Total Variants", total_variants)
            v2.metric("SNPs",           snps)
            v3.metric("Indels",         indels)
            st.success("VCF generated")

        except Exception as e:
            st.error(f"Variant calling failed: {e}")
            wipe_temp()
            continue

        # ─────────────────────────────────────────────────────────────
        # STEP 4 — PDF report (built entirely in memory)
        # ─────────────────────────────────────────────────────────────
        st.subheader("📄 Step 4 — Report")
        pdf_bytes = None
        try:
            pdf_bytes = generate_pdf_bytes(
                sample, qc_data, total_variants, snps, indels,
                flagstat_text, elapsed_align
            )
            st.success("PDF built in memory — ready to download")
        except Exception as e:
            st.error(f"PDF error: {e}")

        # ─────────────────────────────────────────────────────────────
        # DOWNLOADS — read files into memory, then delete from disk
        # ─────────────────────────────────────────────────────────────
        st.subheader("⬇️ Downloads")
        st.info(
            "Download your files now. "
            "Disk copies will be deleted immediately after this section."
        )

        vcf_bytes   = slurp_and_delete(vcf_path)      if os.path.exists(vcf_path)       else None
        snp_bytes   = slurp_and_delete(snp_vcf_path)  if os.path.exists(snp_vcf_path)   else None
        indel_bytes = slurp_and_delete(indel_vcf_path) if os.path.exists(indel_vcf_path) else None

        # Main VCF — primary button
        if vcf_bytes:
            st.download_button(
                f"⬇️ Full VCF — {sample}.vcf",
                data=vcf_bytes,
                file_name=f"{sample}.vcf",
                mime="text/plain",
                key=f"vcf_{sample}",
                use_container_width=True,
                type="primary",
            )

        col1, col2, col3 = st.columns(3)
        with col1:
            if snp_bytes:
                st.download_button(
                    "⬇️ SNPs VCF",
                    data=snp_bytes,
                    file_name=f"{sample}_snps.vcf",
                    mime="text/plain",
                    key=f"snp_{sample}",
                )
        with col2:
            if indel_bytes:
                st.download_button(
                    "⬇️ Indels VCF",
                    data=indel_bytes,
                    file_name=f"{sample}_indels.vcf",
                    mime="text/plain",
                    key=f"indel_{sample}",
                )
        with col3:
            if pdf_bytes:
                st.download_button(
                    "⬇️ PDF Report",
                    data=pdf_bytes,
                    file_name=f"{sample}_report.pdf",
                    mime="application/pdf",
                    key=f"pdf_{sample}",
                )

        # ── Final wipe — nothing left on disk ─────────────────────────
        wipe_temp()
        st.success(
            f"✅ {sample} complete — temp/ wiped. "
            "Zero bytes remain on disk for this sample."
        )

        # Disk usage sanity check
        try:
            remaining = sum(
                os.path.getsize(os.path.join("temp", f))
                for f in os.listdir("temp")
            )
            st.caption(f"temp/ size after wipe: {remaining} bytes")
        except Exception:
            pass

    # ── QC Dashboard (all samples) ────────────────────────────────────
    if qc_dashboard:
        st.markdown("---")
        st.header("📊 QC Dashboard — All Samples")
        df = pd.DataFrame(qc_dashboard)
        st.dataframe(df, use_container_width=True)
        col1, col2 = st.columns(2)
        with col1:
            st.subheader("Reads per Sample")
            st.bar_chart(df.set_index("Sample")["Reads"])
        with col2:
            st.subheader("GC% per Sample")
            st.bar_chart(df.set_index("Sample")["GC"])

    # ── Pipeline Log ──────────────────────────────────────────────────
    with st.expander("📋 Pipeline Log"):
        if os.path.exists(LOG_FILE):
            with open(LOG_FILE) as f:
                st.code(f.read(), language="bash")

else:
    st.info("Upload one or more FASTQ / FASTQ.GZ files to begin.")
    st.markdown(f"""
**Pipeline steps (all in temp/, nothing saved permanently):**

| Step | Tool | What happens to files |
|------|------|-----------------------|
| 1 | FastQC | Run → parse metrics → zip/html deleted immediately |
| 2 | BWA MEM + samtools | Pipe: no SAM ever written; FASTQ deleted after BAM is ready |
| 3 | bcftools | VCF written to temp; BAM deleted immediately after calling |
| 4 | PDF | Built entirely in memory — never touches disk |
| ⬇️ | Downloads | VCF/PDF read into RAM → served to browser → temp/ wiped |

**Performance:** {THREADS} of {TOTAL_CPUS} CPU cores used across all steps.
BWA + view + sort run as a single pipe — no intermediate SAM file.
    """)