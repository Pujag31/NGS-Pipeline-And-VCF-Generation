import streamlit as st
import subprocess, os, shutil, zipfile, platform, multiprocessing, threading, io, time
import pandas as pd
from datetime import datetime
from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer, Table,
                                 TableStyle, HRFlowable)
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib import colors
from reportlab.lib.units import mm
from reportlab.lib.pagesizes import A4

# ═══════════════════════════════════════════════════════════════════════
#  CONFIGURATION — edit this block only
# ═══════════════════════════════════════════════════════════════════════
BASE_DIR    = os.path.dirname(os.path.abspath(__file__))

# Reference genome (GRCh38)
REFERENCE   = os.path.join(BASE_DIR, "reference", "genome", "GrCh38", "genome.fa")

# Known-sites VCFs required for BQSR (download from GATK resource bundle)
KNOWN_DBSNP  = os.path.join(BASE_DIR, "reference", "known_sites", "dbsnp_146.hg38.vcf.gz")
KNOWN_MILLS  = os.path.join(BASE_DIR, "reference", "known_sites", "Mills_and_1000G_gold_standard.indels.hg38.vcf.gz")
KNOWN_1000G  = os.path.join(BASE_DIR, "reference", "known_sites", "1000G_phase1.snps.high_confidence.hg38.vcf.gz")

# 1 TB server mount point
SERVER_ROOT  = os.path.join(BASE_DIR, "server_storage")   # e.g. "/mnt/ngs_1tb"

# ═══════════════════════════════════════════════════════════════════════
#  STORAGE LAYOUT
#
#  SERVER (1 TB) — auto-saved, organised by date:
#    SERVER_ROOT/
#      QC_Reports/YYYY-MM-DD/   <sample>_QC_Report.pdf   ← auto
#      VCF/YYYY-MM-DD/          <sample>.g.vcf.gz         ← auto (GVCF)
#                               <sample>.vcf.gz           ← auto (genotyped)
#                               <sample>.filtered.vcf.gz  ← auto (filtered)
#      logs/                    pipeline.log
#
#  TEMP — wiped after each sample:
#    FASTQ → trimmed FASTQ → BAM → dedup BAM → BQSR BAM → GVCF → VCF
#    (each file deleted the moment it is no longer needed)
#
#  DOWNLOADS (browser buttons):
#    Filtered VCF              — always shown
#    Raw GVCF                  — always shown
#    BAM (post-BQSR)           — always shown, download-only, not saved to server
#    QC PDF                    — always shown
# ═══════════════════════════════════════════════════════════════════════

TODAY      = datetime.now().strftime("%Y-%m-%d")
SERVER_QC  = os.path.join(SERVER_ROOT, "QC_Reports", TODAY)
SERVER_VCF = os.path.join(SERVER_ROOT, "VCF",        TODAY)
SERVER_LOG = os.path.join(SERVER_ROOT, "logs")

for d in [SERVER_QC, SERVER_VCF, SERVER_LOG, "temp"]:
    os.makedirs(d, exist_ok=True)

LOG_FILE   = os.path.join(SERVER_LOG, "pipeline.log")
IS_WINDOWS = platform.system() == "Windows"

TOTAL_CPUS = multiprocessing.cpu_count()
THREADS    = max(1, TOTAL_CPUS - 1)
SORT_MEM   = "4G"

# ═══════════════════════════════════════════════════════════════════════
#  PAGE
# ═══════════════════════════════════════════════════════════════════════
st.set_page_config(page_title="NGS WES/WGS Pipeline", layout="wide")
st.title("🧬 NGS WES/WGS Pipeline  —  GATK Best Practices")

# ── Tool check ────────────────────────────────────────────────────────
def tool_available(name):
    try:
        return subprocess.run(
            ["where" if IS_WINDOWS else "which", name],
            capture_output=True).returncode == 0
    except Exception:
        return False

TOOLS = {t: tool_available(t)
         for t in ["fastqc", "fastp", "bwa", "samtools", "gatk", "bcftools"]}

with st.expander("🔧 Tool Status", expanded=True):
    c1,c2,c3,c4,c5,c6 = st.columns(6)
    c1.metric("FastQC",   "✅" if TOOLS["fastqc"]   else "❌")
    c2.metric("fastp",    "✅" if TOOLS["fastp"]    else "❌")
    c3.metric("BWA",      "✅" if TOOLS["bwa"]      else "❌")
    c4.metric("Samtools", "✅" if TOOLS["samtools"] else "❌")
    c5.metric("GATK4",    "✅" if TOOLS["gatk"]     else "❌")
    c6.metric("bcftools", "✅" if TOOLS["bcftools"] else "❌")
    missing = [t for t,ok in TOOLS.items() if not ok]
    if missing:
        st.warning(f"Missing: {', '.join(missing)}")
        st.code(
            "conda install -c bioconda fastqc fastp bwa samtools bcftools\n"
            "conda install -c bioconda gatk4"
        )

# Reference & known-sites check
ref_ok = os.path.exists(REFERENCE)
dbsnp_ok  = os.path.exists(KNOWN_DBSNP)
mills_ok  = os.path.exists(KNOWN_MILLS)
kg_ok     = os.path.exists(KNOWN_1000G)
bqsr_ready = dbsnp_ok and mills_ok and kg_ok

if not ref_ok:
    st.error(f"❌ Reference not found: `{REFERENCE}`")
else:
    st.info(f"✅ Reference: `{REFERENCE}`")

with st.expander("📚 Known-sites (BQSR)", expanded=not bqsr_ready):
    st.markdown(f"{'✅' if dbsnp_ok  else '❌'} dbSNP:   `{KNOWN_DBSNP}`")
    st.markdown(f"{'✅' if mills_ok  else '❌'} Mills:   `{KNOWN_MILLS}`")
    st.markdown(f"{'✅' if kg_ok     else '❌'} 1000G:   `{KNOWN_1000G}`")
    if not bqsr_ready:
        st.warning("BQSR will be skipped until all three known-sites VCFs are present. "
                   "Download from: https://gatk.broadinstitute.org/hc/en-us/articles/360035890811")

# Server storage status
try:
    su = shutil.disk_usage(SERVER_ROOT)
    pct = su.used/su.total*100
    s1,s2,s3,s4 = st.columns(4)
    s1.metric("Server Used",  f"{su.used/1e12:.3f} TB")
    s2.metric("Server Free",  f"{su.free/1e12:.3f} TB")
    s3.metric("Server Total", f"{su.total/1e12:.2f} TB")
    s4.metric("Usage",        f"{pct:.1f}%")
    if pct > 90: st.error("🔴 Server >90% full — archive immediately.")
    elif pct > 75: st.warning("🟡 Server >75% full — plan archiving soon.")
except Exception:
    st.caption(f"📁 Server root: `{SERVER_ROOT}`")

# Sidebar
st.sidebar.header("⚙️ Pipeline Options")
st.sidebar.caption(f"CPUs: {TOTAL_CPUS}  |  Threads: {THREADS}")
st.sidebar.markdown("---")
st.sidebar.markdown("**Auto-saved to server:**")
st.sidebar.markdown(f"- 📄 QC PDF → `QC_Reports/{TODAY}/`")
st.sidebar.markdown(f"- 🧬 VCF (GVCF + filtered) → `VCF/{TODAY}/`")
st.sidebar.markdown("---")
st.sidebar.markdown("**Download buttons (on demand):**")
st.sidebar.markdown("- Filtered VCF + GVCF")
st.sidebar.markdown("- QC PDF")
st.sidebar.markdown("- BAM (post-BQSR) — not saved to server")

# ═══════════════════════════════════════════════════════════════════════
#  HELPERS
# ═══════════════════════════════════════════════════════════════════════
def log(msg):
    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    with open(LOG_FILE, "a") as fh:
        fh.write(f"[{ts}] {msg}\n")

def run_cmd(cmd, step):
    log(f"START  {step}")
    r = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if r.returncode != 0:
        log(f"FAILED {step}\n{r.stderr[-400:]}")
        raise Exception(f"{step} failed:\n{r.stderr[-500:]}")
    log(f"DONE   {step}")
    return r

def wipe_temp():
    shutil.rmtree("temp", ignore_errors=True)
    os.makedirs("temp", exist_ok=True)

def mb(path):
    try:    return os.path.getsize(path) / 1e6
    except: return 0.0

def safe_remove(*paths):
    for p in paths:
        try: os.remove(p)
        except Exception: pass

def server_copy(src, dst, retries=3):
    for attempt in range(1, retries+1):
        try:
            shutil.copy2(src, dst)
            log(f"Saved to server: {dst}")
            return
        except Exception as e:
            log(f"server_copy attempt {attempt} failed: {e}")
            if attempt < retries: time.sleep(2)
    raise IOError(f"Failed to copy {src} → {dst} after {retries} attempts")

def t(path):
    """temp path helper"""
    return os.path.join("temp", path)

# ═══════════════════════════════════════════════════════════════════════
#  FASTQC PARSER
# ═══════════════════════════════════════════════════════════════════════
def parse_fastqc(txt_path, sample):
    d = {"Sample": sample, "Reads": 0, "GC": 0, "Length": "N/A",
         "Duplication": "N/A", "Adapter": "N/A"}
    section = ""
    with open(txt_path, "r", encoding="utf-8", errors="ignore") as fh:
        for raw in fh:
            line = raw.strip()
            if line.startswith(">>"): section = line
            if   "Total Sequences" in line: d["Reads"]      = int(line.split("\t")[-1])
            elif "%GC"             in line: d["GC"]         = int(line.split("\t")[-1])
            elif "Sequence length" in line: d["Length"]     = line.split("\t")[-1]
            elif "%Deduplicated"   in line: d["Duplication"] = line.split("\t")[-1].strip()+"% unique"
            elif "Adapter Content" in section:
                if   "PASS" in section: d["Adapter"] = "PASS"
                elif "FAIL" in section: d["Adapter"] = "FAIL"
                elif "WARN" in section: d["Adapter"] = "WARN"
    return d

# ═══════════════════════════════════════════════════════════════════════
#  ALIGNMENT — BWA MEM fully piped, no SAM on disk
# ═══════════════════════════════════════════════════════════════════════
def run_alignment(input_path, sample, bam_dest):
    log(f"START  BWA MEM piped (threads={THREADS})")
    bwa_e=[]; view_e=[]; sort_e=[]

    def _drain(s, b):
        try: b.append(s.read().decode(errors="replace"))
        except Exception: pass

    bwa = subprocess.Popen(
        ["bwa", "mem", "-t", str(THREADS),
         "-R", f"@RG\\tID:{sample}\\tSM:{sample}\\tPL:ILLUMINA\\tLB:{sample}\\tPU:{sample}",
         REFERENCE, input_path],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    threading.Thread(target=_drain, args=(bwa.stderr, bwa_e), daemon=True).start()

    view = subprocess.Popen(
        ["samtools", "view", "-@", str(THREADS), "-Sb", "-"],
        stdin=bwa.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    bwa.stdout.close()
    threading.Thread(target=_drain, args=(view.stderr, view_e), daemon=True).start()

    sort = subprocess.Popen(
        ["samtools", "sort", "-@", str(THREADS), "-m", SORT_MEM, "-o", bam_dest, "-"],
        stdin=view.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    view.stdout.close()
    ts = threading.Thread(target=_drain, args=(sort.stderr, sort_e), daemon=True)
    ts.start()

    sort.wait(); view.wait(); bwa.wait()
    ts.join(timeout=5)

    if bwa.returncode  != 0: raise Exception(f"BWA MEM failed:\n{''.join(bwa_e)[-600:]}")
    if view.returncode != 0: raise Exception(f"samtools view failed:\n{''.join(view_e)[-400:]}")
    if sort.returncode != 0: raise Exception(f"samtools sort failed:\n{''.join(sort_e)[-400:]}")

    log(f"DONE   BWA+sort → {bam_dest}")
    run_cmd(["samtools", "index", bam_dest], "samtools index")
    return bam_dest

# ═══════════════════════════════════════════════════════════════════════
#  MEDGENOME-STYLE QC PDF
# ═══════════════════════════════════════════════════════════════════════
def build_qc_pdf(sample, qc_raw, qc_trim, flagstat_text, elapsed_align,
                 dup_pct, bqsr_done, total_variants, snps, indels,
                 ts_filter, server_pdf_path):

    buf = io.BytesIO()
    doc = SimpleDocTemplate(buf, pagesize=A4,
                            leftMargin=15*mm, rightMargin=15*mm,
                            topMargin=15*mm, bottomMargin=15*mm)
    W = A4[0] - 30*mm

    NAVY  = colors.HexColor("#003366")
    GREEN = colors.HexColor("#006600")
    AMBER = colors.HexColor("#cc6600")
    RED   = colors.HexColor("#cc0000")
    LBLUE = colors.HexColor("#e8eef5")
    STRIP = colors.HexColor("#f5f8fc")

    Sty = lambda n, **kw: ParagraphStyle(n, **kw)
    title_s = Sty("T1", fontSize=16, fontName="Helvetica-Bold", textColor=NAVY, spaceAfter=3)
    sub_s   = Sty("T2", fontSize=10, fontName="Helvetica-Bold", textColor=NAVY, spaceAfter=2)
    body_s  = Sty("T3", fontSize=9,  fontName="Helvetica",      spaceAfter=2)
    mono_s  = Sty("T4", fontSize=7.5,fontName="Courier",        spaceAfter=3)
    foot_s  = Sty("T5", fontSize=7.5,fontName="Helvetica",      textColor=colors.grey)

    def HR(): return HRFlowable(width=W, thickness=1.5, color=NAVY, spaceAfter=4)
    def sec(txt): return [HR(), Paragraph(txt, sub_s)]

    def make_tbl(rows, widths, centre_from=1):
        tbl = Table(rows, colWidths=widths)
        tbl.setStyle(TableStyle([
            ("BACKGROUND",    (0,0),(-1,0),  NAVY),
            ("TEXTCOLOR",     (0,0),(-1,0),  colors.white),
            ("FONTNAME",      (0,0),(-1,0),  "Helvetica-Bold"),
            ("FONTSIZE",      (0,0),(-1,-1), 9),
            ("ROWBACKGROUNDS",(0,1),(-1,-1), [colors.white, STRIP]),
            ("GRID",          (0,0),(-1,-1), 0.4, colors.lightgrey),
            ("ALIGN",         (centre_from,0),(-1,-1), "CENTER"),
            ("TOPPADDING",    (0,0),(-1,-1), 4),
            ("BOTTOMPADDING", (0,0),(-1,-1), 4),
        ]))
        return tbl

    def colour_col(tbl_obj, rows, col=2):
        for i, row in enumerate(rows[1:], start=1):
            v = row[col]
            c = GREEN if v=="PASS" else AMBER if v=="WARN" else RED if v in("FAIL","CHECK") else None
            if c:
                tbl_obj.setStyle(TableStyle([
                    ("TEXTCOLOR",(col,i),(col,i),c),
                    ("FONTNAME", (col,i),(col,i),"Helvetica-Bold"),
                ]))

    # parse flagstat
    total_aln=mapped_n=mapped_pct_str="N/A"
    for line in flagstat_text.split("\n"):
        if "in total" in line: total_aln = line.split()[0]
        if "mapped (" in line and "primary" not in line:
            mapped_n = line.split()[0]
            try: mapped_pct_str = line.split("(")[1].split("%")[0].strip()+"%"
            except Exception: pass
    try: mp_f=float(mapped_pct_str.replace("%",""))
    except: mp_f=0.0
    mp_status = "PASS" if mp_f>=90 else ("WARN" if mp_f>=70 else "FAIL")

    story = []
    story += [
        Paragraph("NGS Quality Control Report", title_s),
        Paragraph("MedGenome-style Clinical Genomics Report — GATK Best Practices", sub_s),
        HRFlowable(width=W, thickness=2, color=NAVY, spaceAfter=6),
    ]

    # metadata
    meta = [
        ["Sample ID",     sample,              "Report Date",  datetime.now().strftime("%d-%b-%Y %H:%M")],
        ["Reference",     "GRCh38 / hg38",     "Pipeline",     "GATK4 Best Practices"],
        ["Aligner",       f"BWA MEM ({THREADS}t)","Caller",    "GATK HaplotypeCaller"],
        ["Align Time",    f"{elapsed_align}s",  "Threads",     f"{THREADS}/{TOTAL_CPUS}"],
        ["BQSR",          "Yes" if bqsr_done else "Skipped (no known-sites)",
         "Variant Filter","GATK VariantFiltration"],
    ]
    mt = Table(meta, colWidths=[32*mm,58*mm,32*mm,58*mm])
    mt.setStyle(TableStyle([
        ("BACKGROUND",    (0,0),(0,-1),LBLUE),
        ("BACKGROUND",    (2,0),(2,-1),LBLUE),
        ("FONTNAME",      (0,0),(-1,-1),"Helvetica"),
        ("FONTSIZE",      (0,0),(-1,-1),8),
        ("GRID",          (0,0),(-1,-1),0.4,colors.lightgrey),
        ("TOPPADDING",    (0,0),(-1,-1),3),
        ("BOTTOMPADDING", (0,0),(-1,-1),3),
    ]))
    story += [mt, Spacer(1,8)]

    # Sec 1: FastQC raw
    story += sec("1.  Raw Read Quality  (FastQC — pre-trim)")
    qr = [
        ["Metric","Value","Status"],
        ["Total Reads",    f"{qc_raw['Reads']:,}",  "—"],
        ["GC Content",     f"{qc_raw['GC']}%",
             "PASS" if 30<=qc_raw["GC"]<=65 else "WARN"],
        ["Read Length",    qc_raw["Length"],         "—"],
        ["Duplication",    qc_raw["Duplication"],    "—"],
        ["Adapter Content",qc_raw["Adapter"],
             "PASS" if qc_raw["Adapter"]=="PASS" else
             "WARN" if qc_raw["Adapter"]=="WARN" else "CHECK"],
    ]
    qt = make_tbl(qr, [60*mm,60*mm,60*mm])
    colour_col(qt, qr)
    story += [qt, Spacer(1,6)]

    # Sec 1b: fastp trimming summary
    story += sec("1b. Trimming Summary  (fastp)")
    if qc_trim:
        tr = [
            ["Metric","Value"],
            ["Reads Before Trim", f"{qc_trim.get('reads_before',0):,}"],
            ["Reads After Trim",  f"{qc_trim.get('reads_after',0):,}"],
            ["Reads Removed",     f"{qc_trim.get('reads_removed',0):,}"],
            ["Q20 Rate",          f"{qc_trim.get('q20',0):.1f}%"],
            ["Q30 Rate",          f"{qc_trim.get('q30',0):.1f}%"],
            ["Adapter Trimmed",   f"{qc_trim.get('adapter_trimmed',0):,}"],
        ]
        tt = make_tbl(tr, [90*mm,90*mm], centre_from=1)
        story += [tt, Spacer(1,6)]
    else:
        story.append(Paragraph("fastp not available — trimming was skipped.", body_s))
        story.append(Spacer(1,6))

    # Sec 2: Alignment
    story += sec("2.  Alignment Summary  (BWA MEM → samtools sort)")
    ar = [
        ["Metric","Value","Status"],
        ["Total Reads",    total_aln,          "—"],
        ["Mapped Reads",   mapped_n,           "—"],
        ["Mapping Rate",   mapped_pct_str,     mp_status],
        ["Alignment Time", f"{elapsed_align}s","—"],
    ]
    at = make_tbl(ar, [60*mm,60*mm,60*mm])
    colour_col(at, ar)
    story += [at, Spacer(1,6)]
    story += sec("2a. Flagstat Detail")
    story.append(Paragraph(flagstat_text.replace("\n","<br/>"), mono_s))
    story.append(Spacer(1,6))

    # Sec 3: Mark Duplicates
    story += sec("3.  Mark Duplicates  (GATK MarkDuplicates)")
    try:
        dup_f = float(str(dup_pct).replace("%",""))
        dup_status = "PASS" if dup_f<20 else ("WARN" if dup_f<30 else "FAIL")
    except Exception:
        dup_status = "—"
    dr = [
        ["Metric","Value","Status"],
        ["Duplication Rate", f"{dup_pct}", dup_status],
    ]
    dt = make_tbl(dr, [60*mm,60*mm,60*mm])
    colour_col(dt, dr)
    story += [dt, Spacer(1,6)]

    # Sec 4: BQSR
    story += sec("4.  Base Quality Score Recalibration  (GATK BQSR)")
    story.append(Paragraph(
        "BQSR applied using dbSNP + Mills + 1000G known sites." if bqsr_done
        else "⚠ BQSR skipped — known-sites VCFs not found. "
             "Download from the GATK resource bundle for accurate results.",
        body_s))
    story.append(Spacer(1,6))

    # Sec 5: Variant Calling
    story += sec("5.  Variant Calling  (GATK HaplotypeCaller)")
    ratio = f"{snps/indels:.2f}" if indels else "N/A"
    vr = [
        ["Metric","Count"],
        ["Total Variants",   f"{total_variants:,}"],
        ["SNPs",             f"{snps:,}"],
        ["Indels",           f"{indels:,}"],
        ["SNP/Indel Ratio",  ratio],
    ]
    vt2 = make_tbl(vr, [90*mm,90*mm], centre_from=1)
    story += [vt2, Spacer(1,6)]

    # Sec 6: Variant Filtering
    story += sec("6.  Variant Filtering  (GATK VariantFiltration — hard filters)")
    ff = [
        ["Filter","Expression","Applied"],
        ["SNP QD",     "QD < 2.0",          "✅"],
        ["SNP FS",     "FS > 60.0",         "✅"],
        ["SNP MQ",     "MQ < 40.0",         "✅"],
        ["SNP MQRankSum","MQRankSum < -12.5","✅"],
        ["SNP ReadPosRankSum","ReadPosRankSum < -8.0","✅"],
        ["Indel QD",   "QD < 2.0",          "✅"],
        ["Indel FS",   "FS > 200.0",        "✅"],
        ["Indel ReadPosRankSum","ReadPosRankSum < -20.0","✅"],
    ]
    ft = make_tbl(ff, [45*mm,75*mm,60*mm])
    story += [ft, Spacer(1,6)]

    # Sec 7: Filtering stats
    story += sec("7.  Post-filter Summary")
    pf = [
        ["Metric","Count"],
        ["Variants passing filter", f"{ts_filter.get('pass',0):,}"],
        ["Variants filtered out",   f"{ts_filter.get('filtered',0):,}"],
        ["Pass rate",               f"{ts_filter.get('pass_pct','N/A')}"],
    ]
    pft = make_tbl(pf, [90*mm,90*mm], centre_from=1)
    story += [pft, Spacer(1,6)]

    # Sec 8: Storage
    story += sec("8.  Storage Summary")
    for line in [
        f"QC PDF        : {server_pdf_path}",
        f"GVCF          : {os.path.join(SERVER_VCF, sample+'.g.vcf.gz')}",
        f"Genotyped VCF : {os.path.join(SERVER_VCF, sample+'.vcf.gz')}",
        f"Filtered VCF  : {os.path.join(SERVER_VCF, sample+'.filtered.vcf.gz')}",
        "BAM (BQSR)    : download-only, not saved to server",
        f"Pipeline log  : {LOG_FILE}",
        "Temp folder   : wiped after each sample",
    ]:
        story.append(Paragraph(f"• {line}", body_s))
    story.append(Spacer(1,6))

    story += [
        HRFlowable(width=W, thickness=0.5, color=colors.lightgrey, spaceAfter=3),
        Paragraph(
            f"NGS WES/WGS Pipeline  |  GATK Best Practices  |  GRCh38  |  "
            f"Generated {datetime.now().strftime('%d-%b-%Y %H:%M:%S')}",
            foot_s),
    ]

    doc.build(story)
    pdf_bytes = buf.getvalue()
    with open(server_pdf_path, "wb") as fh:
        fh.write(pdf_bytes)
    log(f"QC PDF saved: {server_pdf_path}")
    return pdf_bytes

# ═══════════════════════════════════════════════════════════════════════
#  SERVER BROWSER
# ═══════════════════════════════════════════════════════════════════════
def show_server_browser():
    with st.expander("📁 Server browser — today's files", expanded=False):
        for label, folder in [("QC Reports", SERVER_QC), ("VCF files", SERVER_VCF)]:
            try:   files_found = sorted(os.listdir(folder))
            except: files_found = []
            if files_found:
                st.markdown(f"**{label}** (`{folder}`)")
                rows = []
                for fn in files_found:
                    fp = os.path.join(folder, fn)
                    rows.append({"File": fn,
                                 "Size": f"{mb(fp):.1f} MB" if os.path.isfile(fp) else "—",
                                 "Modified": datetime.fromtimestamp(os.path.getmtime(fp)).strftime("%Y-%m-%d %H:%M") if os.path.isfile(fp) else "—"})
                st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)
            else:
                st.caption(f"No files yet in {label} for today.")

# ═══════════════════════════════════════════════════════════════════════
#  UPLOAD
# ═══════════════════════════════════════════════════════════════════════
files = st.file_uploader(
    "📂 Upload FASTQ / FASTQ.GZ  (single-end or paired R1+R2)",
    type=["fastq", "fq", "gz"],
    accept_multiple_files=True
)

# ═══════════════════════════════════════════════════════════════════════
#  MAIN PIPELINE
# ═══════════════════════════════════════════════════════════════════════
if files:
    # Pair R1 / R2 if both uploaded
    file_map = {}
    for f in files:
        base = f.name.replace("_R1","").replace("_R2","").replace("_1","").replace("_2","").split(".")[0]
        file_map.setdefault(base, []).append(f)

    qc_dashboard = []

    for sample, sample_files in file_map.items():
        st.markdown("---")
        st.header(f"🧬  {sample}")

        wipe_temp()

        # Save uploaded files
        r1_path = r2_path = None
        for sf in sample_files:
            dest = t(sf.name)
            with open(dest, "wb") as fh: fh.write(sf.getbuffer())
            if "_R2" in sf.name or "_2." in sf.name: r2_path = dest
            else:                                     r1_path = dest
        if r1_path is None: r1_path = t(sample_files[0].name)
        is_paired = r2_path is not None
        mode_str = "Paired-end" if is_paired else "Single-end"
        st.info(f"Mode: {mode_str}  |  R1: `{os.path.basename(r1_path)}`"
                + (f"  R2: `{os.path.basename(r2_path)}`" if is_paired else ""))

        # ─────────────────────────────────────────────────────────────
        # STEP 1 — FastQC (raw reads)
        # ─────────────────────────────────────────────────────────────
        st.subheader("🧪 Step 1 — FastQC  (raw reads)")
        qc_raw = {"Sample": sample, "Reads": 0, "GC": 0, "Length": "N/A",
                  "Duplication": "N/A", "Adapter": "N/A"}

        if TOOLS["fastqc"]:
            try:
                t0 = time.time()
                fqc_inputs = [r1_path] + ([r2_path] if is_paired else [])
                run_cmd(["fastqc"] + fqc_inputs + ["-o", "temp",
                         "-t", str(THREADS), "--nogroup"], "FastQC raw")
                qc_zip = next((t(f) for f in os.listdir("temp")
                               if f.endswith("_fastqc.zip")), None)
                if qc_zip:
                    with zipfile.ZipFile(qc_zip) as z:
                        txt_e = next(x for x in z.namelist() if "fastqc_data.txt" in x)
                        z.extract(txt_e, "temp")
                        qc_raw = parse_fastqc(os.path.join("temp", txt_e), sample)
                    for fn in os.listdir("temp"):
                        if fn.endswith("_fastqc.zip") or fn.endswith("_fastqc.html"):
                            safe_remove(t(fn))
                qc_dashboard.append(qc_raw)
                r1c, r2c, r3c, r4c = st.columns(4)
                r1c.metric("Total Reads",  f"{qc_raw['Reads']:,}")
                r2c.metric("GC Content",   f"{qc_raw['GC']}%")
                r3c.metric("Read Length",  qc_raw["Length"])
                r4c.metric("Duplication",  qc_raw["Duplication"])
                st.success(f"FastQC done ({time.time()-t0:.0f}s)")
            except Exception as e:
                st.error(f"FastQC failed: {e}")
        else:
            st.warning("FastQC not installed — skipping")

        # ─────────────────────────────────────────────────────────────
        # STEP 2 — Trimming (fastp)
        # ─────────────────────────────────────────────────────────────
        st.subheader("✂️ Step 2 — Adapter & Quality Trimming  (fastp)")
        qc_trim    = None
        trim_r1    = t(f"{sample}_trimmed_R1.fastq.gz")
        trim_r2    = t(f"{sample}_trimmed_R2.fastq.gz") if is_paired else None
        trim_json  = t(f"{sample}_fastp.json")
        trim_html  = t(f"{sample}_fastp.html")

        if TOOLS["fastp"]:
            try:
                t0 = time.time()
                fastp_cmd = [
                    "fastp",
                    "-i", r1_path,
                    "-o", trim_r1,
                    "--thread", str(min(THREADS, 16)),
                    "--json", trim_json,
                    "--html", trim_html,
                    "--detect_adapter_for_pe" if is_paired else "--detect_adapter",
                    "--qualified_quality_phred", "20",
                    "--length_required", "50",
                    "--correction",
                ]
                if is_paired:
                    fastp_cmd += ["-I", r2_path, "-O", trim_r2]
                run_cmd(fastp_cmd, "fastp trimming")

                # Parse JSON summary
                import json
                with open(trim_json) as jf:
                    jd = json.load(jf)
                smry = jd.get("summary", {})
                b4   = smry.get("before_filtering", {})
                af   = smry.get("after_filtering",  {})
                qc_trim = {
                    "reads_before":    b4.get("total_reads", 0),
                    "reads_after":     af.get("total_reads", 0),
                    "reads_removed":   b4.get("total_reads",0) - af.get("total_reads",0),
                    "q20":             af.get("q20_rate",0)*100,
                    "q30":             af.get("q30_rate",0)*100,
                    "adapter_trimmed": jd.get("filtering_result",{}).get("adapter_trimmed",0),
                }
                safe_remove(r1_path)
                if is_paired: safe_remove(r2_path)
                safe_remove(trim_json, trim_html)

                # Use trimmed reads downstream
                r1_path = trim_r1
                if is_paired: r2_path = trim_r2

                tc1,tc2,tc3,tc4 = st.columns(4)
                tc1.metric("Reads Before", f"{qc_trim['reads_before']:,}")
                tc2.metric("Reads After",  f"{qc_trim['reads_after']:,}")
                tc3.metric("Q20 Rate",     f"{qc_trim['q20']:.1f}%")
                tc4.metric("Q30 Rate",     f"{qc_trim['q30']:.1f}%")
                st.success(f"Trimming done ({time.time()-t0:.0f}s)")
            except Exception as e:
                st.error(f"fastp failed: {e}")
                st.warning("Continuing with untrimmed reads")
        else:
            st.warning("fastp not installed — skipping trimming (not recommended for WES/WGS)")

        # ─────────────────────────────────────────────────────────────
        # STEP 3 — Alignment (BWA MEM, piped)
        # ─────────────────────────────────────────────────────────────
        st.subheader(f"🔗 Step 3 — Alignment  (BWA MEM, {THREADS} threads)")

        if not ref_ok:
            st.error("Reference not found — cannot continue"); wipe_temp(); continue
        if not (TOOLS["bwa"] and TOOLS["samtools"]):
            st.error(f"Missing tools: bwa / samtools"); wipe_temp(); continue

        bam_sorted = t(f"{sample}.sorted.bam")
        flagstat_text = ""
        elapsed_align = 0

        try:
            t0 = time.time()
            # Build input list for BWA (handle paired-end)
            bwa_inputs = [r1_path] + ([r2_path] if is_paired and r2_path else [])

            # Extend run_alignment to accept extra inputs
            log(f"START  BWA MEM piped (threads={THREADS})")
            bwa_e=[]; view_e=[]; sort_e=[]
            def _drain(s, b):
                try: b.append(s.read().decode(errors="replace"))
                except Exception: pass

            bwa_proc = subprocess.Popen(
                ["bwa", "mem", "-t", str(THREADS),
                 "-R", f"@RG\\tID:{sample}\\tSM:{sample}\\tPL:ILLUMINA\\tLB:{sample}\\tPU:{sample}",
                 REFERENCE] + bwa_inputs,
                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            threading.Thread(target=_drain, args=(bwa_proc.stderr, bwa_e), daemon=True).start()

            view_proc = subprocess.Popen(
                ["samtools", "view", "-@", str(THREADS), "-Sb", "-"],
                stdin=bwa_proc.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            bwa_proc.stdout.close()
            threading.Thread(target=_drain, args=(view_proc.stderr, view_e), daemon=True).start()

            sort_proc = subprocess.Popen(
                ["samtools", "sort", "-@", str(THREADS), "-m", SORT_MEM,
                 "-o", bam_sorted, "-"],
                stdin=view_proc.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            view_proc.stdout.close()
            ts2 = threading.Thread(target=_drain, args=(sort_proc.stderr, sort_e), daemon=True)
            ts2.start()

            with st.spinner(f"BWA MEM → view → sort  ({THREADS} threads)…"):
                sort_proc.wait(); view_proc.wait(); bwa_proc.wait()
            ts2.join(timeout=5)

            if bwa_proc.returncode  != 0: raise Exception(f"BWA failed:\n{''.join(bwa_e)[-500:]}")
            if view_proc.returncode != 0: raise Exception(f"samtools view failed:\n{''.join(view_e)[-400:]}")
            if sort_proc.returncode != 0: raise Exception(f"samtools sort failed:\n{''.join(sort_e)[-400:]}")

            run_cmd(["samtools", "index", bam_sorted], "samtools index sorted")
            elapsed_align = int(time.time()-t0)
            safe_remove(r1_path)
            if is_paired and r2_path: safe_remove(r2_path)

            flagstat_text = subprocess.run(
                ["samtools", "flagstat", "-@", str(THREADS), bam_sorted],
                capture_output=True, text=True).stdout

            mapped_pct = "N/A"
            for line in flagstat_text.split("\n"):
                if "mapped (" in line and "primary" not in line:
                    try: mapped_pct = line.split("(")[1].split("%")[0].strip()+"%"
                    except Exception: pass
                    break

            a1,a2,a3 = st.columns(3)
            a1.metric("Alignment Time", f"{elapsed_align}s")
            a2.metric("Reads Mapped",   mapped_pct)
            a3.metric("BAM Size",       f"{mb(bam_sorted):.0f} MB")
            with st.expander("Flagstat details"):
                st.code(flagstat_text, language="bash")
            st.success("Alignment complete")

        except Exception as e:
            st.error(f"Alignment failed: {e}"); wipe_temp(); continue

        # ─────────────────────────────────────────────────────────────
        # STEP 4 — Mark Duplicates (GATK)
        # ─────────────────────────────────────────────────────────────
        st.subheader("🔁 Step 4 — Mark Duplicates  (GATK MarkDuplicates)")
        bam_dedup  = t(f"{sample}.dedup.bam")
        dedup_metrics = t(f"{sample}.dedup_metrics.txt")
        dup_pct    = "N/A"

        if TOOLS["gatk"]:
            try:
                t0 = time.time()
                with st.spinner("Marking duplicates…"):
                    run_cmd([
                        "gatk", "MarkDuplicates",
                        "-I", bam_sorted,
                        "-O", bam_dedup,
                        "-M", dedup_metrics,
                        "--VALIDATION_STRINGENCY", "SILENT",
                        "--CREATE_INDEX", "true",
                    ], "MarkDuplicates")
                safe_remove(bam_sorted, bam_sorted+".bai")

                # Parse duplication rate
                with open(dedup_metrics) as mf:
                    for line in mf:
                        if line.startswith("Unknown") or (line[0].isdigit() if line else False):
                            cols = line.strip().split("\t")
                            if len(cols) > 8:
                                try: dup_pct = f"{float(cols[8])*100:.2f}%"
                                except Exception: pass
                            break
                safe_remove(dedup_metrics)

                d1,d2 = st.columns(2)
                d1.metric("Dedup BAM Size",   f"{mb(bam_dedup):.0f} MB")
                d2.metric("Duplication Rate", dup_pct)
                st.success(f"Mark Duplicates done ({time.time()-t0:.0f}s)")
            except Exception as e:
                st.error(f"MarkDuplicates failed: {e}")
                st.warning("Continuing with sorted BAM (no dedup)")
                bam_dedup = bam_sorted
        else:
            st.warning("GATK not installed — skipping MarkDuplicates")
            bam_dedup = bam_sorted

        # ─────────────────────────────────────────────────────────────
        # STEP 5 — BQSR (Base Quality Score Recalibration)
        # ─────────────────────────────────────────────────────────────
        st.subheader("📊 Step 5 — BQSR  (Base Quality Score Recalibration)")
        bam_bqsr  = t(f"{sample}.bqsr.bam")
        bqsr_done = False

        if TOOLS["gatk"] and bqsr_ready:
            try:
                t0 = time.time()
                recal_table = t(f"{sample}.recal.table")

                with st.spinner("BaseRecalibrator — building recalibration table…"):
                    run_cmd([
                        "gatk", "BaseRecalibrator",
                        "-I", bam_dedup,
                        "-R", REFERENCE,
                        "--known-sites", KNOWN_DBSNP,
                        "--known-sites", KNOWN_MILLS,
                        "--known-sites", KNOWN_1000G,
                        "-O", recal_table,
                    ], "BaseRecalibrator")

                with st.spinner("ApplyBQSR…"):
                    run_cmd([
                        "gatk", "ApplyBQSR",
                        "-I", bam_dedup,
                        "-R", REFERENCE,
                        "--bqsr-recal-file", recal_table,
                        "-O", bam_bqsr,
                    ], "ApplyBQSR")

                safe_remove(bam_dedup, bam_dedup.replace(".bam",".bai"), recal_table)
                bqsr_done = True
                st.success(f"BQSR done ({time.time()-t0:.0f}s)  →  `{os.path.basename(bam_bqsr)}`  ({mb(bam_bqsr):.0f} MB)")

            except Exception as e:
                st.error(f"BQSR failed: {e}")
                st.warning("Continuing with dedup BAM (no BQSR)")
                bam_bqsr = bam_dedup
        else:
            if not TOOLS["gatk"]:
                st.warning("GATK not installed — skipping BQSR")
            else:
                st.warning("Known-sites VCFs missing — skipping BQSR  "
                           "(results will be less accurate)")
            bam_bqsr = bam_dedup

        # Final BAM for download
        final_bam = bam_bqsr

        # ─────────────────────────────────────────────────────────────
        # STEP 6 — HaplotypeCaller (GATK)
        # ─────────────────────────────────────────────────────────────
        st.subheader("🧬 Step 6 — Variant Calling  (GATK HaplotypeCaller)")
        gvcf_t  = t(f"{sample}.g.vcf.gz")
        vcf_raw = t(f"{sample}.vcf.gz")
        total_variants = snps = indels = 0

        if TOOLS["gatk"]:
            try:
                t0 = time.time()
                with st.spinner("HaplotypeCaller — GVCF mode…"):
                    run_cmd([
                        "gatk", "HaplotypeCaller",
                        "-I", final_bam,
                        "-R", REFERENCE,
                        "-O", gvcf_t,
                        "-ERC", "GVCF",
                        "--tmp-dir", "temp",
                        "--native-pair-hmm-threads", str(THREADS),
                    ], "HaplotypeCaller GVCF")

                with st.spinner("GenotypeGVCFs…"):
                    run_cmd([
                        "gatk", "GenotypeGVCFs",
                        "-R", REFERENCE,
                        "-V", gvcf_t,
                        "-O", vcf_raw,
                        "--tmp-dir", "temp",
                    ], "GenotypeGVCFs")

                st.success(f"HaplotypeCaller done ({time.time()-t0:.0f}s)")

            except Exception as e:
                st.error(f"HaplotypeCaller failed: {e}"); wipe_temp(); continue
        else:
            st.error("GATK not installed — cannot call variants"); wipe_temp(); continue

        # ─────────────────────────────────────────────────────────────
        # STEP 7 — Variant Filtering (GATK hard filters)
        # ─────────────────────────────────────────────────────────────
        st.subheader("🔬 Step 7 — Variant Filtering  (GATK hard filters)")
        vcf_snps_t   = t(f"{sample}.snps.vcf.gz")
        vcf_indels_t = t(f"{sample}.indels.vcf.gz")
        vcf_snps_f   = t(f"{sample}.snps.filtered.vcf.gz")
        vcf_indels_f = t(f"{sample}.indels.filtered.vcf.gz")
        vcf_filtered = t(f"{sample}.filtered.vcf.gz")
        ts_filter    = {"pass": 0, "filtered": 0, "pass_pct": "N/A"}

        try:
            t0 = time.time()

            # Select SNPs
            run_cmd(["gatk","SelectVariants","-R",REFERENCE,
                     "-V",vcf_raw,"--select-type-to-include","SNP","-O",vcf_snps_t], "SelectVariants SNPs")
            # Filter SNPs
            run_cmd(["gatk","VariantFiltration","-R",REFERENCE,"-V",vcf_snps_t,
                     "--filter-expression","QD < 2.0",    "--filter-name","QD2",
                     "--filter-expression","FS > 60.0",   "--filter-name","FS60",
                     "--filter-expression","MQ < 40.0",   "--filter-name","MQ40",
                     "--filter-expression","MQRankSum < -12.5","--filter-name","MQRankSum-12.5",
                     "--filter-expression","ReadPosRankSum < -8.0","--filter-name","ReadPosRankSum-8",
                     "-O",vcf_snps_f], "FilterVariants SNPs")

            # Select Indels
            run_cmd(["gatk","SelectVariants","-R",REFERENCE,
                     "-V",vcf_raw,"--select-type-to-include","INDEL","-O",vcf_indels_t], "SelectVariants Indels")
            # Filter Indels
            run_cmd(["gatk","VariantFiltration","-R",REFERENCE,"-V",vcf_indels_t,
                     "--filter-expression","QD < 2.0",    "--filter-name","QD2",
                     "--filter-expression","FS > 200.0",  "--filter-name","FS200",
                     "--filter-expression","ReadPosRankSum < -20.0","--filter-name","ReadPosRankSum-20",
                     "-O",vcf_indels_f], "FilterVariants Indels")

            # Merge SNPs + Indels
            run_cmd(["gatk","MergeVcfs",
                     "-I",vcf_snps_f,"-I",vcf_indels_f,
                     "-O",vcf_filtered], "MergeVcfs")

            safe_remove(vcf_snps_t, vcf_indels_t, vcf_snps_f, vcf_indels_f)

            # Count PASS variants
            count_r = subprocess.run(
                ["bcftools","stats",vcf_filtered] if TOOLS["bcftools"] else
                ["grep","-v","^#",vcf_filtered],
                capture_output=True, text=True)
            if TOOLS["bcftools"]:
                for line in count_r.stdout.split("\n"):
                    if line.startswith("SN") and "number of SNPs" in line:
                        snps = int(line.split("\t")[-1])
                    if line.startswith("SN") and "number of indels" in line:
                        indels = int(line.split("\t")[-1])
            total_variants = snps + indels

            # Count PASS vs FILTER
            pass_count = filt_count = 0
            check = subprocess.run(["bcftools","view","-f","PASS",vcf_filtered],
                                   capture_output=True, text=True)
            for line in check.stdout.split("\n"):
                if not line.startswith("#") and line.strip(): pass_count += 1
            filt_count = total_variants - pass_count
            pass_pct = f"{pass_count/total_variants*100:.1f}%" if total_variants else "N/A"
            ts_filter = {"pass": pass_count, "filtered": filt_count, "pass_pct": pass_pct}

            # Save VCFs to server
            with st.spinner("Saving VCF files to server…"):
                server_copy(gvcf_t,     os.path.join(SERVER_VCF, f"{sample}.g.vcf.gz"))
                server_copy(vcf_raw,    os.path.join(SERVER_VCF, f"{sample}.vcf.gz"))
                server_copy(vcf_filtered,os.path.join(SERVER_VCF, f"{sample}.filtered.vcf.gz"))

            v1,v2,v3,v4,v5 = st.columns(5)
            v1.metric("Total Variants",  f"{total_variants:,}")
            v2.metric("SNPs",            f"{snps:,}")
            v3.metric("Indels",          f"{indels:,}")
            v4.metric("PASS",            f"{pass_count:,}")
            v5.metric("Filtered Out",    f"{filt_count:,}")
            st.success(f"Variant filtering done ({time.time()-t0:.0f}s) — VCFs saved to server")

        except Exception as e:
            st.error(f"Variant filtering failed: {e}"); wipe_temp(); continue

        # ─────────────────────────────────────────────────────────────
        # STEP 8 — QC PDF
        # ─────────────────────────────────────────────────────────────
        st.subheader("📄 Step 8 — QC Report  (MedGenome-style PDF)")
        server_pdf = os.path.join(SERVER_QC, f"{sample}_QC_Report.pdf")
        pdf_bytes  = None
        try:
            with st.spinner("Building QC PDF…"):
                pdf_bytes = build_qc_pdf(
                    sample, qc_raw, qc_trim, flagstat_text, elapsed_align,
                    dup_pct, bqsr_done, total_variants, snps, indels,
                    ts_filter, server_pdf)
            st.success(f"✅ QC PDF auto-saved to server: `{server_pdf}`")
        except Exception as e:
            st.error(f"PDF failed: {e}")

        # ─────────────────────────────────────────────────────────────
        # DOWNLOADS
        # ─────────────────────────────────────────────────────────────
        st.subheader("⬇️ Downloads")

        st.markdown("**🧬 VCF Files** *(auto-saved to server)*")
        filt_bytes = open(vcf_filtered,"rb").read() if os.path.exists(vcf_filtered) else None
        gvcf_bytes = open(gvcf_t,      "rb").read() if os.path.exists(gvcf_t)       else None
        raw_bytes  = open(vcf_raw,     "rb").read() if os.path.exists(vcf_raw)      else None

        if filt_bytes:
            st.download_button(f"⬇️  Filtered VCF — {sample}.filtered.vcf.gz",
                data=filt_bytes, file_name=f"{sample}.filtered.vcf.gz",
                mime="application/gzip", key=f"fvcf_{sample}",
                use_container_width=True, type="primary")
        dv1, dv2 = st.columns(2)
        with dv1:
            if gvcf_bytes:
                st.download_button("⬇️  GVCF (HaplotypeCaller raw)",
                    data=gvcf_bytes, file_name=f"{sample}.g.vcf.gz",
                    mime="application/gzip", key=f"gvcf_{sample}")
        with dv2:
            if raw_bytes:
                st.download_button("⬇️  Genotyped VCF (pre-filter)",
                    data=raw_bytes, file_name=f"{sample}.vcf.gz",
                    mime="application/gzip", key=f"rvcf_{sample}")

        st.markdown("**📄 QC Report** *(auto-saved to server)*")
        if pdf_bytes:
            st.download_button("⬇️  QC Report PDF (MedGenome-style)",
                data=pdf_bytes, file_name=f"{sample}_QC_Report.pdf",
                mime="application/pdf", key=f"pdf_{sample}",
                use_container_width=True)

        st.markdown("**💾 BAM File** *(download only — not saved to server)*")
        st.caption(f"Post-{'BQSR' if bqsr_done else 'dedup'} BAM  "
                   f"({mb(final_bam):.0f} MB) — download if needed before closing.")
        bam_bytes = open(final_bam,      "rb").read() if os.path.exists(final_bam)      else None
        bai_bytes = open(final_bam+".bai","rb").read() if os.path.exists(final_bam+".bai") else None
        b1,b2 = st.columns(2)
        with b1:
            if bam_bytes:
                st.download_button(f"⬇️  BAM  ({mb(final_bam):.0f} MB)",
                    data=bam_bytes, file_name=f"{sample}.bam",
                    mime="application/octet-stream", key=f"bam_{sample}",
                    use_container_width=True)
        with b2:
            if bai_bytes:
                st.download_button("⬇️  BAI index",
                    data=bai_bytes, file_name=f"{sample}.bam.bai",
                    mime="application/octet-stream", key=f"bai_{sample}",
                    use_container_width=True)

        # Wipe temp
        wipe_temp()
        st.success(f"✅  {sample} complete — temp wiped.  QC PDF + VCF auto-saved to server.")

        with st.expander("📁 Server file locations"):
            st.code(
                f"QC Report     : {server_pdf}\n"
                f"GVCF          : {os.path.join(SERVER_VCF, sample+'.g.vcf.gz')}\n"
                f"Genotyped VCF : {os.path.join(SERVER_VCF, sample+'.vcf.gz')}\n"
                f"Filtered VCF  : {os.path.join(SERVER_VCF, sample+'.filtered.vcf.gz')}\n"
                f"BAM           : download-only, not saved to server",
                language="bash")

    # QC Dashboard
    if qc_dashboard:
        st.markdown("---")
        st.header("📊 QC Dashboard — All Samples")
        df = pd.DataFrame(qc_dashboard)
        st.dataframe(df, use_container_width=True, hide_index=True)
        c1,c2 = st.columns(2)
        with c1:
            st.subheader("Reads per Sample")
            st.bar_chart(df.set_index("Sample")["Reads"])
        with c2:
            st.subheader("GC% per Sample")
            st.bar_chart(df.set_index("Sample")["GC"])

    show_server_browser()
    with st.expander("📋 Pipeline Log"):
        if os.path.exists(LOG_FILE):
            with open(LOG_FILE) as fh: st.code(fh.read(), language="bash")

else:
    st.info("📂 Upload FASTQ / FASTQ.GZ files to begin.")
    show_server_browser()
    st.markdown(f"""
---
### Full GATK Best Practices Pipeline

| Step | Tool | What it does |
|------|------|-------------|
| 1 | **FastQC** | Raw read QC — GC%, adapter, duplication |
| 2 | **fastp** | Adapter trimming, quality filtering, Q20/Q30 stats |
| 3 | **BWA MEM** | Align to GRCh38, piped to samtools sort (no SAM on disk) |
| 4 | **GATK MarkDuplicates** | Flag PCR/optical duplicates |
| 5 | **GATK BQSR** | Recalibrate base quality scores using dbSNP + Mills + 1000G |
| 6 | **GATK HaplotypeCaller** | Germline variant calling in GVCF mode |
| 6b | **GATK GenotypeGVCFs** | Genotype the GVCF → raw VCF |
| 7 | **GATK VariantFiltration** | Hard-filter SNPs and Indels separately |
| 8 | **QC PDF** | MedGenome-style report covering all steps |

**Saved to server automatically:** Filtered VCF, GVCF, raw genotyped VCF, QC PDF  
**Download-only (not saved to server):** BAM (post-BQSR)  
**Threads:** {THREADS} of {TOTAL_CPUS} cores
    """)