from pathlib import Path

from reportlab.lib.pagesizes import A4
from reportlab.pdfgen import canvas


def file_size(path):
    path = Path(path)
    if not path.exists():
        return "missing"
    size = path.stat().st_size
    for unit in ("B", "KB", "MB", "GB"):
        if size < 1024:
            return f"{size:.2f} {unit}"
        size /= 1024
    return f"{size:.2f} TB"


def write_pdf_report(sample, job_dir, vcf, bam, metrics, log_path, output_pdf):
    output_pdf.parent.mkdir(parents=True, exist_ok=True)
    pdf = canvas.Canvas(str(output_pdf), pagesize=A4)
    width, height = A4
    y = height - 50

    def line(text):
        nonlocal y
        if y < 60:
            pdf.showPage()
            y = height - 50
        pdf.drawString(50, y, str(text)[:110])
        y -= 18

    pdf.setFont("Helvetica-Bold", 16)
    line("NGS Variant Calling Report")
    pdf.setFont("Helvetica", 11)
    line(f"Sample: {sample}")
    line(f"Job folder: {job_dir}")
    line("")
    line("Outputs")
    line(f"VCF: {vcf.name} ({file_size(vcf)})")
    line(f"BAM: {bam.name} ({file_size(bam)})")
    line(f"Duplicate metrics: {metrics.name} ({file_size(metrics)})")
    line("")
    line("Pipeline")
    for step in (
        "FastQC quality check",
        "fastp trimming",
        "BWA-MEM alignment to GRCh38",
        "samtools BAM sorting and indexing",
        "GATK MarkDuplicates",
        "GATK BaseRecalibrator and ApplyBQSR if known-sites are configured",
        "GATK HaplotypeCaller VCF generation",
    ):
        line(f"- {step}")
    line("")
    line("Recent Log Lines")
    if Path(log_path).exists():
        lines = Path(log_path).read_text(encoding="utf-8", errors="replace").splitlines()[-18:]
        for item in lines:
            line(item)

    pdf.save()
