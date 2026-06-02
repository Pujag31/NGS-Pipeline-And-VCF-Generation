import json
import subprocess
import uuid
from pathlib import Path

from flask import Flask, redirect, render_template, request, send_file, url_for
from werkzeug.utils import secure_filename


BASE_DIR = Path(__file__).resolve().parent
CONFIG_PATH = BASE_DIR / "config.json"
INPUT_DIR = BASE_DIR / "input_fastq"
JOBS_DIR = BASE_DIR / "jobs"
RUNNER = BASE_DIR / "pipeline_runner.py"

app = Flask(__name__)
app.config["MAX_CONTENT_LENGTH"] = 120 * 1024 * 1024 * 1024


def load_config():
    with CONFIG_PATH.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def safe_name(value):
    cleaned = "".join(char if char.isalnum() or char in "._-" else "_" for char in value)
    return cleaned.strip("._-") or "sample"


def list_fastqs():
    INPUT_DIR.mkdir(exist_ok=True)
    files = []
    for pattern in ("*.fastq", "*.fq", "*.fastq.gz", "*.fq.gz"):
        files.extend(INPUT_DIR.glob(pattern))
    return sorted(files, key=lambda path: path.name.lower())


def save_upload(uploaded_file, path):
    chunk_size = 64 * 1024 * 1024
    with path.open("wb") as handle:
        while True:
            chunk = uploaded_file.stream.read(chunk_size)
            if not chunk:
                break
            handle.write(chunk)


def read_text(path, default=""):
    if not path.exists():
        return default
    return path.read_text(encoding="utf-8", errors="replace")


@app.route("/", methods=["GET"])
def index():
    JOBS_DIR.mkdir(exist_ok=True)
    INPUT_DIR.mkdir(exist_ok=True)
    fastqs = list_fastqs()
    jobs = sorted([path for path in JOBS_DIR.iterdir() if path.is_dir()], reverse=True)
    return render_template("index.html", fastqs=fastqs, jobs=jobs, read_text=read_text)


@app.route("/start", methods=["POST"])
def start():
    sample = safe_name(request.form.get("sample", "sample"))
    mode = request.form.get("input_mode", "select")
    job_id = f"{sample}_{uuid.uuid4().hex[:12]}"
    job_dir = JOBS_DIR / job_id
    raw_dir = job_dir / "raw_fastq"
    raw_dir.mkdir(parents=True, exist_ok=True)

    if mode == "upload":
        r1_file = request.files.get("r1_upload")
        r2_file = request.files.get("r2_upload")
        if not r1_file or not r2_file:
            return "Upload both R1 and R2 files.", 400
        r1 = raw_dir / secure_filename(r1_file.filename)
        r2 = raw_dir / secure_filename(r2_file.filename)
        save_upload(r1_file, r1)
        save_upload(r2_file, r2)
    else:
        r1 = INPUT_DIR / request.form["r1_select"]
        r2 = INPUT_DIR / request.form["r2_select"]
        if r1 == r2:
            return "R1 and R2 must be different files.", 400

    command = [
        "python3",
        str(RUNNER),
        "--sample",
        sample,
        "--r1",
        str(r1),
        "--r2",
        str(r2),
        "--job-dir",
        str(job_dir),
        "--config",
        str(CONFIG_PATH),
    ]

    with (job_dir / "launcher.log").open("w", encoding="utf-8") as log:
        subprocess.Popen(
            command,
            cwd=BASE_DIR,
            stdout=log,
            stderr=subprocess.STDOUT,
            start_new_session=True,
        )

    return redirect(url_for("job", job_id=job_id))


@app.route("/job/<job_id>")
def job(job_id):
    job_dir = JOBS_DIR / secure_filename(job_id)
    if not job_dir.exists():
        return "Job not found", 404

    outputs = sorted([path for path in job_dir.rglob("*") if path.is_file()])
    return render_template("job.html", job_dir=job_dir, outputs=outputs, read_text=read_text)


@app.route("/download/<job_id>/<path:filename>")
def download(job_id, filename):
    job_dir = JOBS_DIR / secure_filename(job_id)
    target = (job_dir / filename).resolve()
    if job_dir.resolve() not in target.parents and target != job_dir.resolve():
        return "Invalid path", 400
    if not target.exists():
        return "File not found", 404
    return send_file(target, as_attachment=True)


if __name__ == "__main__":
    load_config()
    INPUT_DIR.mkdir(exist_ok=True)
    JOBS_DIR.mkdir(exist_ok=True)
    app.run(host="0.0.0.0", port=5000, debug=False)
