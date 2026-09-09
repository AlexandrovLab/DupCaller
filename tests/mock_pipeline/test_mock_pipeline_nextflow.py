"""Install-time regression test for the containerized Nextflow pipeline
(nextflow/DupCaller.nf): index -> trim -> bwa mem -> gatk MarkDuplicates ->
call -> estimate, run through Nextflow against the published
yuhecheng62/dupcaller image instead of DupCaller.py directly on PATH (see
test_mock_pipeline.py for that variant).

Skipped automatically unless nextflow, bwa, and samtools are on PATH and
either a reachable docker daemon or singularity is available.
"""
import shutil
import subprocess
from pathlib import Path

import pytest

MOCK_DIR = Path(__file__).resolve().parent


def _container_profile():
    """Returns 'docker' or 'singularity' -- whichever is actually usable --
    or None if neither is."""
    if shutil.which("docker") is not None:
        try:
            subprocess.run(
                ["docker", "info"],
                check=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                timeout=15,
            )
            return "docker"
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired):
            pass
    if shutil.which("singularity") is not None:
        return "singularity"
    return None


def _missing_tools():
    missing = [t for t in ("nextflow", "bwa", "samtools") if shutil.which(t) is None]
    if _container_profile() is None:
        missing.append("docker (running) or singularity")
    return missing


@pytest.fixture(scope="module")
def pipeline_output(tmp_path_factory):
    missing = _missing_tools()
    if missing:
        pytest.skip(f"required tool(s) not available: {', '.join(missing)}")

    outdir = tmp_path_factory.mktemp("mock_pipeline_nextflow_run")
    profile = _container_profile()
    subprocess.run(
        [
            "bash",
            str(MOCK_DIR / "run_nextflow_pipeline.sh"),
            str(outdir),
            profile,
        ],
        check=True,
        timeout=1200,
    )
    return outdir


@pytest.mark.parametrize(
    "rel_path",
    [
        "results/mock/SBS/mock_sbs.vcf",
        "results/mock/INDEL/mock_indel.vcf",
        "results/mock/DBS/mock_dbs.vcf",
        "results/mock/SBS/mock_sbs_burden.txt",
        "results/mock/INDEL/mock_indel_burden.txt",
        "results/mock/mock_coverage.bed.gz",
    ],
)
def test_output_file_produced(pipeline_output, rel_path):
    path = pipeline_output / rel_path
    assert path.exists(), f"nextflow pipeline did not produce {rel_path}"
    assert path.stat().st_size > 0, f"{rel_path} is empty"
