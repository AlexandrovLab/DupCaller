"""Install-time regression test for the containerized Nextflow pipeline
(nextflow/DupCaller.nf): DupCaller index + bwa-mem2 index -> trim ->
bwa-mem2 mem -> gatk MarkDuplicates -> call -> estimate, run through
Nextflow against the published yuhecheng62/dupcaller image instead of
DupCaller.py directly on PATH (see test_mock_pipeline.py for that variant).

Runs every case of run_nextflow_pipeline.sh: the default configuration
with every optional resource unset, then one run per optional input
(germline VCF, noise masks, target/indel/gene BED, estimate_dilute, shared
normal BAM), each checked for the DupCaller.py option it should (or, by
default, should not) produce.

The default case's deterministic outputs (VCFs, burdens, spectra, rate
tables, error profiles, stats) must also match expected_matched_normal/,
the same files the direct-CLI regression (test_mock_pipeline.py,
MATCHED_NORMAL=1 run) is held to: a container image or workflow change that
alters calls fails here.

Skipped automatically unless nextflow is on PATH and either a reachable
docker daemon or singularity is available.
"""
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

MOCK_DIR = Path(__file__).resolve().parent
EXPECTED_MATCHED_DIR = MOCK_DIR / "expected_matched_normal"

sys.path.insert(0, str(MOCK_DIR))
from compare_outputs import compare_file  # noqa: E402

CASES = [
    "default",
    "germline",
    "noise_mask",
    "target_bed",
    "indel_bed",
    "gene_bed",
    "dilute",
    "normal_bam",
]
OUTPUTS = [
    "mock/SBS/mock_sbs.vcf",
    "mock/INDEL/mock_indel.vcf",
    "mock/DBS/mock_dbs.vcf",
    "mock/SBS/mock_sbs_burden.txt",
    "mock/INDEL/mock_indel_burden.txt",
    "mock/mock_coverage.bed.gz",
    "mock_estimate_params.log",
]


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
    missing = [t for t in ("nextflow",) if shutil.which(t) is None]
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
            *CASES,
        ],
        check=True,
        timeout=3600,
    )
    return outdir / "results"


def _params(path):
    """{key: value} from a DupCaller *_params.log's resolved-parameter list."""
    out = {}
    for line in path.read_text().splitlines():
        if line.startswith("  ") and ": " in line:
            key, value = line.strip().split(": ", 1)
            out[key] = value
    return out


@pytest.mark.parametrize("case", CASES)
@pytest.mark.parametrize("rel_path", OUTPUTS)
def test_output_file_produced(pipeline_output, case, rel_path):
    path = pipeline_output / case / rel_path
    assert path.exists(), f"case {case}: nextflow pipeline did not produce {rel_path}"
    assert path.stat().st_size > 0, f"case {case}: {rel_path} is empty"


def test_default_passes_no_optional_resource(pipeline_output):
    call = _params(pipeline_output / "default/mock/mock_call_params.log")
    est = _params(pipeline_output / "default/mock_estimate_params.log")
    assert call["germline"] == "None"
    assert call["noise"] == "None"
    assert call["regionfile"] == "None"
    assert call["indelbed"] == "False"
    assert est["genebed"] == "None"
    assert est["dilute"] == "False"


@pytest.mark.parametrize(
    "case, log, key, expected",
    [
        ("germline", "call", "germline", "germline.vcf.gz"),
        ("noise_mask", "call", "noise", "['snp_mask.bed.gz', 'noise_mask.bed.gz']"),
        ("target_bed", "call", "regionfile", "target.bed.gz"),
        ("indel_bed", "call", "indelbed", "indel_pon.bed.gz"),
        ("gene_bed", "estimate", "genebed", "genes.bed.gz"),
        ("dilute", "estimate", "dilute", "True"),
        ("normal_bam", "call", "normalBams", "['mock_normal.mkdped.bam']"),
    ],
)
def test_optional_resource_reaches_dupcaller(pipeline_output, case, log, key, expected):
    path = (
        pipeline_output / case / "mock/mock_call_params.log"
        if log == "call"
        else pipeline_output / case / "mock_estimate_params.log"
    )
    assert _params(path)[key] == expected


def test_gene_bed_writes_gene_coverage(pipeline_output):
    path = pipeline_output / "gene_bed/mock/mock_gene_coverage.txt"
    assert path.exists() and path.stat().st_size > 0


def test_estimate_leaves_call_output_unchanged(pipeline_output):
    """ESTIMATE_BURDEN must not write into CALL_VARIANTS' cached output:
    the gene_bed case runs before dilute and shares its CALL_VARIANTS task
    via -resume, so its gene coverage file must not show up in dilute's
    results."""
    assert not (pipeline_output / "dilute/mock/mock_gene_coverage.txt").exists()


@pytest.mark.parametrize("case", CASES)
def test_stats_has_one_coverage_block(pipeline_output, case):
    """estimate appends base coverage to _stats.txt. Cases sharing one cached
    CALL_VARIANTS task via -resume must each publish exactly one block, i.e.
    no estimate run appended through to the call output."""
    lines = (pipeline_output / case / "mock/mock_stats.txt").read_text().splitlines()
    for key in ("SBS Base Coverage", "Indel Base Coverage", "DBS Base Coverage"):
        assert sum(line.startswith(key + "\t") for line in lines) == 1, key


@pytest.mark.parametrize(
    "rel_path",
    sorted(
        str(p.relative_to(EXPECTED_MATCHED_DIR))
        for p in EXPECTED_MATCHED_DIR.rglob("*")
        if p.is_file()
    ),
)
def test_default_case_matches_direct_cli_expected(pipeline_output, rel_path):
    actual_path = pipeline_output / "default" / rel_path
    assert actual_path.exists(), f"nextflow pipeline did not produce {rel_path}"
    diffs = compare_file(EXPECTED_MATCHED_DIR / rel_path, actual_path)
    assert not diffs, "\n".join(diffs)
