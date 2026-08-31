"""Install-time regression test for the full DupCaller pipeline.

Runs index -> trim -> bwa mem -> gatk MarkDuplicates -> call -> estimate
against the synthetic reference/read set in data/ and diffs every
deterministic output file against the premade results in expected/.

Requires DupCaller.py, bwa, samtools, and gatk; skipped automatically if
any of those aren't on PATH.
"""
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

MOCK_DIR = Path(__file__).resolve().parent
REPO_ROOT = MOCK_DIR.parent.parent
EXPECTED_DIR = MOCK_DIR / "expected"

sys.path.insert(0, str(MOCK_DIR))
from compare_outputs import compare_file  # noqa: E402


def _resolve_dupcaller():
    on_path = shutil.which("DupCaller.py")
    if on_path:
        return [on_path]
    fallback = REPO_ROOT / "src" / "DupCaller.py"
    if fallback.exists():
        return [sys.executable, str(fallback)]
    return None


def _missing_tools():
    missing = []
    if _resolve_dupcaller() is None:
        missing.append("DupCaller.py")
    for tool in ("bwa", "samtools", "gatk"):
        if shutil.which(tool) is None:
            missing.append(tool)
    return missing


@pytest.fixture(scope="module")
def pipeline_output(tmp_path_factory):
    missing = _missing_tools()
    if missing:
        pytest.skip(f"required tool(s) not on PATH: {', '.join(missing)}")

    outdir = tmp_path_factory.mktemp("mock_pipeline_run")
    env = dict(os.environ)
    dupcaller_cmd = _resolve_dupcaller()
    if len(dupcaller_cmd) == 1:
        env["DUPCALLER"] = dupcaller_cmd[0]
    else:
        # Wrap the "python path/to/DupCaller.py" fallback in a shim script
        # so run_pipeline.sh can invoke $DUPCALLER as a single command.
        shim = outdir / "DupCaller.py"
        shim.write_text(
            f'#!/usr/bin/env bash\nexec "{sys.executable}" "{dupcaller_cmd[1]}" "$@"\n'
        )
        shim.chmod(0o755)
        env["DUPCALLER"] = str(shim)

    subprocess.run(
        ["bash", str(MOCK_DIR / "run_pipeline.sh"), str(outdir)],
        check=True,
        env=env,
    )
    return outdir


def _expected_files():
    return sorted(p for p in EXPECTED_DIR.rglob("*") if p.is_file())


@pytest.mark.parametrize(
    "rel_path",
    [str(p.relative_to(EXPECTED_DIR)) for p in _expected_files()],
)
def test_output_matches_expected(pipeline_output, rel_path):
    expected_path = EXPECTED_DIR / rel_path
    actual_path = pipeline_output / rel_path
    assert actual_path.exists(), f"pipeline did not produce {rel_path}"
    diffs = compare_file(expected_path, actual_path)
    assert not diffs, "\n".join(diffs)
