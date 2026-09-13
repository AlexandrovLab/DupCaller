"""MC validation must fail before Caller creates outputs or starts work."""
from types import SimpleNamespace

import pysam
import pytest

from DupCaller_sub import Caller


def _write_bam(path, missing_flag=None, count=3, missing_index=2):
    header = {"SQ": [{"SN": "chr1", "LN": 1000000}]}
    with pysam.AlignmentFile(str(path), "wb", header=header) as bam:
        for index in range(count):
            rec = pysam.AlignedSegment(bam.header)
            rec.query_name = f"read{index}"
            rec.flag = 3
            rec.query_sequence = "A" * 10
            rec.reference_id = 0
            rec.reference_start = index * 20
            rec.cigarstring = "10M"
            rec.set_tag("MC", "10M")
            if index == missing_index and missing_flag is not None:
                rec.set_tag("MC", None)
                rec.flag = missing_flag
                if rec.is_unmapped:
                    rec.reference_id = -1
                    rec.reference_start = -1
                    rec.cigarstring = None
            bam.write(rec)


def test_accepts_fully_tagged_bam(tmp_path):
    path = tmp_path / "tagged.bam"
    _write_bam(path)
    Caller.check_mate_cigar_tags(str(path))


@pytest.mark.parametrize("flag", [3, 99, 147])
def test_missing_mc_stops_caller_before_creating_outputs(
    tmp_path, monkeypatch, capsys, flag
):
    path = tmp_path / "missing.bam"
    _write_bam(path, missing_flag=flag)
    output = tmp_path / "output"
    # threads > 1: only multi-threaded runs hand reads off across chunk
    # boundaries via MC tags, so only they need the guard (see -p1 test
    # below).
    args = SimpleNamespace(bam=str(path), reference=None, output=str(output), threads=2)
    monkeypatch.setattr(Caller, "check_input_files_exist", lambda args: None)

    with pytest.raises(SystemExit) as exc:
        Caller.do_call(args)

    assert exc.value.code == 1
    message = capsys.readouterr().err
    assert "read2" in message
    assert str(path) in message
    assert "MC (mate CIGAR)" in message
    assert not output.exists()


def test_single_threaded_call_skips_mc_check(monkeypatch, tmp_path):
    """-p 1 has no chunk boundaries to reconcile, so it must not be
    blocked by a BAM that legitimately lacks MC tags."""
    calls = []
    monkeypatch.setattr(Caller, "check_input_files_exist", lambda args: None)
    monkeypatch.setattr(
        Caller, "check_mate_cigar_tags", lambda *a, **k: calls.append((a, k))
    )
    args = SimpleNamespace(
        bam=str(tmp_path / "missing.bam"),
        reference=None,
        output=str(tmp_path / "output"),
        threads=1,
    )

    with pytest.raises(AttributeError):
        # do_call proceeds past the (skipped) MC check into seed
        # resolution, which this minimal fixture doesn't provide -- we
        # only care that the MC check itself was never invoked.
        Caller.do_call(args)

    assert calls == []


def test_multi_threaded_call_runs_mc_check(monkeypatch, tmp_path):
    calls = []
    monkeypatch.setattr(Caller, "check_input_files_exist", lambda args: None)
    monkeypatch.setattr(
        Caller, "check_mate_cigar_tags", lambda *a, **k: calls.append((a, k))
    )
    args = SimpleNamespace(
        bam=str(tmp_path / "missing.bam"),
        reference=None,
        output=str(tmp_path / "output"),
        threads=2,
    )

    with pytest.raises(AttributeError):
        Caller.do_call(args)

    assert len(calls) == 1


@pytest.mark.parametrize("flag", [0, 4, 9, 256, 2048])
def test_ignores_missing_mc_without_proper_pair_flag(tmp_path, flag):
    path = tmp_path / "not_proper.bam"
    _write_bam(path, missing_flag=flag)
    Caller.check_mate_cigar_tags(str(path))


@pytest.mark.parametrize("flag", [3 | 256, 3 | 2048])
def test_ignores_missing_mc_on_secondary_or_supplementary_proper_pair(tmp_path, flag):
    """Matches _eligible_overflow_read: a proper-pair read that is also
    secondary/supplementary is never the overflow machinery's concern, so
    a missing MC tag on it must not abort the run."""
    path = tmp_path / "secondary_supplementary.bam"
    _write_bam(path, missing_flag=flag)
    Caller.check_mate_cigar_tags(str(path))


def test_scans_entire_file_not_just_first_10000_reads(tmp_path):
    """The guard must not sample only a prefix -- a merged BAM can have
    MC tags on an early lane/read-group and not a later one."""
    path = tmp_path / "large.bam"
    _write_bam(path, missing_flag=3, count=10001, missing_index=10000)
    with pytest.raises(SystemExit):
        Caller.check_mate_cigar_tags(str(path))
