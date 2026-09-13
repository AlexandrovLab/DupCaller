"""--learnOnly must stop do_call right after error-rate estimation,
before any variant-calling work (or its SBS/INDEL/DBS output dirs) starts."""
from types import SimpleNamespace

import numpy as np
import pysam
import pytest

from DupCaller_sub import Caller
from DupCaller_sub.funcs.learn import NUM_BQ

# do_call's num2trinuc covers 2 (4 minus x 2 ref x 4 plus) x 2 loops = 64
# trinucleotide contexts (C/T- and G/A-centered), not the full 96 SBS
# substitution-type x context count.
NUM_TRINUC = 64


def _write_bam(path):
    header = {"SQ": [{"SN": "chr1", "LN": 1000}]}
    with pysam.AlignmentFile(str(path), "wb", header=header):
        pass


def _base_args(tmp_path, bam_path, learn_only):
    return SimpleNamespace(
        bam=str(bam_path),
        reference="unused.fa",
        output=str(tmp_path / "sample"),
        regions=["chr1"],
        regionfile=None,
        regionst=None,
        threads=1,
        normalBams=None,
        germline=None,
        germlineAfCutoff=0.001,
        noise=None,
        indelbed=None,
        trimF=0,
        trimR=0,
        minNdepth=1,
        minBq=18,
        maxAF=1.0,
        nmflt=5,
        windowSize=100000,
        minMeanASXS=50,
        naf=0.05,
        rescue=False,
        maxZeroQualFrac=1.0,
        maxPileupDepth=100000,
        mapq=40,
        barcode="DB,1,-",
        seed=1,
        errprefix=None,
        muterateprefix=None,
        minAltQual=0,
        minRef=0,
        minAlt=0,
        pseudocount=0.5,
        learnOnly=learn_only,
    )


def _learn_stage_result():
    mismatch_profile = np.zeros((NUM_TRINUC, 4))
    hp_alt_profile = np.zeros((10, 12), dtype=int)
    str_alt_profile = np.zeros((5, 11), dtype=int)
    mismatch_dmg_profile = np.zeros((NUM_TRINUC, 4))
    hp_dmg_profile = np.zeros((10, 12), dtype=int)
    str_dmg_profile = np.zeros((5, 11), dtype=int)
    sbs_alt_bq_hist = np.zeros((NUM_TRINUC, 4, NUM_BQ), dtype=int)
    return (
        mismatch_profile,
        hp_alt_profile,
        str_alt_profile,
        mismatch_dmg_profile,
        hp_dmg_profile,
        str_dmg_profile,
        sbs_alt_bq_hist,
    )


def test_learn_only_stops_before_calling(tmp_path, monkeypatch):
    bam_path = tmp_path / "in.bam"
    _write_bam(bam_path)
    calls = []

    def fake_callBam(params, nn):
        calls.append(params.get("isLearn"))
        return _learn_stage_result()

    monkeypatch.setattr(Caller, "check_input_files_exist", lambda args: None)
    monkeypatch.setattr(Caller, "callBam", fake_callBam)

    args = _base_args(tmp_path, bam_path, learn_only=True)
    Caller.do_call(args)

    # Only the learn-stage callBam ran; calling was never dispatched.
    assert calls == [True]
    error_dir = tmp_path / "sample" / "ERROR"
    assert (error_dir / "sample.amp.tn.srd.txt").exists()
    assert not (tmp_path / "sample" / "SBS").exists()


def test_without_learn_only_calling_is_attempted(tmp_path, monkeypatch):
    bam_path = tmp_path / "in.bam"
    _write_bam(bam_path)
    calls = []

    def fake_callBam(params, nn):
        calls.append(params.get("isLearn"))
        return _learn_stage_result()

    monkeypatch.setattr(Caller, "check_input_files_exist", lambda args: None)
    monkeypatch.setattr(Caller, "callBam", fake_callBam)

    args = _base_args(tmp_path, bam_path, learn_only=False)
    # The calling round expects a much longer result tuple than the
    # learn-shaped stub above -- unpacking it failing is exactly the
    # signal that do_call went on to dispatch the calling round.
    with pytest.raises(ValueError):
        Caller.do_call(args)

    assert calls == [True, None]
