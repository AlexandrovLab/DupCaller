"""--doNotOutputMaskedSite: estimate reads the flag from _call_params.log."""

from DupCaller_sub.Estimate import call_flag_set, call_skipped_coverage_pass


def _write_params(tmp_path, lines):
    (tmp_path / "s_call_params.log").write_text(
        "DupCaller call — parameter log\n"
        "Command:  DupCaller.py call --doNotOutputMaskedSite\n"
        "\n--- All parameters (resolved values) ---\n" + "".join(lines)
    )


def test_call_flag_set_reads_resolved_values(tmp_path):
    _write_params(
        tmp_path,
        ["  doNotOutputMaskedSite: True\n", "  skipCoveragePass: False\n"],
    )
    assert call_flag_set(str(tmp_path), "s", "doNotOutputMaskedSite")
    assert not call_skipped_coverage_pass(str(tmp_path), "s")


def test_call_flag_unset_or_missing(tmp_path):
    _write_params(tmp_path, ["  doNotOutputMaskedSite: False\n"])
    assert not call_flag_set(str(tmp_path), "s", "doNotOutputMaskedSite")
    # A call run from before the option has no line for it.
    _write_params(tmp_path, ["  skipCoveragePass: True\n"])
    assert not call_flag_set(str(tmp_path), "s", "doNotOutputMaskedSite")
    assert call_skipped_coverage_pass(str(tmp_path), "s")
    assert not call_flag_set(str(tmp_path / "absent"), "s", "doNotOutputMaskedSite")
