from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.reports import write_alignment_filtering_methods_summary_text


def fixture(name: str) -> Path:
    direct = Path(__file__).parent / "fixtures" / "alignments" / name
    if direct.exists():
        return direct
    return Path(__file__).parent / "fixtures" / "metadata" / name


def test_write_alignment_filtering_methods_summary_text_reports_profile_and_retention(
    tmp_path: Path,
) -> None:
    output_path = tmp_path / "alignment-filtering-methods-summary.md"

    result = write_alignment_filtering_methods_summary_text(
        output_path,
        alignment_path=fixture("example_alignment_filtering.fasta"),
        profile_name="moderate",
        group_table_path=fixture("example_alignment_groups.tsv"),
        group_columns=["region"],
    )

    assert result.output_path == output_path
    assert result.removed_site_count == 4
    assert result.removed_sequence_count == 1
    assert result.retained_sequence_count == 3
    assert result.retained_alignment_length == 8
    assert "Alignment Filtering Methods Summary" in result.text
    assert "- named profile: `moderate`" in result.text
    assert "- site missingness threshold: `0.5`" in result.text
    assert "- sequence missingness threshold: `0.4`" in result.text
    assert "- removed sites: `4` (missingness-threshold=4)" in result.text
    assert "- removed sequence identities: `D` (missingness-threshold)" in result.text
    assert "- retained alignment length: `8` of `12`" in result.text
    assert (
        "group retention `region=island`: original `1`, retained `0`, removed `1`, removed fraction `1`"
        in result.text
    )
    assert (
        "cleaning removed most taxa from one or more metadata or trait groups"
        in result.text
    )
    assert output_path.read_text(encoding="utf-8") == result.text


def test_write_alignment_filtering_methods_summary_text_reports_no_removals(
    tmp_path: Path,
) -> None:
    output_path = tmp_path / "alignment-filtering-methods-summary-clean.md"

    result = write_alignment_filtering_methods_summary_text(
        output_path,
        alignment_path=fixture("example_alignment.fasta"),
        profile_name="phylogenomics-scale",
    )

    assert result.removed_site_count == 0
    assert result.removed_sequence_count == 0
    assert "- removed sites: `0`" in result.text
    assert "- removed sequences: `0`" in result.text
    assert "no metadata or trait group-retention audit was requested" in result.text
