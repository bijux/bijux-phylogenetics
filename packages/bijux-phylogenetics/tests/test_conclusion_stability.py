from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.diagnostics.conclusion_stability import (
    AncestralStateStabilityRow,
    ComparativeCoefficientStabilityRow,
    ConclusionStabilityConclusionRow,
    KeyCladeStabilityRow,
    SupportValueStabilityRow,
    build_conclusion_stability_report,
    write_conclusion_stability_report_html,
    write_conclusion_stability_summary_table,
)


def test_build_conclusion_stability_report_counts_classes() -> None:
    report = build_conclusion_stability_report(
        key_clade_rows=[
            KeyCladeStabilityRow(
                clade_id="bat|canid",
                descendant_taxa=("bat", "canid"),
                bootstrap_frequency=0.95,
                method_presence_fraction=0.9,
                combined_score=0.925,
                stability_class="stable",
                evidence="bootstrap_frequency=0.95; method_presence_fraction=0.9",
            )
        ],
        support_value_rows=[
            SupportValueStabilityRow(
                clade_id="canid|fox",
                descendant_taxa=("canid", "fox"),
                baseline_support_fraction=0.9,
                bootstrap_frequency=0.7,
                method_presence_fraction=0.8,
                mean_method_support_fraction=0.85,
                maximum_support_delta=0.2,
                method_support_consistency=0.64,
                combined_score=0.67,
                stability_class="weak",
                evidence="bootstrap_frequency=0.7; method_support_consistency=0.64",
            )
        ],
        ancestral_state_rows=[
            AncestralStateStabilityRow(
                trait="host_group",
                clade_id="bat|canid|fox",
                descendant_taxa=("bat", "canid", "fox"),
                baseline_state="bat",
                bootstrap_dominant_state="canid",
                bootstrap_state_consistency=0.25,
                method_dominant_state="canid",
                method_state_consistency=0.25,
                combined_score=0.25,
                stability_class="unstable",
                evidence="baseline_state=bat; bootstrap_state_consistency=0.25; method_state_consistency=0.25",
            )
        ],
        comparative_coefficient_rows=[
            ComparativeCoefficientStabilityRow(
                term="host_group[canid]",
                baseline_estimate=12.0,
                baseline_direction="positive",
                baseline_significant=True,
                bootstrap_direction_consistency=1.0,
                bootstrap_significance_consistency=1.0,
                method_direction_consistency=1.0,
                method_significance_consistency=1.0,
                combined_score=1.0,
                stability_class="stable",
                evidence="baseline_direction=positive; bootstrap_direction_consistency=1; method_direction_consistency=1",
            )
        ],
    )

    assert report.summary.stable_count == 2
    assert report.summary.weak_count == 1
    assert report.summary.unstable_count == 1
    assert report.summary.key_clade_count == 1
    assert report.summary.support_value_count == 1
    assert report.summary.ancestral_state_count == 1
    assert report.summary.comparative_coefficient_count == 1
    assert any(
        row
        == ConclusionStabilityConclusionRow(
            category="comparative_coefficient",
            conclusion_id="host_group[canid]",
            label="host_group[canid]",
            combined_score=1.0,
            stability_class="stable",
            evidence="baseline_direction=positive; bootstrap_direction_consistency=1; method_direction_consistency=1",
        )
        for row in report.conclusion_rows
    )


def test_write_conclusion_stability_outputs_separate_classes(tmp_path: Path) -> None:
    report = build_conclusion_stability_report(
        key_clade_rows=[
            KeyCladeStabilityRow(
                clade_id="stable-clade",
                descendant_taxa=("a", "b"),
                bootstrap_frequency=0.9,
                method_presence_fraction=0.9,
                combined_score=0.9,
                stability_class="stable",
                evidence="stable evidence",
            )
        ],
        support_value_rows=[
            SupportValueStabilityRow(
                clade_id="weak-clade",
                descendant_taxa=("a", "c"),
                baseline_support_fraction=0.8,
                bootstrap_frequency=0.6,
                method_presence_fraction=0.7,
                mean_method_support_fraction=0.75,
                maximum_support_delta=0.25,
                method_support_consistency=0.6,
                combined_score=0.6,
                stability_class="weak",
                evidence="weak evidence",
            )
        ],
        ancestral_state_rows=[
            AncestralStateStabilityRow(
                trait="region_group",
                clade_id="unstable-clade",
                descendant_taxa=("a", "d"),
                baseline_state="north",
                bootstrap_dominant_state="south",
                bootstrap_state_consistency=0.1,
                method_dominant_state="south",
                method_state_consistency=0.2,
                combined_score=0.15,
                stability_class="unstable",
                evidence="unstable evidence",
            )
        ],
        comparative_coefficient_rows=[],
    )

    summary_path = write_conclusion_stability_summary_table(
        tmp_path / "conclusion-stability-summary.tsv",
        report,
    )
    html_path = write_conclusion_stability_report_html(
        tmp_path / "conclusion-stability-report.html",
        report,
    )

    summary_text = summary_path.read_text(encoding="utf-8")
    html_text = html_path.read_text(encoding="utf-8")
    embedded_manifest = html_text.split(
        '<script id="bijux-report-manifest" type="application/json">',
        1,
    )[1].split("</script>", 1)[0]
    payload = json.loads(embedded_manifest)

    assert "stable_count\tweak_count\tunstable_count" in summary_text
    assert "stable-conclusions" in html_text
    assert "weak-conclusions" in html_text
    assert "unstable-conclusions" in html_text
    assert payload["summary"] == {
        "stable_count": 1,
        "unstable_count": 1,
        "weak_count": 1,
    }
