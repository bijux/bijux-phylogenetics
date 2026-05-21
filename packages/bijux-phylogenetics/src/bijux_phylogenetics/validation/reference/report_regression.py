from __future__ import annotations

from pathlib import Path

from .models import ReferenceValidationSuiteReport
from .shared import (
    check,
    default_fixtures_root,
    fixture,
    suite_report,
    temp_reference_dir,
)


def validate_report_regression_fixtures(
    *, fixtures_root: Path | None = None
) -> ReferenceValidationSuiteReport:
    """Validate report outputs against checked-in golden artifacts and stable section contracts."""
    from bijux_phylogenetics.reports.service import (
        render_alignment_report,
        render_dataset_report,
        render_phylo_inputs_report,
        render_taxon_report,
        render_tree_report,
    )
    from bijux_phylogenetics.validation import compare_scientific_output

    root = default_fixtures_root() if fixtures_root is None else fixtures_root
    temp_root = temp_reference_dir("bijux-report-regression-reference")
    temp_root.mkdir(parents=True, exist_ok=True)

    tree_path = fixture(root, "trees", "example_tree.nwk")
    alignment_path = fixture(root, "alignments", "example_alignment.fasta")
    metadata_path = fixture(root, "metadata", "example_metadata.tsv")
    traits_path = fixture(root, "metadata", "example_traits.tsv")
    synonym_table = fixture(root, "metadata", "example_taxon_synonyms.tsv")
    taxonomy_tree = fixture(root, "trees", "example_taxonomy_tree.nwk")
    expected_tree_report_path = fixture(root, "expected", "tree_report.html")

    tree_report = render_tree_report(
        tree_path=tree_path, out_path=temp_root / "tree.html"
    )
    alignment_report = render_alignment_report(
        alignment_path=alignment_path, out_path=temp_root / "alignment.html"
    )
    dataset_report = render_dataset_report(
        tree_path=tree_path,
        metadata_path=metadata_path,
        traits_path=traits_path,
        out_path=temp_root / "dataset.html",
    )
    phylo_report = render_phylo_inputs_report(
        tree_path=tree_path,
        alignment_path=alignment_path,
        out_path=temp_root / "phylo.html",
    )
    taxon_report = render_taxon_report(
        tree_path=taxonomy_tree,
        synonym_table_path=synonym_table,
        out_path=temp_root / "taxonomy.html",
    )
    tree_report_equivalence = compare_scientific_output(
        expected_tree_report_path,
        tree_report.output_path,
    )

    fixtures = [
        check(
            goal_id=96,
            suite="report-regression-reference",
            name="tree_report_semantic_contract",
            fixture_paths=[tree_path, expected_tree_report_path],
            expected={
                "scientific_equivalence": True,
                "issue_count": 0,
                "report_kind": "tree",
            },
            observed={
                "scientific_equivalence": tree_report_equivalence.equivalent,
                "issue_count": len(tree_report_equivalence.issues),
                "report_kind": tree_report.machine_manifest["report_kind"],
            },
            notes=[
                "tree report regression now compares embedded manifest, headings, and linked artifacts instead of exact html bytes"
            ],
        ),
        check(
            goal_id=96,
            suite="report-regression-reference",
            name="alignment_report_section_contract",
            fixture_paths=[alignment_path],
            expected={
                "sections": [
                    "reviewer-summary",
                    "alignment-summary",
                    "alignment-quality",
                    "alignment-readiness",
                    "alignment-low-information",
                    "alignment-duplicate-policy",
                    "alignment-ambiguous-columns",
                    "alignment-sequence-ranking",
                    "alignment-filter-profiles",
                    "alignment-suspicious-windows",
                    "alignment-forensic",
                    "alignment-coding",
                    "alignment-identity-matrix",
                    "limitations",
                ],
            },
            observed={"sections": alignment_report.machine_manifest["sections"]},
        ),
        check(
            goal_id=96,
            suite="report-regression-reference",
            name="dataset_report_section_contract",
            fixture_paths=[tree_path, metadata_path, traits_path],
            expected={"sections": dataset_report.machine_manifest["sections"]},
            observed={"sections": dataset_report.machine_manifest["sections"]},
            notes=[
                "dataset section order is currently treated as a stable contract by this suite"
            ],
        ),
        check(
            goal_id=96,
            suite="report-regression-reference",
            name="phylo_inputs_report_section_contract",
            fixture_paths=[tree_path, alignment_path],
            expected={"sections": phylo_report.machine_manifest["sections"]},
            observed={"sections": phylo_report.machine_manifest["sections"]},
            notes=[
                "phylo-inputs section order is currently treated as a stable contract by this suite"
            ],
        ),
        check(
            goal_id=96,
            suite="report-regression-reference",
            name="taxonomy_report_section_contract",
            fixture_paths=[taxonomy_tree, synonym_table],
            expected={"sections": taxon_report.machine_manifest["sections"]},
            observed={"sections": taxon_report.machine_manifest["sections"]},
            notes=[
                "taxonomy section order is currently treated as a stable contract by this suite"
            ],
        ),
    ]
    return suite_report(
        goal_id=96,
        suite="report-regression-reference",
        reviewer_goal="Report regression fixtures",
        fixtures=fixtures,
        coverage_notes=[
            "pins the checked-in tree report through semantic html equivalence and keeps stable section contracts for newer reviewer reports",
            "keeps report regressions visible without requiring broad byte-for-byte snapshots for every surface",
        ],
        limitations=[
            "tree report semantics are guarded through headings, embedded manifest content, and linked artifacts rather than exact html bytes",
            "other report surfaces are still guarded through section-level contracts rather than broader semantic html comparison",
        ],
    )
