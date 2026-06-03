from __future__ import annotations

from dataclasses import asdict
import json
from pathlib import Path
import shutil
import tempfile

from bijux_phylogenetics.io.fasta import (
    infer_alignment_alphabet,
    load_fasta_alignment,
)
from bijux_phylogenetics.io.newick import dumps_newick, write_newick
from bijux_phylogenetics.trees import (
    compute_clade_frequency_table,
    compute_consensus_tree,
    write_clade_frequency_table,
)

from .genetic_distance_matrix import (
    compute_pairwise_genetic_distance_matrix,
    write_genetic_distance_matrix,
)
from .models import (
    AmbiguityPolicy,
    DistanceGapPolicyDeltaRow,
    DistanceGapPolicySensitivityReport,
    DistanceMethodMaturityCheck,
    DistanceMethodMaturityGateReport,
    DistanceMethodReport,
    DistanceModel,
    DistanceModelComparisonReport,
    DistanceModelComparisonRow,
    DistanceReproducibilityBundleReport,
    GapHandlingMode,
)
from .quality import assess_distance_method_assumptions, inspect_distance_matrix_quality
from .shared import (
    _allowed_models_for_alphabet,
    _file_sha256,
    _pair_key,
    _require_supported_distance_tree_method,
    _unique_genetic_distance_pairs,
)
from .tree_inference import (
    bootstrap_distance_trees,
    build_distance_tree,
    compare_distance_tree_topologies,
    summarize_distance_bootstrap_support,
    write_distance_bootstrap_support,
)
from .validation import validate_distance_reference_examples


def compare_distance_models(
    path: Path,
    *,
    gap_handling: GapHandlingMode = "pairwise-deletion",
    ambiguity_policy: AmbiguityPolicy = "ignore",
) -> DistanceModelComparisonReport:
    """Compare all supported distance models for one alignment without overclaiming a winner."""
    records = load_fasta_alignment(path)
    inferred_alphabet = infer_alignment_alphabet(records)
    models = sorted(_allowed_models_for_alphabet(inferred_alphabet))
    rows: list[DistanceModelComparisonRow] = []
    for model in models:
        quality = inspect_distance_matrix_quality(
            path,
            model=model,
            gap_handling=gap_handling,
            ambiguity_policy=ambiguity_policy,
        )
        matrix = compute_pairwise_genetic_distance_matrix(
            path,
            model=model,
            gap_handling=gap_handling,
            ambiguity_policy=ambiguity_policy,
        )
        unique_pairs = _unique_genetic_distance_pairs(matrix)
        defined_distances = sorted(
            float(pair.distance) for pair in unique_pairs if pair.distance is not None
        )
        rows.append(
            DistanceModelComparisonRow(
                model=model,
                defined_pair_count=len(defined_distances),
                saturated_pair_count=len(quality.saturated_pairs),
                low_information_pair_count=len(quality.low_information_pairs),
                mean_distance=None
                if not defined_distances
                else round(sum(defined_distances) / len(defined_distances), 15),
                maximum_distance=None
                if not defined_distances
                else max(defined_distances),
                decision=quality.method_assessment.decision,
                reasons=quality.method_assessment.reasons,
            )
        )
    warnings: list[str] = []
    if any(row.saturated_pair_count > 0 for row in rows):
        warnings.append(
            "one or more supported models enter a saturation regime on this alignment"
        )
    if len(rows) < 2:
        warnings.append(
            "the inferred alphabet supports only one distance model, so cross-model sensitivity is limited"
        )
    return DistanceModelComparisonReport(
        alignment_path=path,
        inferred_alphabet=inferred_alphabet,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
        rows=rows,
        warnings=warnings,
    )


def compare_distance_gap_policies(
    path: Path,
    *,
    model: DistanceModel = "p-distance",
    ambiguity_policy: AmbiguityPolicy = "ignore",
) -> DistanceGapPolicySensitivityReport:
    """Summarize how pairwise versus complete deletion changes distance estimates."""
    pairwise = compute_pairwise_genetic_distance_matrix(
        path,
        model=model,
        gap_handling="pairwise-deletion",
        ambiguity_policy=ambiguity_policy,
    )
    complete = compute_pairwise_genetic_distance_matrix(
        path,
        model=model,
        gap_handling="complete-deletion",
        ambiguity_policy=ambiguity_policy,
    )
    pairwise_by_key = {
        _pair_key(pair.left_identifier, pair.right_identifier): pair
        for pair in _unique_genetic_distance_pairs(pairwise)
    }
    complete_by_key = {
        _pair_key(pair.left_identifier, pair.right_identifier): pair
        for pair in _unique_genetic_distance_pairs(complete)
    }
    rows: list[DistanceGapPolicyDeltaRow] = []
    for pair_key in sorted(set(pairwise_by_key) | set(complete_by_key)):
        pairwise_pair = pairwise_by_key.get(pair_key)
        complete_pair = complete_by_key.get(pair_key)
        if pairwise_pair is None or complete_pair is None:
            continue
        distance_delta = None
        if pairwise_pair.distance is not None and complete_pair.distance is not None:
            distance_delta = round(complete_pair.distance - pairwise_pair.distance, 15)
        if (
            pairwise_pair.distance != complete_pair.distance
            or pairwise_pair.comparable_sites != complete_pair.comparable_sites
        ):
            rows.append(
                DistanceGapPolicyDeltaRow(
                    left_identifier=pairwise_pair.left_identifier,
                    right_identifier=pairwise_pair.right_identifier,
                    pairwise_distance=pairwise_pair.distance,
                    complete_distance=complete_pair.distance,
                    pairwise_comparable_sites=pairwise_pair.comparable_sites,
                    complete_comparable_sites=complete_pair.comparable_sites,
                    distance_delta=distance_delta,
                    comparable_site_delta=complete_pair.comparable_sites
                    - pairwise_pair.comparable_sites,
                )
            )
    warnings: list[str] = []
    if rows:
        warnings.append(
            "distance estimates change when gap handling switches between pairwise and complete deletion"
        )
    if any(row.complete_comparable_sites == 0 for row in rows):
        warnings.append(
            "complete deletion removes all comparable sites for one or more taxon pairs"
        )
    return DistanceGapPolicySensitivityReport(
        alignment_path=path,
        model=model,
        ambiguity_policy=ambiguity_policy,
        changed_pair_count=len(rows),
        pair_count=len(pairwise_by_key),
        rows=rows,
        warnings=warnings,
    )


def assess_distance_method_maturity(
    path: Path,
    *,
    method: str = "neighbor-joining",
    model: DistanceModel = "p-distance",
    gap_handling: GapHandlingMode = "pairwise-deletion",
    ambiguity_policy: AmbiguityPolicy = "ignore",
    bootstrap_replicates: int = 25,
    bootstrap_seed: int = 1,
    validate_bundle: bool = False,
) -> DistanceMethodMaturityGateReport:
    """Evaluate whether distance-analysis surfaces are available and validated for one alignment."""
    method_policy = _require_supported_distance_tree_method(method)
    reference_validation = validate_distance_reference_examples()
    quality = inspect_distance_matrix_quality(
        path,
        model=model,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
    )
    assumptions = assess_distance_method_assumptions(
        path,
        model=model,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
    )
    bootstrap_summary = summarize_distance_bootstrap_support(
        bootstrap_distance_trees(
            path,
            method=method_policy.method,
            model=model,
            gap_handling=gap_handling,
            ambiguity_policy=ambiguity_policy,
            replicates=bootstrap_replicates,
            seed=bootstrap_seed,
        )[1]
    )
    model_comparison = compare_distance_models(
        path,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
    )
    gap_sensitivity = compare_distance_gap_policies(
        path,
        model=model,
        ambiguity_policy=ambiguity_policy,
    )

    bundle_satisfied = False
    bundle_details = "bundle validation was skipped"
    if validate_bundle:
        bundle_dir = (
            Path(tempfile.mkdtemp(prefix="bijux-distance-maturity-")) / "bundle"
        )
        bundle = write_distance_reproducibility_bundle(
            bundle_dir,
            alignment_path=path,
            method=method_policy.method,
            model=model,
            gap_handling=gap_handling,
            ambiguity_policy=ambiguity_policy,
            replicates=min(bootstrap_replicates, 10),
            seed=bootstrap_seed,
        )
        manifest = json.loads(
            (bundle.out_dir / "distance-analysis.manifest.json").read_text(
                encoding="utf-8"
            )
        )
        bundle_satisfied = (
            "output_checksums" in manifest and "input_checksums" in manifest
        )
        bundle_details = (
            "bundle manifest includes input and output checksums"
            if bundle_satisfied
            else "bundle manifest is missing checksum provenance"
        )

    checks = [
        DistanceMethodMaturityCheck(
            name="reference_validation",
            satisfied=reference_validation.all_passed,
            details="built-in p-distance, JC69, K2P, protein, NJ, and UPGMA reference cases all pass"
            if reference_validation.all_passed
            else "one or more built-in distance reference cases failed",
        ),
        DistanceMethodMaturityCheck(
            name="upgma_assumption_visibility",
            satisfied=bool(assumptions.upgma_assumptions)
            and assumptions.ultrametric_tolerance > 0.0,
            details="UPGMA assumptions and ultrametric violations are surfaced explicitly",
        ),
        DistanceMethodMaturityCheck(
            name="suitability_decision",
            satisfied=quality.method_assessment.decision
            in {"allowed", "risky", "blocked"},
            details=f"distance workflow received explicit decision '{quality.method_assessment.decision}'",
        ),
        DistanceMethodMaturityCheck(
            name="bootstrap_support_summary",
            satisfied=bootstrap_summary.clade_count >= 0,
            details="bootstrap support frequencies are summarized over site-resampled replicate trees",
        ),
        DistanceMethodMaturityCheck(
            name="model_comparison",
            satisfied=bool(model_comparison.rows),
            details="all supported distance models were compared on the same alignment",
        ),
        DistanceMethodMaturityCheck(
            name="gap_policy_sensitivity",
            satisfied=True,
            details=f"{gap_sensitivity.changed_pair_count} taxon pairs changed across pairwise versus complete deletion",
        ),
        DistanceMethodMaturityCheck(
            name="bundle_provenance",
            satisfied=bundle_satisfied if validate_bundle else True,
            details=bundle_details,
        ),
    ]
    if not all(check.satisfied for check in checks):
        decision = "blocked"
    elif quality.method_assessment.decision == "allowed":
        decision = "production_candidate"
    else:
        decision = "validated_with_limits"
    return DistanceMethodMaturityGateReport(
        alignment_path=path,
        method=method_policy.method,
        method_policy=method_policy,
        model=model,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
        decision=decision,
        checks=checks,
        warnings=quality.warnings
        + bootstrap_summary.warnings
        + model_comparison.warnings
        + gap_sensitivity.warnings,
    )


def build_distance_method_report(
    path: Path,
    *,
    method: str = "neighbor-joining",
    model: DistanceModel = "p-distance",
    gap_handling: GapHandlingMode = "pairwise-deletion",
    ambiguity_policy: AmbiguityPolicy = "ignore",
    bootstrap_replicates: int = 25,
    bootstrap_seed: int = 1,
) -> DistanceMethodReport:
    """Build a structured distance-method report that can back JSON or HTML renderers."""
    method_policy = _require_supported_distance_tree_method(method)
    matrix = compute_pairwise_genetic_distance_matrix(
        path,
        model=model,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
    )
    quality = inspect_distance_matrix_quality(
        path,
        model=model,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
    )
    assumptions = assess_distance_method_assumptions(
        path,
        model=model,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
    )
    reference_validation = validate_distance_reference_examples()
    tree, _ = build_distance_tree(
        path,
        method=method_policy.method,
        model=model,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
    )
    alternative_method = (
        "upgma" if method_policy.method == "neighbor-joining" else "neighbor-joining"
    )
    alternative_tree, _ = build_distance_tree(
        path,
        method=alternative_method,
        model=model,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
    )
    topology_comparison = compare_distance_tree_topologies(
        path,
        model=model,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
    )
    bootstrap_summary = summarize_distance_bootstrap_support(
        bootstrap_distance_trees(
            path,
            method=method_policy.method,
            model=model,
            gap_handling=gap_handling,
            ambiguity_policy=ambiguity_policy,
            replicates=bootstrap_replicates,
            seed=bootstrap_seed,
        )[1]
    )
    model_comparison = compare_distance_models(
        path,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
    )
    gap_policy_sensitivity = compare_distance_gap_policies(
        path,
        model=model,
        ambiguity_policy=ambiguity_policy,
    )
    maturity_gate = assess_distance_method_maturity(
        path,
        method=method_policy.method,
        model=model,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
        bootstrap_replicates=bootstrap_replicates,
        bootstrap_seed=bootstrap_seed,
        validate_bundle=False,
    )
    return DistanceMethodReport(
        alignment_path=path,
        method=method_policy.method,
        method_policy=method_policy,
        model=model,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
        matrix=matrix,
        quality=quality,
        assumptions=assumptions,
        reference_validation=reference_validation,
        built_tree_newick=dumps_newick(tree),
        alternative_tree_newick=dumps_newick(alternative_tree),
        topology_comparison=topology_comparison,
        bootstrap_summary=bootstrap_summary,
        model_comparison=model_comparison,
        gap_policy_sensitivity=gap_policy_sensitivity,
        maturity_gate=maturity_gate,
    )


def write_distance_reproducibility_bundle(
    out_dir: Path,
    *,
    alignment_path: Path,
    method: str,
    model: DistanceModel = "p-distance",
    gap_handling: GapHandlingMode = "pairwise-deletion",
    ambiguity_policy: AmbiguityPolicy = "ignore",
    replicates: int = 100,
    seed: int = 1,
) -> DistanceReproducibilityBundleReport:
    """Write a reproducibility bundle for one distance analysis."""
    method_policy = _require_supported_distance_tree_method(method)
    out_dir.mkdir(parents=True, exist_ok=True)
    bundled_alignment_path = out_dir / "input-alignment.fasta"
    shutil.copy2(alignment_path, bundled_alignment_path)
    matrix = compute_pairwise_genetic_distance_matrix(
        alignment_path,
        model=model,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
    )
    quality = inspect_distance_matrix_quality(
        alignment_path,
        model=model,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
    )
    tree, _ = build_distance_tree(
        alignment_path,
        method=method_policy.method,
        model=model,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
    )
    bootstrap_trees, bootstrap = bootstrap_distance_trees(
        alignment_path,
        method=method_policy.method,
        model=model,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
        replicates=replicates,
        seed=seed,
    )
    matrix_path = write_genetic_distance_matrix(out_dir / "distance-matrix.tsv", matrix)
    tree_path = write_newick(out_dir / "distance-tree.nwk", tree)
    support_path = write_distance_bootstrap_support(
        out_dir / "bootstrap-support.tsv", bootstrap
    )
    from bijux_phylogenetics.simulation import write_tree_set

    tree_set_path = write_tree_set(
        out_dir / "bootstrap-replicates.trees", bootstrap_trees
    )
    consensus_tree, _ = compute_consensus_tree(tree_set_path)
    consensus_path = write_newick(out_dir / "bootstrap-consensus.nwk", consensus_tree)
    clade_frequency_path = write_clade_frequency_table(
        out_dir / "bootstrap-clades.tsv",
        compute_clade_frequency_table(tree_set_path),
    )
    summary = build_distance_method_report(
        alignment_path,
        method=method_policy.method,
        model=model,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
        bootstrap_replicates=min(replicates, 25),
        bootstrap_seed=seed,
    )
    summary_path = out_dir / "distance-summary.json"
    summary_path.write_text(
        json.dumps(asdict(summary), default=str, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    manifest_path = out_dir / "distance-analysis.manifest.json"
    manifest_path.write_text(
        json.dumps(
            {
                "alignment_path": str(alignment_path),
                "method": method_policy.method,
                "model": model,
                "gap_handling": gap_handling,
                "ambiguity_policy": ambiguity_policy,
                "replicates": replicates,
                "seed": seed,
                "quality_decision": quality.method_assessment.decision,
                "quality_reasons": quality.method_assessment.reasons,
                "files": [
                    "input-alignment.fasta",
                    "distance-matrix.tsv",
                    "distance-tree.nwk",
                    "bootstrap-support.tsv",
                    "bootstrap-replicates.trees",
                    "bootstrap-consensus.nwk",
                    "bootstrap-clades.tsv",
                    "distance-summary.json",
                ],
                "input_checksums": {
                    str(alignment_path): _file_sha256(alignment_path),
                },
                "output_checksums": {
                    "input-alignment.fasta": _file_sha256(bundled_alignment_path),
                    "distance-matrix.tsv": _file_sha256(matrix_path),
                    "distance-tree.nwk": _file_sha256(tree_path),
                    "bootstrap-support.tsv": _file_sha256(support_path),
                    "bootstrap-replicates.trees": _file_sha256(tree_set_path),
                    "bootstrap-consensus.nwk": _file_sha256(consensus_path),
                    "bootstrap-clades.tsv": _file_sha256(clade_frequency_path),
                    "distance-summary.json": _file_sha256(summary_path),
                },
                "reference_validation_passed": validate_distance_reference_examples().all_passed,
                "caveats": [
                    "distance methods reduce site-by-site evidence to pairwise summaries before tree building",
                    "bootstrap support measures repeatability under site resampling, not absolute correctness",
                    *quality.warnings,
                ],
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    files = [
        bundled_alignment_path,
        matrix_path,
        tree_path,
        support_path,
        tree_set_path,
        consensus_path,
        clade_frequency_path,
        summary_path,
        manifest_path,
    ]
    return DistanceReproducibilityBundleReport(
        out_dir=out_dir,
        alignment_path=alignment_path,
        method=method_policy.method,
        model=model,
        gap_handling=gap_handling,
        ambiguity_policy=ambiguity_policy,
        replicates=replicates,
        files=files,
    )
