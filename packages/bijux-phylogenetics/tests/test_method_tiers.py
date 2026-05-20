from __future__ import annotations

import pytest

from bijux_phylogenetics.evidence.provenance.method_tiers import (
    MethodTierAssessment,
    bayesian_report_method_tier,
    comparative_report_method_tier,
    fasta_to_tree_method_tier,
    method_tier_metrics,
    method_tier_warnings,
    phylogenetic_logistic_method_tier,
    release_method_tier_inventory,
    tree_report_method_tier,
)


def test_supported_method_tiers_require_real_validation_basis() -> None:
    with pytest.raises(ValueError, match="supported method tiers require"):
        MethodTierAssessment(
            tier="supported",
            surface="unsupported-test-surface",
            summary="bad supported tier",
            validation_basis=("governed-fixture:only",),
            inference_mode="inference",
        )


def test_experimental_method_tiers_require_clear_warning() -> None:
    with pytest.raises(ValueError, match="experimental method tiers require"):
        MethodTierAssessment(
            tier="experimental",
            surface="experimental-test-surface",
            summary="bad experimental tier",
            validation_basis=(),
            inference_mode="inference",
        )


def test_parser_only_method_tiers_reject_inference_claims() -> None:
    with pytest.raises(ValueError, match="parser-only method tiers must declare"):
        MethodTierAssessment(
            tier="parser-only",
            surface="parser-only-test-surface",
            summary="bad parser-only tier",
            validation_basis=("parser-contract:test",),
            inference_mode="inference",
        )


def test_surface_specific_method_tiers_expose_expected_runtime_contracts() -> None:
    supported = fasta_to_tree_method_tier()
    experimental = phylogenetic_logistic_method_tier(
        "phylogenetic-working-correlation-gee"
    )
    advisory = tree_report_method_tier()
    parser_only = bayesian_report_method_tier("bayesian-posterior")
    comparative = comparative_report_method_tier()

    assert supported.tier == "supported"
    assert supported.inference_mode == "inference"
    assert any(
        basis.startswith("real-engine-validation:")
        for basis in supported.validation_basis
    )

    assert experimental.tier == "experimental"
    assert experimental.approximation == "phylogenetic-working-correlation-gee"
    assert experimental.excluded_reference_surfaces == ("ape::compar.gee",)
    assert method_tier_warnings(experimental)

    assert advisory.tier == "advisory"
    assert advisory.inference_mode == "review-only"

    assert parser_only.tier == "parser-only"
    assert parser_only.inference_mode == "parser-only"
    assert "does not claim" in parser_only.summary

    assert comparative.tier == "supported"
    assert any(
        basis.startswith("reference-parity:") for basis in comparative.validation_basis
    )


def test_method_tier_metrics_emit_cli_ready_fields() -> None:
    metrics = method_tier_metrics(phylogenetic_logistic_method_tier("gee"))

    assert metrics == {
        "method_tier": "experimental",
        "method_surface": "phylogenetic-logistic",
        "method_inference_mode": "inference",
        "method_validation_basis": [],
        "method_approximation": "gee",
        "method_excluded_reference_surfaces": ["ape::compar.gee"],
    }


def test_phylogenetic_logistic_method_tier_warnings_include_explicit_non_claim() -> (
    None
):
    warnings = method_tier_warnings(
        phylogenetic_logistic_method_tier("phylogenetic-working-correlation-gee")
    )

    assert warnings[0].startswith("experimental method tier:")
    assert any("ape::compar.gee parity" in warning for warning in warnings)


def test_release_method_tier_inventory_lists_governed_release_surfaces() -> None:
    inventory = release_method_tier_inventory()

    assert any(item.surface == "fasta-to-tree" for item in inventory)
    assert any(item.surface == "phylogenetic-logistic" for item in inventory)
    assert any(item.surface == "bayesian-posterior-report" for item in inventory)
