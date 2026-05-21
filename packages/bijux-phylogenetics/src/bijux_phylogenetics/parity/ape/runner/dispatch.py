from __future__ import annotations

from ..registry import ApeParityCase
from .ancestral_payloads import (
    _build_bijux_continuous_ancestral_rows,
    _build_bijux_discrete_ancestral_rows,
)
from .comparative_payloads import (
    _build_bijux_branching_time_rows,
    _build_bijux_brownian_covariance_rows,
    _build_bijux_diversification_gamma_rows,
    _build_bijux_independent_contrast_rows,
    _build_bijux_tree_node_depth_rows,
    _build_bijux_tree_simulation_envelope_rows,
    _build_bijux_tree_ultrametric_rows,
)
from .sequence_payloads import (
    _build_bijux_base_frequency_summary,
    _build_bijux_distance_rows,
    _build_bijux_dnabin_rows,
    _build_bijux_neighbor_joining_structure,
    _build_bijux_segregating_site_rows,
    _build_bijux_translation_rows,
)
from .tree_payloads import (
    _build_bijux_consensus_rows,
    _build_bijux_drop_tip_structure,
    _build_bijux_extract_clade_structure,
    _build_bijux_keep_tip_structure,
    _build_bijux_monophyly_summary,
    _build_bijux_mrca_summary,
    _build_bijux_prop_clades_rows,
    _build_bijux_root_outgroup_structure,
    _build_bijux_topology_distance_rows,
    _build_bijux_tree_set_structure,
    _build_bijux_tree_structure,
    _build_bijux_tree_tip_distance_rows,
    _build_bijux_unroot_structure,
)


def _build_bijux_case_payload(
    case: ApeParityCase,
) -> tuple[dict[str, object], list[dict[str, object]] | None, str | None]:
    if case.operation in {"read-tree-structure", "write-tree-structure"}:
        summary, rows, normalized_text = _build_bijux_tree_structure(case.input_fixture)
        return summary, rows, normalized_text
    if case.operation == "root-tree-outgroup":
        summary, rows, normalized_text = _build_bijux_root_outgroup_structure(
            case.input_fixture,
            outgroup_taxa=case.outgroup_taxa,
        )
        return summary, rows, normalized_text
    if case.operation == "unroot-tree":
        summary, rows, normalized_text = _build_bijux_unroot_structure(
            case.input_fixture,
        )
        return summary, rows, normalized_text
    if case.operation == "drop-tree-taxa":
        summary, rows, normalized_text = _build_bijux_drop_tip_structure(
            case.input_fixture,
            excluded_taxa=case.excluded_taxa,
        )
        return summary, rows, normalized_text
    if case.operation == "keep-tree-taxa":
        summary, rows, normalized_text = _build_bijux_keep_tip_structure(
            case.input_fixture,
            requested_taxa=case.requested_taxa,
        )
        return summary, rows, normalized_text
    if case.operation == "extract-tree-clade":
        if case.node_id is None:
            raise ValueError(
                f"ape parity case '{case.case_id}' is missing an extraction node id"
            )
        summary, rows, normalized_text = _build_bijux_extract_clade_structure(
            case.input_fixture,
            node_id=case.node_id,
        )
        return summary, rows, normalized_text
    if case.operation == "get-tree-mrca":
        if not case.mrca_taxa:
            raise ValueError(f"ape parity case '{case.case_id}' is missing MRCA taxa")
        summary = _build_bijux_mrca_summary(
            case.input_fixture,
            mrca_taxa=case.mrca_taxa,
        )
        return summary, None, None
    if case.operation == "assess-tree-monophyly":
        if case.monophyly_reroot is None:
            raise ValueError(
                f"ape parity case '{case.case_id}' is missing a monophyly reroot policy"
            )
        summary = _build_bijux_monophyly_summary(
            case.input_fixture,
            requested_taxa=case.requested_taxa,
            reroot=case.monophyly_reroot,
        )
        return summary, None, None
    if case.operation in {"read-tree-set-structure", "write-tree-set-structure"}:
        summary, rows, normalized_text = _build_bijux_tree_set_structure(
            case.input_fixture
        )
        return summary, rows, normalized_text
    if case.operation == "tree-consensus":
        if case.consensus_method is None:
            raise ValueError(
                f"ape parity case '{case.case_id}' is missing a consensus method"
            )
        summary, rows, normalized_text = _build_bijux_consensus_rows(
            case.input_fixture,
            consensus_method=case.consensus_method,
        )
        return summary, rows, normalized_text
    if case.operation == "tree-clade-support":
        if case.reference_tree_path is None:
            raise ValueError(
                f"ape parity case '{case.case_id}' is missing a reference tree path"
            )
        summary, rows = _build_bijux_prop_clades_rows(
            case.reference_tree_path,
            case.input_fixture,
        )
        return summary, rows, None
    if case.operation == "tree-tip-distance":
        summary, rows = _build_bijux_tree_tip_distance_rows(case.input_fixture)
        return summary, rows, None
    if case.operation == "distance-matrix-neighbor-joining":
        return _build_bijux_neighbor_joining_structure(case.input_fixture)
    if case.operation == "tree-topology-distance":
        if case.rf_mode is None:
            raise ValueError(
                f"ape parity case '{case.case_id}' is missing a topology rf mode"
            )
        summary, rows = _build_bijux_topology_distance_rows(
            case.input_fixture,
            rf_mode=case.rf_mode,
        )
        return summary, rows, None
    if case.operation == "tree-brownian-covariance":
        summary, rows = _build_bijux_brownian_covariance_rows(case.input_fixture)
        return summary, rows, None
    if case.operation == "tree-continuous-ancestral-states":
        if (
            case.trait_table_path is None
            or case.trait_name is None
            or case.trait_taxon_column is None
        ):
            raise ValueError(
                f"ape parity case '{case.case_id}' is missing a trait table path, trait name, or taxon column"
            )
        summary, rows = _build_bijux_continuous_ancestral_rows(
            case.input_fixture,
            trait_table_path=case.trait_table_path,
            trait_name=case.trait_name,
            trait_taxon_column=case.trait_taxon_column,
        )
        return summary, rows, None
    if case.operation == "tree-discrete-ancestral-states":
        if (
            case.trait_table_path is None
            or case.trait_name is None
            or case.trait_taxon_column is None
        ):
            raise ValueError(
                f"ape parity case '{case.case_id}' is missing a trait table path, trait name, or taxon column"
            )
        summary, rows = _build_bijux_discrete_ancestral_rows(
            case.input_fixture,
            trait_table_path=case.trait_table_path,
            trait_name=case.trait_name,
            trait_taxon_column=case.trait_taxon_column,
            ancestral_model=case.ancestral_model or "equal-rates",
        )
        return summary, rows, None
    if case.operation == "tree-independent-contrasts":
        if case.trait_table_path is None or case.trait_name is None:
            raise ValueError(
                f"ape parity case '{case.case_id}' is missing a trait table path or trait name"
            )
        summary, rows = _build_bijux_independent_contrast_rows(
            case.input_fixture,
            trait_table_path=case.trait_table_path,
            trait_name=case.trait_name,
        )
        return summary, rows, None
    if case.operation == "tree-node-depth":
        summary, rows = _build_bijux_tree_node_depth_rows(case.input_fixture)
        return summary, rows, None
    if case.operation == "tree-branching-times":
        summary, rows = _build_bijux_branching_time_rows(case.input_fixture)
        return summary, rows, None
    if case.operation == "tree-diversification-gamma-statistic":
        summary, rows = _build_bijux_diversification_gamma_rows(case.input_fixture)
        return summary, rows, None
    if case.operation == "tree-simulation-envelope":
        summary, rows = _build_bijux_tree_simulation_envelope_rows(case.fixture_id)
        return summary, rows, None
    if case.operation == "tree-ultrametricity":
        if case.ultrametric_option is None:
            raise ValueError(
                f"ape parity case '{case.case_id}' is missing an ultrametric option"
            )
        summary, rows = _build_bijux_tree_ultrametric_rows(
            case.input_fixture,
            tolerance=case.tolerance,
            option=case.ultrametric_option,
        )
        return summary, rows, None
    if case.operation == "dna-dnabin-structure":
        summary, rows = _build_bijux_dnabin_rows(case.input_fixture)
        return summary, rows, None
    if case.operation == "dna-base-frequency":
        summary, rows = _build_bijux_base_frequency_summary(case.input_fixture)
        return summary, rows, None
    if case.operation == "dna-segregating-sites":
        summary, rows = _build_bijux_segregating_site_rows(case.input_fixture)
        return summary, rows, None
    if case.operation == "dna-distance":
        if case.pairwise_deletion is None:
            raise ValueError(
                f"ape parity case '{case.case_id}' is missing pairwise deletion policy"
            )
        if case.distance_model is None:
            raise ValueError(
                f"ape parity case '{case.case_id}' is missing a distance model"
            )
        summary, rows = _build_bijux_distance_rows(
            case.input_fixture,
            pairwise_deletion=case.pairwise_deletion,
            distance_model=case.distance_model,
        )
        return summary, rows, None
    if case.operation == "dna-translation":
        if case.genetic_code_id is None:
            raise ValueError(
                f"ape parity case '{case.case_id}' is missing a genetic code id"
            )
        summary, rows = _build_bijux_translation_rows(
            case.input_fixture,
            genetic_code_id=case.genetic_code_id,
        )
        return summary, rows, None
    raise ValueError(f"unsupported ape parity operation '{case.operation}'")
