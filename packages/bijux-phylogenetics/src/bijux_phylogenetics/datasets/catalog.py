from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class PublicDatasetSurface:
    """One governed public dataset surface exposed through the demo runtime."""

    dataset_id: str
    label: str
    demo_command: str
    category: str
    summary: str


def list_public_dataset_surfaces() -> list[PublicDatasetSurface]:
    """Return the owned public dataset inventory for release reporting."""
    return [
        PublicDatasetSurface(
            dataset_id="primate_comparative",
            label="Primate comparative dataset",
            demo_command="primate-comparative",
            category="flagship",
            summary="Comparative trait workflow over one packaged primate tree and trait table.",
        ),
        PublicDatasetSurface(
            dataset_id="avian_reproductive_traits",
            label="Avian reproductive trait dataset",
            demo_command="avian-reproductive-traits",
            category="flagship",
            summary="Comparative trait workflow over one packaged avian life-history panel.",
        ),
        PublicDatasetSurface(
            dataset_id="central_european_seashore_flora",
            label="Central European seashore flora dataset",
            demo_command="central-european-seashore-flora",
            category="flagship",
            summary="Packaged plant comparative workflow with governed tree and trait outputs.",
        ),
        PublicDatasetSurface(
            dataset_id="influenza_a_ha_reference_panel",
            label="Influenza A HA reference panel",
            demo_command="influenza-a-ha-reference-panel",
            category="flagship",
            summary="Sequence-to-tree reference workflow over one packaged influenza HA panel.",
        ),
        PublicDatasetSurface(
            dataset_id="gnathostome_ortholog_protein_benchmark",
            label="Gnathostome ortholog protein benchmark",
            demo_command="gnathostome-ortholog-protein-benchmark",
            category="flagship",
            summary="Amino-acid sequence-to-tree benchmark with governed protein model selection and support review.",
        ),
        PublicDatasetSurface(
            dataset_id="pleistocene_bear_cytb_fragments",
            label="Pleistocene bear CYTB fragment panel",
            demo_command="pleistocene-bear-cytb-fragments",
            category="flagship",
            summary="Ancient-DNA-style degraded sequence workflow over one packaged bear fragment panel.",
        ),
        PublicDatasetSurface(
            dataset_id="rabies_cross_host_panel",
            label="Rabies cross-host panel",
            demo_command="rabies-cross-host-panel",
            category="flagship",
            summary="Host-transition review workflow over one packaged rabies host panel.",
        ),
        PublicDatasetSurface(
            dataset_id="rabies_geographic_transition_panel",
            label="Rabies geographic transition panel",
            demo_command="rabies-geographic-transition-panel",
            category="flagship",
            summary="Geographic transition review workflow over one packaged rabies geography panel.",
        ),
        PublicDatasetSurface(
            dataset_id="rabies_cross_host_geography_panel",
            label="Rabies cross-host geography panel",
            demo_command="rabies-cross-host-geography-panel",
            category="flagship",
            summary="Integrated flagship workflow connecting sequence-to-tree, host, geography, and comparative review.",
        ),
        PublicDatasetSurface(
            dataset_id="rabies_method_sensitivity_panel",
            label="Rabies method-sensitivity panel",
            demo_command="rabies-method-sensitivity-panel",
            category="flagship",
            summary="Governed method-sensitivity workflow spanning preprocessing and engine comparison over one rabies panel.",
        ),
        PublicDatasetSurface(
            dataset_id="catarrhine_mitogenome_five_locus_panel",
            label="Catarrhine mitogenome five-locus panel",
            demo_command="catarrhine-mitogenome-five-locus-panel",
            category="flagship",
            summary="Multi-locus concatenation and partitioned inference workflow over one packaged catarrhine panel.",
        ),
        PublicDatasetSurface(
            dataset_id="catarrhine_data_quality_stress_panel",
            label="Catarrhine data quality stress panel",
            demo_command="catarrhine-data-quality-stress-panel",
            category="stress",
            summary="Dirty-data audit and repair surface with explicit malformed sequence, tree, and trait inputs.",
        ),
        PublicDatasetSurface(
            dataset_id="continuous_mode_recovery_panel",
            label="Continuous-mode recovery panel",
            demo_command="continuous-mode-recovery-panel",
            category="reference",
            summary="Deterministic simulation-recovery panel for Brownian, OU, and early-burst comparative fits.",
        ),
        PublicDatasetSurface(
            dataset_id="discrete_mode_recovery_panel",
            label="Discrete-mode recovery panel",
            demo_command="discrete-mode-recovery-panel",
            category="reference",
            summary="Deterministic simulation-recovery panel for ER, SYM, and governed ARD comparative discrete Mk fits.",
        ),
        PublicDatasetSurface(
            dataset_id="known_answer_reference_panel",
            label="Known-answer reference panel",
            demo_command="known-answer-reference-panel",
            category="reference",
            summary="Owned simulation truth suite spanning topology, parameters, node states, and branch events.",
        ),
        PublicDatasetSurface(
            dataset_id="macroevolution_recovery_suite",
            label="Macroevolution recovery suite",
            demo_command="macroevolution-recovery-suite",
            category="reference",
            summary="Unified governed recovery suite over continuous, discrete, and truth-anchored macroevolution simulations.",
        ),
    ]


def list_flagship_dataset_surfaces() -> list[PublicDatasetSurface]:
    """Return the public flagship datasets advertised in the release surface."""
    return [
        surface
        for surface in list_public_dataset_surfaces()
        if surface.category == "flagship"
    ]
