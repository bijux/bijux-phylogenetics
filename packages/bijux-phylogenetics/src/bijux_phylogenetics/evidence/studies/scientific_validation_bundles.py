from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path

from ..bundle_artifacts import build_bundle_governed_artifacts
from ..study_contracts import DEFAULT_OWNER_PACKAGE, load_study_contract

SCIENTIFIC_VALIDATION_BUNDLE_GENERATOR_PATH = (
    "packages/bijux-phylogenetics/src/bijux_phylogenetics/evidence/studies/"
    "scientific_validation_bundles.py"
)


@dataclass(frozen=True, slots=True)
class ScientificValidationBundleDefinition:
    evidence_id: str
    evidence_title: str
    summary: str
    claim_id: str
    analytical_surfaces: tuple[str, ...]
    claim_tags: tuple[str, ...]
    source_fragments: tuple[str, ...]
    limitations: tuple[str, ...]
    next_step: str


@dataclass(frozen=True, slots=True)
class ScientificValidationStudyDefinition:
    study_id: str
    verdict_summary: str
    comparison_mode: str
    current_repository_posture: str
    governed_code_paths: tuple[str, ...]
    source_basis_locators: tuple[str, ...]
    bundle_definitions: tuple[ScientificValidationBundleDefinition, ...]


SCIENTIFIC_VALIDATION_STUDIES: dict[str, ScientificValidationStudyDefinition] = {
    "primate-continuous-and-categorical-regression": ScientificValidationStudyDefinition(
        study_id="primate-continuous-and-categorical-regression",
        verdict_summary=(
            "This bundle is checked in and reviewable, but lecture-specific PCM3"
            " parity for this slice is not yet claimed by the repository."
        ),
        comparison_mode="direct_parity",
        current_repository_posture=(
            "comparative_regression_surface_available_but_lecture_parity_not_yet_governed"
        ),
        governed_code_paths=(
            SCIENTIFIC_VALIDATION_BUNDLE_GENERATOR_PATH,
            "packages/bijux-phylogenetics/src/bijux_phylogenetics/comparative",
        ),
        source_basis_locators=(
            "evidence-book/studies/primate-continuous-and-categorical-regression/provenance/lund-course-sources.json",
            "evidence-book/studies/primate-continuous-and-categorical-regression/datasets/registry.json",
            "evidence-book/studies/primate-continuous-and-categorical-regression/reference/README.md",
            "packages/bijux-phylogenetics/tests/test_pgls_interaction_coefficients.py",
            "packages/bijux-phylogenetics/tests/test_regression_model_selection.py",
        ),
        bundle_definitions=(
            ScientificValidationBundleDefinition(
                evidence_id="evidence-001",
                evidence_title="Primate source intake and workspace alignment bundle",
                summary=(
                    "Governed study slice for PCM3 source intake, workspace"
                    " reconstruction, and shared primate dataset alignment."
                ),
                claim_id="pcm3-source-intake-and-workspace-alignment",
                analytical_surfaces=("workflow-contracts", "dataset-alignment"),
                claim_tags=(
                    "scientific-validation",
                    "primate",
                    "workflow",
                    "dataset-alignment",
                ),
                source_fragments=(
                    "workspace-reload-contract",
                    "shared-dataset-alignment",
                ),
                limitations=(
                    "This bundle records the checked-in study contract and does not yet claim lecture-specific object-for-object parity.",
                ),
                next_step=(
                    "Add a lecture-backed parity ledger for shared workspace objects"
                    " and exclusion accounting."
                ),
            ),
            ScientificValidationBundleDefinition(
                evidence_id="evidence-002",
                evidence_title="Primate baseline regression structure bundle",
                summary=(
                    "Governed study slice for the PCM3 baseline GLS fit with one"
                    " continuous predictor and one categorical predictor."
                ),
                claim_id="pcm3-baseline-gls-with-categorical-predictor",
                analytical_surfaces=("baseline-regression", "categorical-predictors"),
                claim_tags=(
                    "scientific-validation",
                    "primate",
                    "baseline-regression",
                    "categorical-predictor",
                ),
                source_fragments=("baseline-gls-fit", "categorical-encoding-contract"),
                limitations=(
                    "The repository exposes nearby regression surfaces, but this bundle does not yet claim lecture-specific coefficient parity.",
                ),
                next_step=(
                    "Check in coefficient and model-fit ledgers for the baseline"
                    " categorical regression surface."
                ),
            ),
            ScientificValidationBundleDefinition(
                evidence_id="evidence-003",
                evidence_title="Primate transformed response and model test bundle",
                summary=(
                    "Governed study slice for the PCM3 log-response repair,"
                    " residual review, and marginal-versus-sequential testing logic."
                ),
                claim_id="pcm3-transformed-response-and-model-tests",
                analytical_surfaces=("response-transformation", "model-comparison"),
                claim_tags=(
                    "scientific-validation",
                    "primate",
                    "response-transformation",
                    "model-tests",
                ),
                source_fragments=(
                    "log-response-repair",
                    "marginal-versus-sequential-tests",
                ),
                limitations=(
                    "This bundle documents the test structure and response transformation without claiming checked-in lecture-level parity outputs.",
                ),
                next_step=(
                    "Add governed ledgers for transformed-response coefficients,"
                    " residual summaries, and hypothesis-test rows."
                ),
            ),
            ScientificValidationBundleDefinition(
                evidence_id="evidence-004",
                evidence_title="Primate intercept-free contrast bundle",
                summary=(
                    "Governed study slice for the PCM3 intercept-free"
                    " parameterization used to report group-specific estimates."
                ),
                claim_id="pcm3-intercept-free-group-parameterization",
                analytical_surfaces=("contrast-parameterization",),
                claim_tags=(
                    "scientific-validation",
                    "primate",
                    "contrasts",
                    "parameterization",
                ),
                source_fragments=("intercept-free-parameterization",),
                limitations=(
                    "The repository documents the parameterization boundary clearly, but it does not yet claim lecture-specific group-estimate parity.",
                ),
                next_step=(
                    "Publish governed contrast ledgers for group-specific"
                    " estimates under the intercept-free parameterization."
                ),
            ),
            ScientificValidationBundleDefinition(
                evidence_id="evidence-005",
                evidence_title="Primate interaction coefficient interpretation bundle",
                summary=(
                    "Governed study slice for PCM3 continuous-by-categorical"
                    " interaction coefficients and their interpretation."
                ),
                claim_id="pcm3-interaction-coefficient-interpretation",
                analytical_surfaces=("interaction-coefficients",),
                claim_tags=(
                    "scientific-validation",
                    "primate",
                    "interaction-terms",
                    "coefficient-interpretation",
                ),
                source_fragments=("continuous-by-categorical-interaction",),
                limitations=(
                    "Nearby runtime summaries exist for interactions, but this bundle does not yet claim one lecture-specific parity table.",
                ),
                next_step=(
                    "Check in one governed interaction-coefficient ledger tied"
                    " to the lecture term ordering and interpretation rules."
                ),
            ),
            ScientificValidationBundleDefinition(
                evidence_id="evidence-006",
                evidence_title="Primate scaling and intercept semantics bundle",
                summary=(
                    "Governed study slice for PCM3 predictor scaling and the"
                    " intercept-semantics changes that follow standardization."
                ),
                claim_id="pcm3-scaling-and-intercept-semantics",
                analytical_surfaces=("predictor-scaling", "intercept-semantics"),
                claim_tags=(
                    "scientific-validation",
                    "primate",
                    "scaling",
                    "intercept-semantics",
                ),
                source_fragments=(
                    "predictor-scaling",
                    "standardized-intercept-semantics",
                ),
                limitations=(
                    "This bundle keeps the standardization story explicit without yet claiming checked-in parity for the lecture’s scaled outputs.",
                ),
                next_step=(
                    "Add governed before-and-after ledgers for scaled predictors,"
                    " standardized coefficients, and intercept reinterpretation."
                ),
            ),
            ScientificValidationBundleDefinition(
                evidence_id="evidence-007",
                evidence_title="Primate Martins covariance regression bundle",
                summary=(
                    "Governed study slice for the PCM3 phylogenetic regression"
                    " that uses categorical predictors with Martins covariance."
                ),
                claim_id="pcm3-martins-covariance-phylogenetic-regression",
                analytical_surfaces=("phylogenetic-regression", "martins-covariance"),
                claim_tags=(
                    "scientific-validation",
                    "primate",
                    "phylogenetic-regression",
                    "martins-covariance",
                ),
                source_fragments=(
                    "martins-covariance-fit",
                    "categorical-phylogenetic-regression",
                ),
                limitations=(
                    "The repository documents the covariance and model structure clearly, but it does not yet claim lecture-specific Martins-covariance parity.",
                ),
                next_step=(
                    "Publish governed coefficient, covariance, and fit ledgers"
                    " for the Martins-covariance regression surface."
                ),
            ),
            ScientificValidationBundleDefinition(
                evidence_id="evidence-008",
                evidence_title="Primate post-correction diagnostics bundle",
                summary=(
                    "Governed study slice for PCM3 post-correction diagnostics"
                    " and hypothesis testing after the phylogenetic regression fit."
                ),
                claim_id="pcm3-post-correction-diagnostics-and-tests",
                analytical_surfaces=("diagnostics", "post-correction-tests"),
                claim_tags=(
                    "scientific-validation",
                    "primate",
                    "diagnostics",
                    "hypothesis-testing",
                ),
                source_fragments=(
                    "post-correction-diagnostics",
                    "post-correction-hypothesis-tests",
                ),
                limitations=(
                    "This bundle keeps the diagnostic and testing boundary explicit without claiming a finished lecture-parity package for those outputs.",
                ),
                next_step=(
                    "Add governed residual, fit, and test ledgers for the"
                    " post-correction regression review surface."
                ),
            ),
        ),
    ),
    "cichlid-bayesian-linear-and-phylogenetic-mixed-models": ScientificValidationStudyDefinition(
        study_id="cichlid-bayesian-linear-and-phylogenetic-mixed-models",
        verdict_summary=(
            "This bundle is checked in and reviewable, but the repository does"
            " not yet claim governed parity for the PCM4 posterior workflow."
        ),
        comparison_mode="direct_parity",
        current_repository_posture=(
            "bayesian_mixed_model_parity_open_with_pgls_comparator_surface_only"
        ),
        governed_code_paths=(
            SCIENTIFIC_VALIDATION_BUNDLE_GENERATOR_PATH,
            "packages/bijux-phylogenetics/src/bijux_phylogenetics/comparative",
            "packages/bijux-phylogenetics/src/bijux_phylogenetics/phylo",
        ),
        source_basis_locators=(
            "evidence-book/studies/cichlid-bayesian-linear-and-phylogenetic-mixed-models/provenance/lund-course-sources.json",
            "evidence-book/studies/cichlid-bayesian-linear-and-phylogenetic-mixed-models/datasets/registry.json",
            "evidence-book/studies/cichlid-bayesian-linear-and-phylogenetic-mixed-models/reference/README.md",
        ),
        bundle_definitions=(
            ScientificValidationBundleDefinition(
                evidence_id="evidence-001",
                evidence_title="Cichlid alignment and species summarization bundle",
                summary=(
                    "Governed study slice for PCM4 cichlid data-tree alignment"
                    " and species-level summarization before model fitting."
                ),
                claim_id="pcm4-cichlid-alignment-and-species-summarization",
                analytical_surfaces=("tree-alignment", "species-summarization"),
                claim_tags=(
                    "scientific-validation",
                    "cichlid",
                    "tree-alignment",
                    "species-summarization",
                ),
                source_fragments=(
                    "cichlid-data-tree-alignment",
                    "species-level-summarization",
                ),
                limitations=(
                    "This bundle documents the shared species-level contract without claiming lecture-specific parity for the derived workspace.",
                ),
                next_step=(
                    "Check in governed species-alignment and summarization ledgers"
                    " tied to the lecture preprocessing path."
                ),
            ),
            ScientificValidationBundleDefinition(
                evidence_id="evidence-002",
                evidence_title="Cichlid Gaussian baseline regression bundle",
                summary=(
                    "Governed study slice for the PCM4 non-phylogenetic Gaussian"
                    " baseline regression used before mixed-model sampling."
                ),
                claim_id="pcm4-gaussian-baseline-regression",
                analytical_surfaces=("baseline-regression", "gaussian-model"),
                claim_tags=(
                    "scientific-validation",
                    "cichlid",
                    "gaussian-model",
                    "baseline-regression",
                ),
                source_fragments=("gaussian-baseline-fit",),
                limitations=(
                    "The baseline Gaussian fit is documented as its own slice, but this bundle does not yet claim checked-in lecture-parity coefficients.",
                ),
                next_step=(
                    "Add governed baseline-fit ledgers so the shared response and"
                    " predictor contract is closed before posterior review."
                ),
            ),
            ScientificValidationBundleDefinition(
                evidence_id="evidence-003",
                evidence_title="Cichlid Gaussian posterior diagnostics bundle",
                summary=(
                    "Governed study slice for the PCM4 Bayesian Gaussian"
                    " regression chain diagnostics and posterior summaries."
                ),
                claim_id="pcm4-gaussian-posterior-diagnostics",
                analytical_surfaces=("posterior-diagnostics", "gaussian-posterior"),
                claim_tags=(
                    "scientific-validation",
                    "cichlid",
                    "posterior-diagnostics",
                    "gaussian-posterior",
                ),
                source_fragments=(
                    "gaussian-chain-diagnostics",
                    "gaussian-posterior-summaries",
                ),
                limitations=(
                    "The repository does not yet own the lecture posterior workflow, so this bundle remains a governed boundary record rather than a parity claim.",
                ),
                next_step=(
                    "Check in governed posterior summary ledgers and convergence"
                    " diagnostics for the Gaussian mixed-model chain."
                ),
            ),
            ScientificValidationBundleDefinition(
                evidence_id="evidence-004",
                evidence_title="Cichlid phylogenetic inverse setup bundle",
                summary=(
                    "Governed study slice for the PCM4 tree inverse construction"
                    " and phylogenetic random-effect setup."
                ),
                claim_id="pcm4-tree-inverse-and-random-effect-setup",
                analytical_surfaces=("phylogenetic-inverse", "random-effect-setup"),
                claim_tags=(
                    "scientific-validation",
                    "cichlid",
                    "phylogenetic-inverse",
                    "random-effects",
                ),
                source_fragments=(
                    "tree-inverse-construction",
                    "phylogenetic-random-effect-setup",
                ),
                limitations=(
                    "This bundle keeps the inverse-covariance setup reviewable without claiming that the downstream posterior fit is already governed in-repository.",
                ),
                next_step=(
                    "Add governed matrix and random-effect setup ledgers linked to"
                    " the lecture’s phylogenetic inverse workflow."
                ),
            ),
            ScientificValidationBundleDefinition(
                evidence_id="evidence-005",
                evidence_title="Cichlid phylogenetic fixed-effect posterior bundle",
                summary=(
                    "Governed study slice for the PCM4 Bayesian phylogenetic"
                    " mixed-model fixed-effect posterior summaries."
                ),
                claim_id="pcm4-phylogenetic-fixed-effect-posteriors",
                analytical_surfaces=(
                    "phylogenetic-mixed-model",
                    "fixed-effect-posteriors",
                ),
                claim_tags=(
                    "scientific-validation",
                    "cichlid",
                    "phylogenetic-mixed-model",
                    "fixed-effects",
                ),
                source_fragments=("fixed-effect-posteriors",),
                limitations=(
                    "The lecture posterior fixed-effect summaries remain external to the current governed runtime, so this bundle stays boundary-first.",
                ),
                next_step=(
                    "Check in governed fixed-effect posterior ledgers with explicit"
                    " uncertainty summaries and parameter naming."
                ),
            ),
            ScientificValidationBundleDefinition(
                evidence_id="evidence-006",
                evidence_title="Cichlid variance-component posterior bundle",
                summary=(
                    "Governed study slice for the PCM4 phylogenetic"
                    " variance-component posterior summaries."
                ),
                claim_id="pcm4-phylogenetic-variance-component-posteriors",
                analytical_surfaces=("variance-components", "phylogenetic-posterior"),
                claim_tags=(
                    "scientific-validation",
                    "cichlid",
                    "variance-components",
                    "phylogenetic-posterior",
                ),
                source_fragments=("variance-component-posteriors",),
                limitations=(
                    "This bundle keeps the variance-component question explicit without claiming repository-owned parity for the posterior sampler.",
                ),
                next_step=(
                    "Add governed variance-component ledgers with posterior means,"
                    " intervals, and phylogenetic share calculations."
                ),
            ),
            ScientificValidationBundleDefinition(
                evidence_id="evidence-007",
                evidence_title="Cichlid Brownian comparator bundle",
                summary=(
                    "Governed study slice for the PCM4 comparison back to a"
                    " Brownian PGLS baseline after the mixed-model workflow."
                ),
                claim_id="pcm4-brownian-pgls-comparator",
                analytical_surfaces=("brownian-comparator", "comparative-regression"),
                claim_tags=(
                    "scientific-validation",
                    "cichlid",
                    "brownian-pgls",
                    "comparator",
                ),
                source_fragments=("brownian-pgls-comparison",),
                limitations=(
                    "The repository can support comparator surfaces nearby, but this bundle does not yet claim a finished lecture-parity bridge between the posterior model and the Brownian comparator.",
                ),
                next_step=(
                    "Check in one governed comparator ledger that ties the shared"
                    " response and predictor contract back to Brownian PGLS."
                ),
            ),
        ),
    ),
    "cichlid-phylogenetic-generalized-linear-mixed-models": ScientificValidationStudyDefinition(
        study_id="cichlid-phylogenetic-generalized-linear-mixed-models",
        verdict_summary=(
            "This bundle is checked in and reviewable, but the repository does"
            " not yet claim governed parity for the PCM5 generalized posterior workflow."
        ),
        comparison_mode="direct_parity",
        current_repository_posture=(
            "generalized_phylogenetic_mixed_model_parity_open_with_partial_binary_regression_support"
        ),
        governed_code_paths=(
            SCIENTIFIC_VALIDATION_BUNDLE_GENERATOR_PATH,
            "packages/bijux-phylogenetics/src/bijux_phylogenetics/comparative",
            "packages/bijux-phylogenetics/src/bijux_phylogenetics/phylo",
        ),
        source_basis_locators=(
            "evidence-book/studies/cichlid-phylogenetic-generalized-linear-mixed-models/provenance/lund-course-sources.json",
            "evidence-book/studies/cichlid-phylogenetic-generalized-linear-mixed-models/datasets/registry.json",
            "evidence-book/studies/cichlid-phylogenetic-generalized-linear-mixed-models/reference/README.md",
            "packages/bijux-phylogenetics/tests/test_phylogenetic_logistic.py",
        ),
        bundle_definitions=(
            ScientificValidationBundleDefinition(
                evidence_id="evidence-001",
                evidence_title="Cichlid repeated-measure workspace contract bundle",
                summary=(
                    "Governed study slice for the PCM5 fish workspace reload and"
                    " repeated-measure data contract."
                ),
                claim_id="pcm5-repeated-measure-workspace-contract",
                analytical_surfaces=("workflow-contracts", "repeated-measure-data"),
                claim_tags=(
                    "scientific-validation",
                    "cichlid",
                    "workflow",
                    "repeated-measures",
                ),
                source_fragments=(
                    "fish-workspace-reload",
                    "repeated-measure-data-contract",
                ),
                limitations=(
                    "This bundle records the checked-in study contract and does not yet claim lecture-specific parity for the derived fish workspace.",
                ),
                next_step=(
                    "Check in governed reload and species-observation contract"
                    " ledgers for the repeated-measure study workspace."
                ),
            ),
            ScientificValidationBundleDefinition(
                evidence_id="evidence-002",
                evidence_title="Cichlid posterior phylogenetic signal bundle",
                summary=(
                    "Governed study slice for the PCM5 posterior-derived"
                    " phylogenetic signal calculation on the Gaussian mixed model."
                ),
                claim_id="pcm5-posterior-phylogenetic-signal",
                analytical_surfaces=("phylogenetic-signal", "posterior-signal"),
                claim_tags=(
                    "scientific-validation",
                    "cichlid",
                    "phylogenetic-signal",
                    "posterior",
                ),
                source_fragments=("posterior-phylogenetic-signal",),
                limitations=(
                    "This bundle keeps the posterior signal calculation explicit without claiming repository-owned parity for the variance-component workflow.",
                ),
                next_step=(
                    "Add governed variance-component and posterior-signal ledgers"
                    " with the lecture denominator rules preserved."
                ),
            ),
            ScientificValidationBundleDefinition(
                evidence_id="evidence-003",
                evidence_title="Cichlid ancestral posterior export bundle",
                summary=(
                    "Governed study slice for the PCM5 ancestral-state posterior"
                    " extraction and export rules."
                ),
                claim_id="pcm5-ancestral-posterior-export",
                analytical_surfaces=("ancestral-posteriors", "posterior-export"),
                claim_tags=(
                    "scientific-validation",
                    "cichlid",
                    "ancestral-states",
                    "posterior-export",
                ),
                source_fragments=(
                    "ancestral-posterior-extraction",
                    "ancestral-export-rules",
                ),
                limitations=(
                    "The repository does not yet own the lecture ancestral posterior workflow, so this bundle remains a governed boundary record.",
                ),
                next_step=(
                    "Check in governed node-posterior ledgers and export rules for"
                    " the ancestral-state workflow."
                ),
            ),
            ScientificValidationBundleDefinition(
                evidence_id="evidence-004",
                evidence_title="Cichlid ancestral plotting contract bundle",
                summary=(
                    "Governed study slice for the PCM5 tree-linked ancestral"
                    " plotting contract."
                ),
                claim_id="pcm5-tree-linked-ancestral-plotting-contract",
                analytical_surfaces=("ancestral-plotting", "tree-linked-reporting"),
                claim_tags=(
                    "scientific-validation",
                    "cichlid",
                    "ancestral-plotting",
                    "tree-linked-reporting",
                ),
                source_fragments=("tree-linked-ancestral-plotting",),
                limitations=(
                    "This bundle keeps the plotting contract visible without claiming rendered-figure parity for the lecture visuals.",
                ),
                next_step=(
                    "Add governed plot-input ledgers and node-mapping records for"
                    " the ancestral plotting surface."
                ),
            ),
            ScientificValidationBundleDefinition(
                evidence_id="evidence-005",
                evidence_title="Cichlid Poisson count model bundle",
                summary=(
                    "Governed study slice for the PCM5 Poisson phylogenetic count"
                    " model of parasite abundance."
                ),
                claim_id="pcm5-poisson-phylogenetic-count-model",
                analytical_surfaces=("poisson-model", "parasite-abundance"),
                claim_tags=(
                    "scientific-validation",
                    "cichlid",
                    "poisson-model",
                    "parasite-abundance",
                ),
                source_fragments=("poisson-count-model",),
                limitations=(
                    "The repository does not yet own the lecture Poisson mixed-model posterior workflow, so this bundle remains boundary-first.",
                ),
                next_step=(
                    "Check in governed count-model posterior ledgers with explicit"
                    " link-function, parameter, and uncertainty summaries."
                ),
            ),
            ScientificValidationBundleDefinition(
                evidence_id="evidence-006",
                evidence_title="Cichlid multinomial infection summary bundle",
                summary=(
                    "Governed study slice for the PCM5 multinomial infection"
                    " model and its logit-scale summaries."
                ),
                claim_id="pcm5-multinomial-infection-logit-summaries",
                analytical_surfaces=("multinomial-model", "logit-summaries"),
                claim_tags=(
                    "scientific-validation",
                    "cichlid",
                    "multinomial-model",
                    "logit-scale",
                ),
                source_fragments=(
                    "multinomial-infection-model",
                    "logit-scale-summaries",
                ),
                limitations=(
                    "This bundle makes the multinomial summary boundary explicit without claiming repository-owned parity for the lecture posterior fit.",
                ),
                next_step=(
                    "Add governed logit-scale posterior ledgers and model summary"
                    " tables for the multinomial infection workflow."
                ),
            ),
            ScientificValidationBundleDefinition(
                evidence_id="evidence-007",
                evidence_title="Cichlid probability-scale transform bundle",
                summary=(
                    "Governed study slice for the PCM5 probability-scale"
                    " back-transformation of intercepts and slopes."
                ),
                claim_id="pcm5-probability-scale-back-transforms",
                analytical_surfaces=("probability-transforms",),
                claim_tags=(
                    "scientific-validation",
                    "cichlid",
                    "probability-scale",
                    "back-transform",
                ),
                source_fragments=("probability-scale-back-transforms",),
                limitations=(
                    "The repository does not yet claim governed parity for the lecture’s posterior distribution transforms, so this bundle remains a review boundary.",
                ),
                next_step=(
                    "Check in governed probability-scale transform ledgers tied"
                    " to the corresponding logit-scale posterior summaries."
                ),
            ),
            ScientificValidationBundleDefinition(
                evidence_id="evidence-008",
                evidence_title="Cichlid repeated-measure random effects bundle",
                summary=(
                    "Governed study slice for the PCM5 repeated-measure random"
                    " effects that combine species and treetip structure."
                ),
                claim_id="pcm5-repeated-measure-random-effects",
                analytical_surfaces=("repeated-measure-random-effects",),
                claim_tags=(
                    "scientific-validation",
                    "cichlid",
                    "repeated-measures",
                    "random-effects",
                ),
                source_fragments=("species-plus-treetip-random-effects",),
                limitations=(
                    "This bundle keeps the repeated-measure random-effect question explicit without claiming repository-owned parity for the lecture’s mixed-model sampler.",
                ),
                next_step=(
                    "Add governed repeated-measure random-effect ledgers with"
                    " species-level and treetip-level variance accounting."
                ),
            ),
            ScientificValidationBundleDefinition(
                evidence_id="evidence-009",
                evidence_title="Cichlid variance partition review bundle",
                summary=(
                    "Governed study slice for the PCM5 variance partition"
                    " between phylogeny, species, and residual components."
                ),
                claim_id="pcm5-variance-partition-review",
                analytical_surfaces=("variance-partition", "posterior-variance"),
                claim_tags=(
                    "scientific-validation",
                    "cichlid",
                    "variance-partition",
                    "posterior-variance",
                ),
                source_fragments=("variance-partition-review",),
                limitations=(
                    "The repository does not yet claim parity for the lecture variance-partition posteriors, so this bundle stays a governed review slice.",
                ),
                next_step=(
                    "Check in governed variance partition ledgers with explicit"
                    " denominator rules and posterior interval summaries."
                ),
            ),
        ),
    ),
}


def scientific_validation_study_ids() -> tuple[str, ...]:
    return tuple(sorted(SCIENTIFIC_VALIDATION_STUDIES))


def scientific_validation_supported_evidence_ids(study_id: str) -> tuple[str, ...]:
    study = SCIENTIFIC_VALIDATION_STUDIES[study_id]
    return tuple(bundle.evidence_id for bundle in study.bundle_definitions)


def _render_bundle_source_basis(
    study_definition: ScientificValidationStudyDefinition,
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for index, locator in enumerate(study_definition.source_basis_locators, start=1):
        kind = "repository-reference"
        if locator.startswith("external:"):
            kind = "external-source-descriptor"
        rows.append(
            {
                "kind": kind,
                "label": f"source-{index}",
                "locator": locator,
            }
        )
    return rows


def _render_manifest(
    study_definition: ScientificValidationStudyDefinition,
    bundle_definition: ScientificValidationBundleDefinition,
) -> dict[str, object]:
    return {
        "schema_version": 1,
        "study_id": study_definition.study_id,
        "evidence_id": bundle_definition.evidence_id,
        "evidence_title": bundle_definition.evidence_title,
        "summary": bundle_definition.summary,
        "owner_package": DEFAULT_OWNER_PACKAGE,
        "claim_ids": [bundle_definition.claim_id],
        "source_basis": _render_bundle_source_basis(study_definition),
        "freshness": {
            "last_generated_on": "2026-07-04",
            "governed_code_paths": list(study_definition.governed_code_paths),
            "source_basis_locators": list(study_definition.source_basis_locators),
        },
        "ownership": {
            "owner_package": DEFAULT_OWNER_PACKAGE,
            "analytical_surfaces": list(bundle_definition.analytical_surfaces),
        },
        "claim_tags": list(bundle_definition.claim_tags),
        "comparison_mode": study_definition.comparison_mode,
        "verdict": {
            "status": "not_comparable",
            "summary": study_definition.verdict_summary,
        },
        "limitations": list(bundle_definition.limitations),
    }


def _render_claims_payload(
    study_definition: ScientificValidationStudyDefinition,
    bundle_definition: ScientificValidationBundleDefinition,
) -> dict[str, object]:
    return {
        "schema_version": 1,
        "study_id": study_definition.study_id,
        "evidence_id": bundle_definition.evidence_id,
        "claim_count": 1,
        "claims": [
            {
                "claim_id": bundle_definition.claim_id,
                "claim_title": bundle_definition.evidence_title,
                "summary": bundle_definition.summary,
                "verdict": "not_comparable",
                "evidence_ids": [bundle_definition.evidence_id],
                "source_fragments": list(bundle_definition.source_fragments),
            }
        ],
    }


def _render_bundle_readme(
    study_title: str,
    study_definition: ScientificValidationStudyDefinition,
    bundle_definition: ScientificValidationBundleDefinition,
) -> str:
    lines = [
        f"# {bundle_definition.evidence_title}",
        "",
        bundle_definition.summary,
        "",
        f"- study: `{study_title}`",
        f"- study id: `{study_definition.study_id}`",
        f"- evidence id: `{bundle_definition.evidence_id}`",
        f"- comparison mode: `{study_definition.comparison_mode}`",
        "- current verdict: `not_comparable`",
        "",
        "## Why This Bundle Exists",
        "",
        study_definition.verdict_summary,
        "",
        "## Source Fragments",
        "",
    ]
    for fragment in bundle_definition.source_fragments:
        lines.append(f"- `{fragment}`")
    lines.extend(
        [
            "",
            "## Current Repository Posture",
            "",
            f"- `{study_definition.current_repository_posture}`",
            "",
            "## Source Basis",
            "",
        ]
    )
    for locator in study_definition.source_basis_locators:
        lines.append(f"- `{locator}`")
    lines.extend(
        [
            "",
            "## Current Limits",
            "",
        ]
    )
    for item in bundle_definition.limitations:
        lines.append(f"- {item}")
    lines.extend(
        [
            "",
            "## Next Step",
            "",
            f"- {bundle_definition.next_step}",
            "",
        ]
    )
    return "\n".join(lines)


def _write_json(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def _write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def write_scientific_validation_bundle(
    repo_root: Path | str,
    study_id: str,
    evidence_id: str,
) -> list[str]:
    repo_root = Path(repo_root)
    study_definition = SCIENTIFIC_VALIDATION_STUDIES[study_id]
    bundle_definition = next(
        bundle
        for bundle in study_definition.bundle_definitions
        if bundle.evidence_id == evidence_id
    )
    study_root = repo_root / "evidence-book" / "studies" / study_id
    study_title = str(load_study_contract(study_root)["study_title"])
    bundle_root = study_root / evidence_id
    bundle_root.mkdir(parents=True, exist_ok=True)

    _write_json(
        bundle_root / "manifest.json",
        _render_manifest(study_definition, bundle_definition),
    )
    _write_json(
        bundle_root / "claims.json",
        _render_claims_payload(study_definition, bundle_definition),
    )
    _write_text(
        bundle_root / "README.md",
        _render_bundle_readme(study_title, study_definition, bundle_definition),
    )
    for relative_path, payload in build_bundle_governed_artifacts(
        repo_root, bundle_root
    ).items():
        target = bundle_root / relative_path
        if isinstance(payload, dict):
            _write_json(target, payload)
        else:
            _write_text(target, payload)
    return sorted(
        str(path.relative_to(repo_root))
        for path in bundle_root.rglob("*")
        if path.is_file()
    )


def write_scientific_validation_study(
    repo_root: Path | str,
    study_id: str,
    evidence_ids: list[str] | None = None,
) -> list[str]:
    selected = (
        scientific_validation_supported_evidence_ids(study_id)
        if evidence_ids is None
        else tuple(sorted(set(evidence_ids)))
    )
    updated_paths: list[str] = []
    for evidence_id in selected:
        updated_paths.extend(
            write_scientific_validation_bundle(repo_root, study_id, evidence_id)
        )
    return sorted(set(updated_paths))
