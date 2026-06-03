"""Runtime identity, error, and result contracts for Bijux Phylogenetics."""

from typing import Any

__all__ = [
    "AlignmentTaxonMismatchError",
    "AncestralReconstructionError",
    "CLI_ALIASES",
    "CLI_NAME",
    "CommandResult",
    "ComparativeMethodError",
    "DiversificationAnalysisError",
    "DuplicateTaxonError",
    "EngineUnavailableError",
    "EngineWorkflowError",
    "EvidenceContractError",
    "IDENTITY",
    "IMPORT_NAME",
    "InvalidAlignmentError",
    "InvalidBranchLengthError",
    "InvalidDistanceMatrixError",
    "InvalidPartitionError",
    "MetadataJoinError",
    "NonUltrametricTreeError",
    "PACKAGE_NAME",
    "PRODUCT_NAME",
    "PackageIdentity",
    "PhylogeneticsError",
    "ScientificFailureExplanation",
    "TreeParseError",
    "TreeRootingError",
    "UMBRELLA_COMMAND",
    "UnnamedTipError",
    "UnsupportedDistanceTreeMethodError",
    "UnsupportedTreeFormatError",
    "UnrootedTreeError",
    "WorkflowBudgetError",
    "build_command_result",
    "build_error_result",
    "explain_inference_workflow_failure",
    "explain_phylogenetics_error",
    "explanation_payload",
]


def __getattr__(name: str):
    if name in {
        "AlignmentTaxonMismatchError",
        "AncestralReconstructionError",
        "ComparativeMethodError",
        "DiversificationAnalysisError",
        "DuplicateTaxonError",
        "EngineUnavailableError",
        "EngineWorkflowError",
        "EvidenceContractError",
        "InvalidAlignmentError",
        "InvalidBranchLengthError",
        "InvalidDistanceMatrixError",
        "InvalidPartitionError",
        "MetadataJoinError",
        "NonUltrametricTreeError",
        "PhylogeneticsError",
        "TreeParseError",
        "TreeRootingError",
        "UnnamedTipError",
        "UnsupportedDistanceTreeMethodError",
        "UnsupportedTreeFormatError",
        "UnrootedTreeError",
        "WorkflowBudgetError",
    }:
        from .errors import (
            AlignmentTaxonMismatchError,
            AncestralReconstructionError,
            ComparativeMethodError,
            DiversificationAnalysisError,
            DuplicateTaxonError,
            EngineUnavailableError,
            EngineWorkflowError,
            EvidenceContractError,
            InvalidAlignmentError,
            InvalidBranchLengthError,
            InvalidDistanceMatrixError,
            InvalidPartitionError,
            MetadataJoinError,
            NonUltrametricTreeError,
            PhylogeneticsError,
            TreeParseError,
            TreeRootingError,
            UnnamedTipError,
            UnrootedTreeError,
            UnsupportedDistanceTreeMethodError,
            UnsupportedTreeFormatError,
            WorkflowBudgetError,
        )

        exports = {
            "AlignmentTaxonMismatchError": AlignmentTaxonMismatchError,
            "AncestralReconstructionError": AncestralReconstructionError,
            "ComparativeMethodError": ComparativeMethodError,
            "DiversificationAnalysisError": DiversificationAnalysisError,
            "DuplicateTaxonError": DuplicateTaxonError,
            "EngineUnavailableError": EngineUnavailableError,
            "EngineWorkflowError": EngineWorkflowError,
            "EvidenceContractError": EvidenceContractError,
            "InvalidAlignmentError": InvalidAlignmentError,
            "InvalidBranchLengthError": InvalidBranchLengthError,
            "InvalidDistanceMatrixError": InvalidDistanceMatrixError,
            "InvalidPartitionError": InvalidPartitionError,
            "MetadataJoinError": MetadataJoinError,
            "NonUltrametricTreeError": NonUltrametricTreeError,
            "PhylogeneticsError": PhylogeneticsError,
            "TreeParseError": TreeParseError,
            "TreeRootingError": TreeRootingError,
            "UnnamedTipError": UnnamedTipError,
            "UnsupportedDistanceTreeMethodError": UnsupportedDistanceTreeMethodError,
            "UnsupportedTreeFormatError": UnsupportedTreeFormatError,
            "UnrootedTreeError": UnrootedTreeError,
            "WorkflowBudgetError": WorkflowBudgetError,
        }
        return exports[name]
    if name in {
        "CLI_ALIASES",
        "CLI_NAME",
        "IDENTITY",
        "IMPORT_NAME",
        "PACKAGE_NAME",
        "PRODUCT_NAME",
        "PackageIdentity",
        "UMBRELLA_COMMAND",
    }:
        from .identity import (
            CLI_ALIASES,
            CLI_NAME,
            IDENTITY,
            IMPORT_NAME,
            PACKAGE_NAME,
            PRODUCT_NAME,
            UMBRELLA_COMMAND,
            PackageIdentity,
        )

        exports: dict[str, Any] = {
            "CLI_ALIASES": CLI_ALIASES,
            "CLI_NAME": CLI_NAME,
            "IDENTITY": IDENTITY,
            "IMPORT_NAME": IMPORT_NAME,
            "PACKAGE_NAME": PACKAGE_NAME,
            "PRODUCT_NAME": PRODUCT_NAME,
            "PackageIdentity": PackageIdentity,
            "UMBRELLA_COMMAND": UMBRELLA_COMMAND,
        }
        return exports[name]
    if name in {
        "CommandResult",
        "build_command_result",
        "build_error_result",
    }:
        from .results import CommandResult, build_command_result, build_error_result

        exports: dict[str, Any] = {
            "CommandResult": CommandResult,
            "build_command_result": build_command_result,
            "build_error_result": build_error_result,
        }
        return exports[name]
    if name in {
        "ScientificFailureExplanation",
        "explain_inference_workflow_failure",
        "explain_phylogenetics_error",
        "explanation_payload",
    }:
        from .error_explanations import (
            ScientificFailureExplanation,
            explain_inference_workflow_failure,
            explain_phylogenetics_error,
            explanation_payload,
        )

        exports: dict[str, Any] = {
            "ScientificFailureExplanation": ScientificFailureExplanation,
            "explain_inference_workflow_failure": explain_inference_workflow_failure,
            "explain_phylogenetics_error": explain_phylogenetics_error,
            "explanation_payload": explanation_payload,
        }
        return exports[name]
    raise AttributeError(name)
