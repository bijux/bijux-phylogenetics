from __future__ import annotations


class PhylogeneticsError(Exception):
    """Base exception carrying a stable machine-readable error code."""

    code = "phylogenetics_error"

    def __init__(self, message: str) -> None:
        super().__init__(message)
        self.message = message


class TreeParseError(PhylogeneticsError):
    code = "tree_parse_error"


class UnsupportedTreeFormatError(PhylogeneticsError):
    code = "unsupported_tree_format_error"


class InvalidBranchLengthError(PhylogeneticsError):
    code = "invalid_branch_length_error"


class DuplicateTaxonError(PhylogeneticsError):
    code = "duplicate_taxon_error"


class UnnamedTipError(PhylogeneticsError):
    code = "unnamed_tip_error"


class MetadataJoinError(PhylogeneticsError):
    code = "metadata_join_error"


class InvalidAlignmentError(PhylogeneticsError):
    code = "invalid_alignment_error"


class InvalidDistanceMatrixError(PhylogeneticsError):
    code = "invalid_distance_matrix_error"


class AlignmentTaxonMismatchError(PhylogeneticsError):
    code = "alignment_taxon_mismatch_error"


class NonUltrametricTreeError(PhylogeneticsError):
    code = "non_ultrametric_tree_error"


class UnrootedTreeError(PhylogeneticsError):
    code = "unrooted_tree_error"


class EngineUnavailableError(PhylogeneticsError):
    code = "engine_unavailable_error"


class EvidenceContractError(PhylogeneticsError):
    code = "evidence_contract_error"


class ComparativeMethodError(PhylogeneticsError):
    code = "comparative_method_error"


class AncestralReconstructionError(PhylogeneticsError):
    code = "ancestral_reconstruction_error"
