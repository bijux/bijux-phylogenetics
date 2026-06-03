from __future__ import annotations

from argparse import ArgumentParser

PAIRWISE_DISTANCE_MODELS = (
    "raw",
    "p-distance",
    "jc69",
    "jukes-cantor",
    "k80",
    "kimura-2-parameter",
    "f81",
    "felsenstein-81",
    "tn93",
    "tamura-nei-93",
    "amino-acid-p-distance",
)


def add_distance_model_option(parser: ArgumentParser) -> None:
    parser.add_argument(
        "--model",
        choices=PAIRWISE_DISTANCE_MODELS,
        default="p-distance",
    )


def add_gap_handling_option(parser: ArgumentParser) -> None:
    parser.add_argument(
        "--gap-handling",
        choices=("pairwise-deletion", "complete-deletion"),
        default="pairwise-deletion",
    )


def add_ambiguity_policy_option(parser: ArgumentParser) -> None:
    parser.add_argument(
        "--ambiguity-policy",
        choices=("ignore", "partial-match", "strict-mismatch", "report-only"),
        default="ignore",
    )
