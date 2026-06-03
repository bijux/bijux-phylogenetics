from __future__ import annotations

from .._shared import (
    _XML_IDENTIFIER_PATTERN,
    ElementTree,
    Path,
    XmlElement,
    load_tree,
)


def _xml_element(
    tag: str,
    attributes: dict[str, str] | None = None,
    *,
    text: str | None = None,
    children: tuple[XmlElement, ...] = (),
) -> XmlElement:
    """Build one trusted XML element for BEAST output assembly."""
    element = ElementTree.Element(tag, attributes or {})
    if text is not None:
        element.text = text
    for child in children:
        element.append(child)
    return element


def _xml_identifier(raw: str) -> str:
    normalized = _XML_IDENTIFIER_PATTERN.sub("_", raw.strip())
    normalized = normalized.strip("_")
    return normalized or "identifier"


def _beast_data_type(inferred_alphabet: str) -> str:
    if inferred_alphabet in {"dna", "rna"}:
        return "nucleotide"
    if inferred_alphabet == "protein":
        return "aminoacid"
    raise ValueError(
        "BEAST preparation requires a nucleotide, RNA, or protein alignment"
    )


def _default_beast_substitution_model(beast_data_type: str) -> str:
    if beast_data_type == "nucleotide":
        return "HKY"
    return "JTT"


def _read_newick_text(path: Path) -> str:
    raw = path.read_text(encoding="utf-8").strip()
    if not raw:
        raise ValueError(f"tree file is empty: {path}")
    first_line = raw.splitlines()[0].strip()
    if not first_line.endswith(";"):
        first_line = f"{first_line};"
    return first_line


def _validate_tree_taxa_against_alignment(
    *, tree_path: Path, alignment_taxa: set[str]
) -> None:
    tree_taxa = set(load_tree(tree_path).tip_names)
    missing_from_tree = sorted(alignment_taxa - tree_taxa)
    extra_in_tree = sorted(tree_taxa - alignment_taxa)
    if missing_from_tree or extra_in_tree:
        details: list[str] = []
        if missing_from_tree:
            details.append(
                "alignment taxa missing from tree: " + ", ".join(missing_from_tree)
            )
        if extra_in_tree:
            details.append(
                "tree taxa missing from alignment: " + ", ".join(extra_in_tree)
            )
        raise ValueError(
            "BEAST preparation requires the starting tree and alignment to contain the same taxa: "
            + "; ".join(details)
        )
