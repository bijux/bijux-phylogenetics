from __future__ import annotations

import re

from bijux_phylogenetics.io.newick import quote_newick_label

_TRANSLATE_MARKER_PATTERN = re.compile(r"\btranslate\b", flags=re.IGNORECASE)


def parse_nexus_translate_map(text: str) -> dict[str, str]:
    """Parse a NEXUS translate block into resolved taxon labels."""
    block = _extract_nexus_translate_block(text)
    if block is None:
        return {}
    if re.search(r"(^|\n)\s*(tree|begin|end)\b", block, flags=re.IGNORECASE):
        raise ValueError("nexus translate block extends into another section")
    mapping: dict[str, str] = {}
    for entry in _split_nexus_translate_entries(block):
        parts = entry.split(None, 1)
        if len(parts) != 2:
            continue
        key, value = parts
        mapping[key.strip()] = _unquote_nexus_label(value.strip().rstrip(","))
    return mapping


def translate_nexus_tip_labels(newick: str, mapping: dict[str, str]) -> str:
    """Resolve translated NEXUS tip tokens inside one Newick statement."""
    if not mapping:
        return newick

    def replace(match: re.Match[str]) -> str:
        token = match.group(1)
        translated = mapping.get(token)
        if translated is None:
            return match.group(0)
        return match.group(0).replace(token, quote_newick_label(translated))

    return re.sub(r"(?<=[(,])\s*([A-Za-z0-9_.-]+)(?=\s*[:),])", replace, newick)


def _extract_nexus_translate_block(text: str) -> str | None:
    marker = _TRANSLATE_MARKER_PATTERN.search(text)
    if marker is None:
        return None
    remainder = text[marker.end() :]
    current: list[str] = []
    quote_character: str | None = None
    index = 0
    while index < len(remainder):
        character = remainder[index]
        if quote_character is None and character == ";":
            return "".join(current)
        if character in {"'", '"'}:
            if quote_character is None:
                quote_character = character
            elif character == quote_character:
                if (
                    quote_character == "'"
                    and index + 1 < len(remainder)
                    and remainder[index + 1] == "'"
                ):
                    current.append("''")
                    index += 2
                    continue
                quote_character = None
        current.append(character)
        index += 1
    raise ValueError("nexus translate block is unterminated")


def _split_nexus_translate_entries(raw_block: str) -> list[str]:
    entries: list[str] = []
    current: list[str] = []
    quote_character: str | None = None
    index = 0
    while index < len(raw_block):
        character = raw_block[index]
        if character in {"'", '"'}:
            if quote_character is None:
                quote_character = character
            elif character == quote_character:
                if (
                    quote_character == "'"
                    and index + 1 < len(raw_block)
                    and raw_block[index + 1] == "'"
                ):
                    current.append("''")
                    index += 2
                    continue
                quote_character = None
        if character == "," and quote_character is None:
            candidate = "".join(current).strip()
            if candidate:
                entries.append(candidate)
            current = []
            index += 1
            continue
        current.append(character)
        index += 1
    tail = "".join(current).strip()
    if tail:
        entries.append(tail)
    return entries


def _unquote_nexus_label(value: str) -> str:
    text = value.strip()
    if len(text) >= 2 and text[0] == text[-1] and text[0] in {"'", '"'}:
        unquoted = text[1:-1]
        if text[0] == "'":
            return unquoted.replace("''", "'")
        return unquoted
    return text
