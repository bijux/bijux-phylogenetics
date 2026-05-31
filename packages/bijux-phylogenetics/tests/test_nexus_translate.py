from __future__ import annotations

from bijux_phylogenetics.io.nexus_translate import (
    parse_nexus_translate_map,
    translate_nexus_tip_labels,
)


def test_parse_nexus_translate_map_unquotes_labels_with_spaces_and_apostrophes() -> (
    None
):
    mapping = parse_nexus_translate_map(
        "#NEXUS\n"
        "Begin trees;\n"
        "  Translate\n"
        "    1 'Taxon A',\n"
        "    2 'Taxon B''s sample',\n"
        '    3 "Taxon C"\n'
        "  ;\n"
        "End;\n"
    )

    assert mapping == {
        "1": "Taxon A",
        "2": "Taxon B's sample",
        "3": "Taxon C",
    }


def test_translate_nexus_tip_labels_requotes_resolved_taxa_for_newick() -> None:
    translated = translate_nexus_tip_labels(
        "((1:0.1,2:0.1):0.2,3:0.3)",
        {
            "1": "Taxon A",
            "2": "Taxon B's sample",
            "3": "Taxon_C",
        },
    )

    assert translated == "(('Taxon A':0.1,'Taxon B''s sample':0.1):0.2,Taxon_C:0.3)"
