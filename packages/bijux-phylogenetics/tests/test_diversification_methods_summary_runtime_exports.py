from __future__ import annotations

import bijux_phylogenetics.comparative as comparative_api
from bijux_phylogenetics.comparative import (
    GeigerBirthDeathExclusionReport,
    MedusaExclusionReport,
    summarize_geiger_birth_death_exclusion,
    summarize_medusa_exclusion,
    write_diversification_methods_summary_text,
)


def test_diversification_methods_summary_surfaces_export_publicly() -> None:
    assert (
        comparative_api.GeigerBirthDeathExclusionReport
        is GeigerBirthDeathExclusionReport
    )
    assert (
        comparative_api.summarize_geiger_birth_death_exclusion
        is summarize_geiger_birth_death_exclusion
    )
    assert comparative_api.MedusaExclusionReport is MedusaExclusionReport
    assert comparative_api.summarize_medusa_exclusion is summarize_medusa_exclusion
    assert callable(write_diversification_methods_summary_text)
    assert (
        comparative_api.write_diversification_methods_summary_text
        is write_diversification_methods_summary_text
    )
