from __future__ import annotations

import bijux_phylogenetics.reports.package_comparison as package_comparison_api
from bijux_phylogenetics.reports import (
    PublicationPackageComparisonArtifactRow,
    PublicationPackageComparisonCheckRow,
    PublicationPackageComparisonResult,
    write_publication_package_comparison_report,
)


def test_publication_package_comparison_runtime_exports_are_public() -> None:
    assert (
        package_comparison_api.PublicationPackageComparisonArtifactRow
        is PublicationPackageComparisonArtifactRow
    )
    assert (
        package_comparison_api.PublicationPackageComparisonCheckRow
        is PublicationPackageComparisonCheckRow
    )
    assert (
        package_comparison_api.PublicationPackageComparisonResult
        is PublicationPackageComparisonResult
    )
    assert (
        package_comparison_api.write_publication_package_comparison_report
        is write_publication_package_comparison_report
    )
