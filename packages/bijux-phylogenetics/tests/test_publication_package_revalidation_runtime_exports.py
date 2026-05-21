from __future__ import annotations

from bijux_phylogenetics.reports import (
    PublicationPackageRevalidationArtifactRow,
    PublicationPackageRevalidationCheckRow,
    PublicationPackageRevalidationResult,
    write_publication_package_revalidation_report,
)
import bijux_phylogenetics.reports.publication.revalidation as package_revalidation_api


def test_publication_package_revalidation_runtime_exports_are_public() -> None:
    assert (
        package_revalidation_api.PublicationPackageRevalidationArtifactRow
        is PublicationPackageRevalidationArtifactRow
    )
    assert (
        package_revalidation_api.PublicationPackageRevalidationCheckRow
        is PublicationPackageRevalidationCheckRow
    )
    assert (
        package_revalidation_api.PublicationPackageRevalidationResult
        is PublicationPackageRevalidationResult
    )
    assert (
        package_revalidation_api.write_publication_package_revalidation_report
        is write_publication_package_revalidation_report
    )
