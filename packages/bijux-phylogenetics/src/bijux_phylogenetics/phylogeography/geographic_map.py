from __future__ import annotations

from bijux_phylogenetics.phylogeography.geographic_map_review.artifact_outputs import (
    write_geographic_map_exclusion_table as _write_geographic_map_exclusion_table,
    write_geographic_map_line_table as _write_geographic_map_line_table,
    write_geographic_map_marker_table as _write_geographic_map_marker_table,
    write_geographic_map_summary_table as _write_geographic_map_summary_table,
)
from bijux_phylogenetics.phylogeography.geographic_map_review.continuous_review import (
    summarize_continuous_phylogeography_map as _summarize_continuous_phylogeography_map,
)
from bijux_phylogenetics.phylogeography.geographic_map_review.discrete_review import (
    summarize_discrete_region_map as _summarize_discrete_region_map,
)
from bijux_phylogenetics.phylogeography.geographic_map_review.presentation import (
    render_geographic_map_html as _render_geographic_map_html,
)
from bijux_phylogenetics.phylogeography.geographic_map_review.contracts import (
    GeographicMapArtifact as _GeographicMapArtifact,
    GeographicMapExclusionRow as _GeographicMapExclusionRow,
    GeographicMapLineRow as _GeographicMapLineRow,
    GeographicMapMarkerRow as _GeographicMapMarkerRow,
    GeographicMapReport as _GeographicMapReport,
    GeographicMapSummary as _GeographicMapSummary,
)

GeographicMapArtifact = _GeographicMapArtifact
GeographicMapReport = _GeographicMapReport
render_geographic_map_html = _render_geographic_map_html
summarize_continuous_phylogeography_map = _summarize_continuous_phylogeography_map
summarize_discrete_region_map = _summarize_discrete_region_map
GeographicMapExclusionRow = _GeographicMapExclusionRow
GeographicMapLineRow = _GeographicMapLineRow
GeographicMapMarkerRow = _GeographicMapMarkerRow
GeographicMapSummary = _GeographicMapSummary
write_geographic_map_summary_table = _write_geographic_map_summary_table
write_geographic_map_marker_table = _write_geographic_map_marker_table
write_geographic_map_line_table = _write_geographic_map_line_table
write_geographic_map_exclusion_table = _write_geographic_map_exclusion_table
