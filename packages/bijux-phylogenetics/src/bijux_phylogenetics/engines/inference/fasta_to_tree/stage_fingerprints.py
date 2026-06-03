from __future__ import annotations

import hashlib
import json
from pathlib import Path

from .contracts import FastaToTreeStageFingerprint


def normalize_stage_fingerprint_payload(value: object) -> object:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {
            str(key): normalize_stage_fingerprint_payload(item)
            for key, item in sorted(value.items())
        }
    if isinstance(value, (list, tuple)):
        return [normalize_stage_fingerprint_payload(item) for item in value]
    return value


def build_stage_fingerprint(
    *,
    stage: str,
    input_checksums: dict[str, str],
    output_checksums: dict[str, str],
    config: dict[str, object],
    engine_versions: dict[str, str],
    upstream_fingerprints: dict[str, str],
    resumed: bool,
) -> FastaToTreeStageFingerprint:
    payload = normalize_stage_fingerprint_payload(
        {
            "stage": stage,
            "input_checksums": input_checksums,
            "output_checksums": output_checksums,
            "config": config,
            "engine_versions": engine_versions,
            "upstream_fingerprints": upstream_fingerprints,
        }
    )
    digest = hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()
    return FastaToTreeStageFingerprint(
        stage=stage,
        fingerprint=digest,
        input_checksums=input_checksums,
        output_checksums=output_checksums,
        config=config,
        engine_versions=engine_versions,
        upstream_fingerprints=upstream_fingerprints,
        resumed=resumed,
    )
