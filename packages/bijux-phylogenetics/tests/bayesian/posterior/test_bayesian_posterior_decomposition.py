from __future__ import annotations

import csv
from pathlib import Path

import pytest

from bijux_phylogenetics.bayesian import (
    summarize_beast_posterior_decomposition,
    summarize_mrbayes_posterior_decomposition,
    write_beast_posterior_decomposition_table,
    write_mrbayes_posterior_decomposition_table,
)
from bijux_phylogenetics.runtime.errors import EngineWorkflowError

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_summarize_beast_posterior_decomposition_verifies_logged_prior_identity(
    tmp_path: Path,
) -> None:
    log_path = tmp_path / "classified-beast.log"
    output_path = tmp_path / "classified-beast-decomposition.tsv"
    log_path.write_text(
        "# posterior decomposition fixture\n"
        "state\tposterior\tlikelihood\tprior\tclockRate\n"
        "0\t-510.0\t-490.0\t-20.0\t0.0010\n"
        "1000\t-505.0\t-486.0\t-19.0\t0.0011\n"
        "2000\t-500.0\t-482.0\t-18.0\t0.0012\n"
        "3000\t-497.0\t-479.0\t-18.0\t0.0013\n",
        encoding="utf-8",
    )

    report = summarize_beast_posterior_decomposition(log_path, burnin_fraction=0.25)
    write_beast_posterior_decomposition_table(output_path, report)

    assert report.kept_row_count == 3
    assert report.first_kept_state == 1000
    assert report.last_kept_state == 3000
    assert report.posterior_term_source == "logged"
    assert report.likelihood_term_source == "logged"
    assert report.prior_term_source == "logged"
    assert report.verified is True
    assert report.maximum_absolute_delta == pytest.approx(0.0)
    assert report.rows[0].log_prior == pytest.approx(-19.0)
    with output_path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert rows[0]["state"] == "1000"
    assert rows[0]["decomposition_valid"] == "true"


def test_summarize_beast_posterior_decomposition_derives_missing_prior_term(
    tmp_path: Path,
) -> None:
    log_path = tmp_path / "derived-prior-beast.log"
    log_path.write_text(
        "# missing-prior fixture\n"
        "state\tposterior\tlikelihood\n"
        "0\t-100.0\t-90.0\n"
        "10\t-99.0\t-88.5\n",
        encoding="utf-8",
    )

    report = summarize_beast_posterior_decomposition(log_path)

    assert report.prior_term_source == "derived_from_posterior_and_likelihood"
    assert report.rows[0].log_prior == pytest.approx(-10.0)
    assert report.rows[1].log_prior == pytest.approx(-10.5)
    assert report.verified is True


def test_summarize_mrbayes_posterior_decomposition_uses_logged_prior_term() -> None:
    report = summarize_mrbayes_posterior_decomposition(
        fixture("mrbayes", "partitioned-analysis.run1.p")
    )

    assert report.kept_row_count == 3
    assert report.first_kept_generation == 0
    assert report.last_kept_generation == 20
    assert report.posterior_term_source == "derived_from_likelihood_and_prior"
    assert report.likelihood_term_source == "logged"
    assert report.prior_term_source == "logged"
    assert report.verified is True
    assert report.maximum_absolute_delta == pytest.approx(0.0)
    assert report.rows[0].log_posterior == pytest.approx(0.02016)


def test_summarize_mrbayes_posterior_decomposition_rejects_traces_without_prior_term() -> (
    None
):
    with pytest.raises(EngineWorkflowError) as error:
        summarize_mrbayes_posterior_decomposition(
            fixture("engine_outputs", "mrbayes/trace-warning-heavy.run1.p")
        )

    assert error.value.code == "mrbayes_trace_missing_posterior_terms"
    assert error.value.details["artifact_kind"] == "mrbayes-trace"
    assert error.value.details["missing_columns"] == ["LnPr"]


def test_write_mrbayes_posterior_decomposition_table_writes_expected_header(
    tmp_path: Path,
) -> None:
    report = summarize_mrbayes_posterior_decomposition(
        fixture("mrbayes", "partitioned-analysis.run1.p")
    )
    output_path = tmp_path / "mrbayes-posterior-decomposition.tsv"

    write_mrbayes_posterior_decomposition_table(output_path, report)

    text = output_path.read_text(encoding="utf-8")
    assert text.startswith(
        "generation\tlog_posterior\tlog_likelihood\tlog_prior\tdecomposition_delta\tdecomposition_valid\tposterior_term_source\tlikelihood_term_source\tprior_term_source\tidentity_tolerance\n"
    )
    assert "derived_from_likelihood_and_prior" in text
