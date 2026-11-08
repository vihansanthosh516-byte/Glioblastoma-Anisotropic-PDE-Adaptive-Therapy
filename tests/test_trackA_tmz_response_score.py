from __future__ import annotations

import numpy as np
import pytest
import src.trackA_tmz_response_score as proxy

from src.trackA_tmz_response_score import (
    PRIMARY_DONORS,
    ProxyStop,
    compute_proxy_scores,
    exact_two_sided_sign_flip,
    generate_matched_control_panels,
    rank_percentiles,
    read_panel,
    score_rank_percentiles,
    sha256_file,
    verify_freeze,
)


def _synthetic_panel() -> list[dict[str, str]]:
    rows = []
    for i in range(121):
        rho = 0.5 if i < 48 else -0.5
        rows.append(
            {
                "gene_symbol": f"G{i:03d}",
                "rho_viability_100uM": str(rho),
                "p_value": "0.01",
                "p_adjusted": "0.05",
                "direction": "resistance_associated" if rho > 0 else "sensitivity_associated",
            }
        )
    return rows


def _synthetic_dataset():
    rng = np.random.default_rng(30)
    gene_ids = [f"G{i:03d}" for i in range(121)] + [f"BG{i:03d}" for i in range(180)]
    rows = []
    donors = []
    labels = []
    for donor in PRIMARY_DONORS:
        for label in (1, 0):
            for _ in range(20):
                values = rng.integers(5, 12, size=len(gene_ids), dtype=np.int32)
                if label == 1:
                    values[:48] += 30
                    values[48:121] = np.maximum(1, values[48:121] - 4)
                else:
                    values[:48] = np.maximum(1, values[:48] - 4)
                    values[48:121] += 30
                rows.append(values)
                donors.append(donor)
                labels.append(label)
    return np.asarray(rows), np.asarray(donors), np.asarray(labels), gene_ids


def test_official_panel_has_expected_size_and_directions():
    panel = read_panel()
    assert len(panel) == 121
    assert sum(float(row["rho_viability_100uM"]) > 0 for row in panel) == 48
    assert sum(float(row["rho_viability_100uM"]) < 0 for row in panel) == 73


def test_rank_score_increases_when_positive_group_ranks_higher():
    gene_ids = ["positive", "negative", "other"]
    ranks = rank_percentiles(np.asarray([3.0, 1.0, 2.0]))
    score = score_rank_percentiles(ranks, gene_ids, ["positive"], ["negative"])
    assert score == pytest.approx(1.0)


def test_average_ranks_handle_ties():
    ranks = rank_percentiles(np.asarray([1.0, 1.0, 3.0]))
    assert ranks == pytest.approx(np.asarray([0.25, 0.25, 1.0]))


def test_exact_sign_flip_uses_all_64_donor_sign_patterns():
    result = exact_two_sided_sign_flip([1, 1, 1, 1, 1, 1])
    assert result["sign_assignments"] == 64
    assert result["analysis_unit"] == "donor"
    assert result["p_value"] == pytest.approx(2 / 64)


def test_response_score_returns_six_paired_donor_results_and_matched_controls():
    counts, donors, labels, gene_ids = _synthetic_dataset()
    result = compute_proxy_scores(
        counts,
        donors,
        labels,
        gene_ids,
        _synthetic_panel(),
        n_control_panels=25,
        control_seed=20261107,
    )
    assert result["status"] == "complete_exploratory"
    assert result["cohort"]["n_donors"] == 6
    assert len(result["score"]["donor_results"]) == 6
    assert all(row["core_minus_periphery"] > 0 for row in result["score"]["donor_results"].values())
    assert result["matched_controls"]["n_panels"] == 25
    assert "not a significance test" in result["matched_controls"]["use"]


def test_missing_exact_symbol_stops_without_alias_substitution():
    counts, donors, labels, gene_ids = _synthetic_dataset()
    with pytest.raises(ProxyStop, match="Exact panel mapping failed"):
        compute_proxy_scores(counts[:, 1:], donors, labels, gene_ids[1:], _synthetic_panel(), n_control_panels=1)


def test_small_donor_region_stops():
    counts, donors, labels, gene_ids = _synthetic_dataset()
    keep = ~((donors == "GS1") & (labels == 1))
    keep[np.flatnonzero((donors == "GS1") & (labels == 1))[0]] = True
    with pytest.raises(ProxyStop, match="only 1 QC-passed cells"):
        compute_proxy_scores(counts[keep], donors[keep], labels[keep], gene_ids, _synthetic_panel(), n_control_panels=1)


def test_matched_control_panels_are_unique_direction_preserving_and_reproducible():
    counts, _donors, _labels, gene_ids = _synthetic_dataset()
    panel = _synthetic_panel()
    a = generate_matched_control_panels(counts, gene_ids, panel, n_panels=3, seed=55)
    b = generate_matched_control_panels(counts, gene_ids, panel, n_panels=3, seed=55)
    assert a == b
    for positive, negative in a:
        assert len(positive) == 48
        assert len(negative) == 73
        assert len(set(positive) | set(negative)) == 121
        assert not (set(positive) | set(negative)) & {row["gene_symbol"] for row in panel}


def test_freeze_manifest_hashes_all_declared_files_and_rejects_tampering(tmp_path):
    relative_paths = sorted(
        {
            "paper/trackA_tmz_response_121_gene_panel_2026-11-07.tsv",
            "paper/trackA_tmz_response_score_protocol_2026-11-07.md",
            "src/trackA_tmz_response_score.py",
            "tests/test_trackA_tmz_response_score.py",
            "src/189_trackA_yu2020_lodo.py",
            "paper/trackA_yu2020/PROTOCOL.md",
            "paper/trackA_yu2020/SHA256SUMS",
        }
    )
    lines = []
    for relative_path in relative_paths:
        file_path = tmp_path / relative_path
        file_path.parent.mkdir(parents=True, exist_ok=True)
        file_path.write_text(relative_path, encoding="utf-8")
        lines.append(f"{sha256_file(file_path)}  {relative_path}")
    manifest_path = tmp_path / "freeze.SHA256SUMS"
    manifest_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    verified = verify_freeze(root=tmp_path, freeze_path=manifest_path)
    assert set(verified) == set(relative_paths)

    (tmp_path / relative_paths[0]).write_text("changed", encoding="utf-8")
    with pytest.raises(ProxyStop, match="Frozen file hash mismatch"):
        verify_freeze(root=tmp_path, freeze_path=manifest_path)


def test_missing_freeze_stops_before_patient_loader_and_writes_no_result(tmp_path, monkeypatch):
    monkeypatch.setattr(proxy, "_load_frozen_a5_dataset", lambda: pytest.fail("patient loader called before freeze"))
    output_path = tmp_path / "result.json"
    with pytest.raises(ProxyStop, match="Proxy freeze manifest is missing"):
        proxy.run_proxy_analysis(output_path=output_path)
    assert not output_path.exists()
