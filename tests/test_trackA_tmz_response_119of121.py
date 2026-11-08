from __future__ import annotations

import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

import src.trackA_tmz_response_score as base
import src.trackA_tmz_response_score_119of121 as proxy119
from src.trackA_tmz_response_score import ProxyStop


def test_exact_matrix_mapping_keeps_119_source_genes_and_both_original_directions():
    panel = base.read_panel()
    present_ids = [row["gene_symbol"] for row in panel if row["gene_symbol"] not in proxy119.EXPECTED_UNMAPPED_SYMBOLS]
    present_ids.extend(("PPIA3L", "ENSG00000255310", "BACKGROUND1"))

    mapped, missing = proxy119.map_panel_to_matrix(panel, present_ids)

    assert missing == proxy119.EXPECTED_UNMAPPED_SYMBOLS
    assert len(mapped) == 119
    assert {row["gene_symbol"] for row in mapped} == set(present_ids) - {
        "PPIA3L", "ENSG00000255310", "BACKGROUND1"
    }
    assert sum(float(row["rho_viability_100uM"]) > 0 for row in mapped) == 47
    assert sum(float(row["rho_viability_100uM"]) < 0 for row in mapped) == 72
    assert all(row["direction"] == (
        "resistance_associated" if float(row["rho_viability_100uM"]) > 0 else "sensitivity_associated"
    ) for row in mapped)


def test_mapping_stops_if_any_other_source_symbol_is_missing():
    panel = base.read_panel()
    extra_missing = next(
        row["gene_symbol"] for row in panel
        if row["gene_symbol"] not in proxy119.EXPECTED_UNMAPPED_SYMBOLS
    )
    matrix_ids = [
        row["gene_symbol"] for row in panel
        if row["gene_symbol"] not in (*proxy119.EXPECTED_UNMAPPED_SYMBOLS, extra_missing)
    ]

    with pytest.raises(ProxyStop, match="missing symbols differ from the frozen amendment"):
        proxy119.map_panel_to_matrix(panel, matrix_ids)


def test_mapping_stops_on_duplicate_matrix_gene_ids():
    panel = base.read_panel()
    matrix_ids = [row["gene_symbol"] for row in panel if row["gene_symbol"] not in proxy119.EXPECTED_UNMAPPED_SYMBOLS]
    matrix_ids.append(matrix_ids[0])

    with pytest.raises(ProxyStop, match="duplicate gene identifiers"):
        proxy119.map_panel_to_matrix(panel, matrix_ids)


def test_119_gene_panel_scores_synthetic_donors_and_builds_direction_matched_controls():
    panel = base.read_panel()
    gene_ids = [row["gene_symbol"] for row in panel if row["gene_symbol"] not in proxy119.EXPECTED_UNMAPPED_SYMBOLS]
    gene_ids.extend(f"BG{i:03d}" for i in range(300))
    mapped, _missing = proxy119.map_panel_to_matrix(panel, gene_ids)
    rng = np.random.default_rng(119121)
    rows = []
    donors = []
    labels = []
    positive = {row["gene_symbol"] for row in mapped if float(row["rho_viability_100uM"]) > 0}
    negative = {row["gene_symbol"] for row in mapped if float(row["rho_viability_100uM"]) < 0}
    index = {gene: i for i, gene in enumerate(gene_ids)}
    for donor in base.PRIMARY_DONORS:
        for label in (1, 0):
            for _ in range(20):
                values = rng.integers(2, 10, size=len(gene_ids), dtype=np.int32)
                up, down = (positive, negative) if label == 1 else (negative, positive)
                values[[index[gene] for gene in up]] += 20
                values[[index[gene] for gene in down]] = 1
                rows.append(values)
                donors.append(donor)
                labels.append(label)

    result = base.compute_proxy_scores(
        np.asarray(rows), donors, labels, gene_ids, mapped, n_control_panels=20, control_seed=119121
    )

    assert result["status"] == "complete_exploratory"
    assert result["panel"]["n_genes"] == 119
    assert result["panel"]["positive_resistance_associated"] == 47
    assert result["panel"]["negative_sensitivity_associated"] == 72
    assert result["cohort"]["n_donors"] == 6
    assert result["matched_controls"]["n_panels"] == 20


def test_amendment_freeze_verifies_exact_files_and_detects_tampering(tmp_path):
    lines = []
    for relative in proxy119.AMENDMENT_SEAL_FILES:
        path = tmp_path / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(relative, encoding="utf-8")
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        lines.append(f"{digest}  {relative}")
    manifest = tmp_path / "amendment.SHA256SUMS"
    manifest.write_text("\n".join(lines) + "\n", encoding="utf-8")

    verified = proxy119.verify_amendment_freeze(root=tmp_path, freeze_path=manifest)
    assert set(verified) == set(proxy119.AMENDMENT_SEAL_FILES)

    (tmp_path / proxy119.AMENDMENT_SEAL_FILES[0]).write_text("tampered", encoding="utf-8")
    with pytest.raises(ProxyStop, match="Frozen file hash mismatch"):
        proxy119.verify_amendment_freeze(root=tmp_path, freeze_path=manifest)


def test_project_base_and_amendment_freezes_match_all_declared_files():
    assert set(base.verify_freeze()) == {
        "paper/trackA_tmz_response_121_gene_panel_2026-11-07.tsv",
        "paper/trackA_tmz_response_score_protocol_2026-11-07.md",
        "src/trackA_tmz_response_score.py",
        "tests/test_trackA_tmz_response_score.py",
        "src/189_trackA_yu2020_lodo.py",
        "paper/trackA_yu2020/PROTOCOL.md",
        "paper/trackA_yu2020/SHA256SUMS",
    }
    assert set(proxy119.verify_amendment_freeze()) == set(proxy119.AMENDMENT_SEAL_FILES)


def test_missing_amendment_freeze_stops_before_loading_a5(tmp_path, monkeypatch):
    monkeypatch.setattr(proxy119, "AMENDMENT_FREEZE_PATH", tmp_path / "missing.SHA256SUMS")
    monkeypatch.setattr(base, "_load_frozen_a5_dataset", lambda: pytest.fail("A5 data loaded before freeze"))
    output = tmp_path / "result.json"

    with pytest.raises(ProxyStop, match="Amendment freeze manifest is missing"):
        proxy119.run_proxy_analysis(output_path=output)
    assert not output.exists()


def test_run_records_119_of_121_coverage_and_never_overwrites(tmp_path, monkeypatch):
    panel = base.read_panel()
    gene_ids = [row["gene_symbol"] for row in panel if row["gene_symbol"] not in proxy119.EXPECTED_UNMAPPED_SYMBOLS]
    gene_ids.extend(f"BG{i}" for i in range(130))
    loaded_donors = np.asarray(
        [donor for donor in (*base.PRIMARY_DONORS, *proxy119.EXPECTED_LENIENT_LOADER_DONORS) for _ in (0, 1)]
    )
    dataset = SimpleNamespace(
        counts=np.zeros((len(loaded_donors), len(gene_ids))),
        donors=loaded_donors,
        labels=np.tile(np.asarray((1, 0), dtype=np.int64), len(loaded_donors) // 2),
        gene_ids=gene_ids,
    )
    monkeypatch.setattr(base, "verify_freeze", lambda: {"baseline": "hash"})
    monkeypatch.setattr(proxy119, "verify_amendment_freeze", lambda: {"amendment": "hash"})
    monkeypatch.setattr(base, "_load_frozen_a5_dataset", lambda: (dataset, base.PRIMARY_DONORS, {}, {"A5": "hash"}))
    amendment = tmp_path / "amendment.md"
    amendment.write_text("frozen amendment", encoding="utf-8")
    amendment_freeze = tmp_path / "amendment.SHA256SUMS"
    amendment_freeze.write_text("frozen manifest", encoding="utf-8")
    monkeypatch.setattr(proxy119, "AMENDMENT_PATH", amendment)
    monkeypatch.setattr(proxy119, "AMENDMENT_FREEZE_PATH", amendment_freeze)
    observed = {}

    def fake_compute(_counts, scoring_donors, _labels, _gene_ids, mapped_panel):
        observed["panel"] = mapped_panel
        observed["donors"] = set(np.asarray(scoring_donors, dtype=str).tolist())
        observed["cell_count"] = len(scoring_donors)
        return {"status": "complete_exploratory", "score": {}, "matched_controls": {}}

    monkeypatch.setattr(base, "compute_proxy_scores", fake_compute)
    output = tmp_path / "result.json"
    result = proxy119.run_proxy_analysis(output_path=output)

    assert result["panel"]["n_source_genes"] == 121
    assert result["panel"]["n_scored_genes"] == 119
    assert result["panel"]["coverage"] == "119/121 exact matrix symbols"
    assert result["panel"]["unmapped_source_symbols"] == list(proxy119.EXPECTED_UNMAPPED_SYMBOLS)
    assert len(observed["panel"]) == 119
    assert observed["donors"] == set(base.PRIMARY_DONORS)
    assert observed["cell_count"] == len(base.PRIMARY_DONORS) * 2
    assert result["provenance"]["a5_cell_flow"]["loaded_donors"] == sorted(
        {*base.PRIMARY_DONORS, *proxy119.EXPECTED_LENIENT_LOADER_DONORS}
    )
    assert result["provenance"]["a5_cell_flow"]["scored_donors"] == list(base.PRIMARY_DONORS)
    assert result["provenance"]["a5_cell_flow"]["scored_cells"] == len(base.PRIMARY_DONORS) * 2
    assert json.loads(output.read_text(encoding="utf-8"))["panel"] == result["panel"]

    with pytest.raises(ProxyStop, match="Refusing to overwrite existing result"):
        proxy119.run_proxy_analysis(output_path=output)


def test_unexpected_loader_donor_stops_before_score(tmp_path, monkeypatch):
    panel = base.read_panel()
    gene_ids = [row["gene_symbol"] for row in panel if row["gene_symbol"] not in proxy119.EXPECTED_UNMAPPED_SYMBOLS]
    loaded_donors = np.asarray((*base.PRIMARY_DONORS, "GS99"), dtype=str)
    dataset = SimpleNamespace(
        counts=np.zeros((len(loaded_donors), len(gene_ids))),
        donors=loaded_donors,
        labels=np.zeros(len(loaded_donors), dtype=np.int64),
        gene_ids=gene_ids,
    )
    monkeypatch.setattr(base, "verify_freeze", lambda: {"baseline": "hash"})
    monkeypatch.setattr(proxy119, "verify_amendment_freeze", lambda: {"amendment": "hash"})
    monkeypatch.setattr(base, "_load_frozen_a5_dataset", lambda: (dataset, base.PRIMARY_DONORS, {}, {"A5": "hash"}))
    monkeypatch.setattr(base, "compute_proxy_scores", lambda *_args: pytest.fail("score computed for unexpected donor"))
    output = tmp_path / "result.json"

    with pytest.raises(ProxyStop, match="Unexpected donor IDs from A5 loader"):
        proxy119.run_proxy_analysis(output_path=output)
    assert not output.exists()
