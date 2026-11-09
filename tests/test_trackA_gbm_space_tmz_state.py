import json

import numpy as np
import pandas as pd
import pytest
import anndata as ad
from scipy import sparse

from src.trackA_gbm_space_tmz_state import (
    AnalysisStop,
    aggregate_backed_counts,
    classify_metadata,
    exact_sign_flip,
    generate_matched_controls,
    map_shared_panel,
    rank_percentiles,
    score_counts,
    sha256_file,
    verify_input_fingerprint,
    _validate_sparse_count_block,
)


def _synthetic_panel_and_symbols():
    omitted = ["AF131215.2", "PPIAP22", "CROCCP3", "HERC2P2", "HKR1", "LINCR-0001"]
    common = [f"G{i:03d}" for i in range(115)]
    panel_symbols = common + omitted
    rows = []
    for index, symbol in enumerate(panel_symbols):
        positive = index < 47 or symbol == "PPIAP22"
        if symbol == "AF131215.2" or symbol in {"CROCCP3", "HERC2P2", "HKR1", "LINCR-0001"}:
            positive = False
        rows.append({"gene_symbol": symbol, "rho_viability_100uM": 0.5 if positive else -0.5})
    gene_symbols = common
    gene_ids = [f"ENSG{i:06d}" for i in range(len(gene_symbols))]
    return rows, gene_symbols, gene_ids


def test_input_fingerprint_matches_and_fails_closed(tmp_path):
    input_path = tmp_path / "GBM_space_snRNA.h5ad"
    input_path.write_bytes(b"synthetic input bytes")
    fingerprint_path = tmp_path / "fingerprint.json"
    fingerprint_path.write_text(
        json.dumps({
            "filename": input_path.name,
            "bytes": input_path.stat().st_size,
            "sha256": sha256_file(input_path),
        }),
        encoding="utf-8",
    )
    assert verify_input_fingerprint(input_path, fingerprint_path) == sha256_file(input_path)
    input_path.write_bytes(b"x" * input_path.stat().st_size)
    with pytest.raises(AnalysisStop, match="SHA-256"):
        verify_input_fingerprint(input_path, fingerprint_path)


def test_maps_only_frozen_exact_115_gene_intersection():
    panel, symbols, ids = _synthetic_panel_and_symbols()
    mapped, index = map_shared_panel(panel, symbols, ids)
    assert len(mapped) == 115
    assert len(index) == 115
    assert sum(row["rho_viability_100uM"] > 0 for row in mapped) == 47
    assert sum(row["rho_viability_100uM"] < 0 for row in mapped) == 68


def test_panel_mapping_stops_if_missing_symbols_change():
    panel, symbols, ids = _synthetic_panel_and_symbols()
    symbols = symbols[:-1] + ["OTHER"]
    with pytest.raises(AnalysisStop, match="omissions differ"):
        map_shared_panel(panel, symbols, ids)


def test_classify_metadata_groups_only_frozen_malignant_states():
    expected = tuple(f"D{i}" for i in range(12))
    rows = []
    for donor in expected:
        rows.extend([
            {"donor_id": donor, "cell_status": "Malignant", "annotation_coarse": "OPC-like", "total_counts": 10},
            {"donor_id": donor, "cell_status": "Malignant", "annotation_coarse": "Hypoxic", "total_counts": 12},
            {"donor_id": donor, "cell_status": "Malignant", "annotation_coarse": "Proliferative", "total_counts": 8},
            {"donor_id": donor, "cell_status": "TME", "annotation_coarse": "Myeloid", "total_counts": 6},
        ])
    obs = pd.DataFrame(rows)
    codes, malignant, totals, counts = classify_metadata(obs, expected_donors=expected, minimum_cells=1)
    assert codes[0] == 0
    assert codes[1] == 1
    assert codes[2] == -1
    assert codes[3] == -1
    assert malignant.sum() == 36
    assert totals.shape == (48,)
    assert counts["D0"] == {"developmental_like": 1, "gliosis_hypoxia": 1}


def test_classify_metadata_stops_when_any_donor_group_is_below_minimum():
    expected = tuple(f"D{i}" for i in range(12))
    rows = []
    for donor in expected:
        rows.append({"donor_id": donor, "cell_status": "Malignant", "annotation_coarse": "OPC-like", "total_counts": 10})
        if donor != "D11":
            rows.append({"donor_id": donor, "cell_status": "Malignant", "annotation_coarse": "Hypoxic", "total_counts": 10})
    with pytest.raises(AnalysisStop, match="fewer than 1 gliosis/hypoxia"):
        classify_metadata(pd.DataFrame(rows), expected_donors=expected, minimum_cells=1)


def test_streaming_aggregation_builds_donor_state_pseudobulks_and_match_features():
    counts = sparse.csr_matrix(np.array([
        [1, 0, 2, 0],
        [0, 2, 1, 0],
        [3, 0, 0, 0],
        [0, 0, 1, 0],
    ], dtype=np.int64))
    group_codes = np.array([0, 0, 1, -1], dtype=np.int16)
    malignant = np.array([True, True, True, False])
    totals = np.asarray(counts.sum(axis=1)).reshape(-1)
    pseudobulk, means, detection, n_match = aggregate_backed_counts(
        counts, group_codes, malignant, totals, [0, 2], block_rows=2
    )
    assert pseudobulk.shape == (2, 4)
    assert pseudobulk[0].tolist() == [1, 2, 3, 0]
    assert pseudobulk[1].tolist() == [3, 0, 0, 0]
    assert n_match == 2
    assert means.shape == detection.shape == (4,)
    assert np.all((detection >= 0) & (detection <= 1))


def test_streaming_aggregation_accepts_backed_anndata_raw_matrix(tmp_path):
    counts = sparse.csr_matrix(np.array([[2, 0, 1], [0, 3, 0]], dtype=np.int64))
    obs = pd.DataFrame(index=["cell-a", "cell-b"])
    var = pd.DataFrame({"SYMBOL": ["A", "B", "C"]}, index=["ENSG1", "ENSG2", "ENSG3"])
    data = ad.AnnData(X=counts.copy(), obs=obs, var=var)
    data.raw = data
    h5ad_path = tmp_path / "mini.h5ad"
    data.write_h5ad(h5ad_path)
    backed = ad.read_h5ad(h5ad_path, backed="r")
    try:
        pseudobulk, _, _, n_match = aggregate_backed_counts(
            backed.raw.X,
            np.array([0, 1], dtype=np.int16),
            np.array([True, True]),
            np.array([3.0, 3.0]),
            [0, 1],
            block_rows=1,
        )
        assert pseudobulk.tolist() == [[2, 0, 1], [0, 3, 0]]
        assert n_match == 2
    finally:
        backed.file.close()


def test_explicit_sparse_zeros_do_not_inflate_detection_fraction():
    counts = sparse.csr_matrix(
        (
            np.array([2, 0, 3, 2, 0, 3], dtype=np.int64),
            np.array([0, 1, 2, 0, 1, 2], dtype=np.int32),
            np.array([0, 2, 3, 5, 6], dtype=np.int32),
        ),
        shape=(4, 3),
    )
    group_codes = np.array([0, 1, 0, 1], dtype=np.int16)
    malignant = np.ones(4, dtype=bool)
    totals = np.array([2, 3, 2, 3], dtype=np.float64)
    _, _, detection, n_match = aggregate_backed_counts(
        counts, group_codes, malignant, totals, [0, 1, 2, 3], block_rows=2
    )
    assert n_match == 4
    assert np.allclose(detection, [0.5, 0.0, 0.5])


def test_streaming_aggregation_stops_on_non_integer_or_bad_row_total():
    fractional = sparse.csr_matrix(np.array([[0.5, 0.5], [1.0, 1.0]], dtype=np.float64))
    with pytest.raises(AnalysisStop, match="non-integer"):
        aggregate_backed_counts(fractional, np.array([0, 1]), np.array([True, True]), np.array([1.0, 2.0]), [0, 1], block_rows=1)
    counts = sparse.csr_matrix(np.array([[1, 1], [1, 0]], dtype=np.int64))
    with pytest.raises(AnalysisStop, match="row sums"):
        aggregate_backed_counts(counts, np.array([0, 1]), np.array([True, True]), np.array([3.0, 1.0]), [0, 1], block_rows=1)


def test_count_validation_checks_stored_values_before_duplicate_sum():
    negative_duplicate = sparse.csr_matrix(
        (np.array([2.0, -1.0]), np.array([0, 0]), np.array([0, 2])), shape=(1, 1)
    )
    with pytest.raises(AnalysisStop, match="negative stored values"):
        _validate_sparse_count_block(negative_duplicate, np.array([1.0]))

    fractional_duplicates = sparse.csr_matrix(
        (np.array([0.5, 0.5]), np.array([0, 0]), np.array([0, 2])), shape=(1, 1)
    )
    with pytest.raises(AnalysisStop, match="non-integer stored values"):
        _validate_sparse_count_block(fractional_duplicates, np.array([1.0]))


def test_score_uses_signed_group_mean_rank_difference():
    counts = np.zeros(20, dtype=np.int64)
    counts[0:2] = [500, 450]
    counts[2:4] = [1, 1]
    score = score_counts(counts, positive_indices=[0, 1], negative_indices=[2, 3])
    assert score > 0
    ranks = rank_percentiles(np.ones(4))
    assert np.allclose(ranks, [0.5, 0.5, 0.5, 0.5])


def test_exact_sign_flip_reports_all_4096_assignments():
    result = exact_sign_flip(np.zeros(12, dtype=np.float64))
    assert result["assignments"] == 4096
    assert result["p_value"] == 1.0


def test_matched_control_panels_are_reproducible_and_unique():
    ids = [f"ENSG{i:03d}" for i in range(40)]
    means = np.linspace(0, 2, 40)
    detection = np.linspace(0, 1, 40) ** 2
    panel = [
        {"gene_symbol": "A", "rho_viability_100uM": 0.5},
        {"gene_symbol": "B", "rho_viability_100uM": 0.2},
        {"gene_symbol": "C", "rho_viability_100uM": -0.7},
    ]
    symbol_to_index = {"A": 0, "B": 1, "C": 2}
    first = generate_matched_controls(
        means, detection, ids, panel, symbol_to_index,
        excluded_gene_indices=[3], n_panels=10, seed=7, nearest_pool=4
    )
    second = generate_matched_controls(
        means, detection, ids, panel, symbol_to_index,
        excluded_gene_indices=[3], n_panels=10, seed=7, nearest_pool=4
    )
    assert first == second
    for positive, negative in first:
        assert len(positive) == 2
        assert len(negative) == 1
        assert len(set(positive + negative)) == 3
        assert not set(positive + negative).intersection({0, 1, 2, 3})
