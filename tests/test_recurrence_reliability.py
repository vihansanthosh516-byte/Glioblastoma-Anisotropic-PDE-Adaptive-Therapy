"""Synthetic checks for src/recurrence_reliability (no real data needed)."""
import gzip
import sys
from pathlib import Path

import numpy as np
import pytest
import scipy.sparse as sp
from scipy.stats import rankdata

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from recurrence_reliability import inference as inf  # noqa: E402
from recurrence_reliability import metadata as md  # noqa: E402
from recurrence_reliability import scoring as sc  # noqa: E402


def test_barcode_suffixes_normalise_to_same_key():
    assert md.norm_barcode("AAACCCACAGAGTGAC-1") == md.norm_barcode("AAACCCACAGAGTGAC.1") == "AAACCCACAGAGTGAC"
    assert md.norm_barcode("AAACCCACAGAGTGAC") == "AAACCCACAGAGTGAC"


def test_annotation_match_requires_unique_covering_key():
    lib = {f"B{i}" for i in range(100)}
    ann = {"SF1": {f"B{i}": "Tumor" for i in range(100)}, "SF1v2": {f"B{i}": "Normal" for i in range(200, 260)}}
    assert md.match_annotation(lib, ann, "SF1")[0] == "SF1"
    ann_v2 = {"SF2v2": {f"B{i}": "Tumor" for i in range(100)}}
    assert md.match_annotation(lib, ann_v2, "SF2")[0] == "SF2v2"
    both = {"SF3": {f"B{i}": "Tumor" for i in range(100)}, "SF3v2": {f"B{i}": "Tumor" for i in range(100)}}
    assert md.match_annotation(lib, both, "SF3")[0] is None
    low = {"SF4": {f"B{i}": "Tumor" for i in range(80)}}   # 80% join < 95%
    assert md.match_annotation(lib, low, "SF4")[0] is None
    extra = {"SF5": {**{f"B{i}": "Tumor" for i in range(100)}, "ZZ": "Normal"}}   # label not in library
    assert md.match_annotation(lib, extra, "SF5")[0] is None


def test_conflicting_duplicate_annotation_raises(tmp_path):
    p = tmp_path / "a.txt.gz"
    with gzip.open(p, "wt") as f:
        f.write("Sample#_Barcode Tumor_Normal_annotation\nSF1_AAA Tumor\nSF1_AAA Tumor\nSF1_CCC Normal\n")
    assert md.load_annotation(p) == {"SF1": {"AAA": "Tumor", "CCC": "Normal"}}
    with gzip.open(p, "wt") as f:
        f.write("Sample#_Barcode Tumor_Normal_annotation\nSF1_AAA Tumor\nSF1_AAA Normal\n")
    with pytest.raises(ValueError):
        md.load_annotation(p)


def test_library_relationship_classes():
    a = {f"B{i}" for i in range(1000)}
    assert md.library_relationship(a, {f"B{i}" for i in range(100)})[0] == "same_capture_barcode_subset"
    assert md.library_relationship(a, {f"C{i}" for i in range(500)} | {"B1"})[0] == "independent_capture"
    assert md.library_relationship(a, {f"B{i}" for i in range(50)} | {f"C{i}" for i in range(50)})[0] == "ambiguous_partial_overlap"


def test_log_cpm_hand_values():
    m = sp.csc_matrix(np.array([[1, 0], [3, 10]]))
    le = sc.log_cpm(m).toarray()
    assert le[0, 0] == pytest.approx(np.log2(0.25e6 / 10 + 1))
    assert le[1, 1] == pytest.approx(np.log2(1e6 / 10 + 1))
    assert le[0, 1] == 0


def test_tirosh_score_hand_value_and_controls_exclude_module():
    le = sp.csr_matrix(np.array([[2.0, 0.0], [4.0, 2.0], [1.0, 1.0], [3.0, 3.0]]))
    center = np.array([1.0, 3.0, 1.0, 3.0])
    s = sc.tirosh_score(le, center, np.array([0, 1]), np.array([2, 3]))
    # cell 0: module mean (2+4)/2 - 2 = 1; control mean (1+3)/2 - 2 = 0
    assert s[0] == pytest.approx(1.0)
    assert s[1] == pytest.approx((1.0 - 2.0) - (2.0 - 2.0))
    bins = np.array([0, 0, 0, 1, 1, 1])
    ctrl = sc.control_genes(np.array([0, 3]), bins, np.random.default_rng(0), n_ctrl=2)
    assert not ({0, 3} & set(ctrl.tolist()))


def test_specimen_difference_does_not_depend_on_centering_for_one_module():
    rng = np.random.default_rng(1)
    le = sp.csr_matrix(rng.gamma(1, 1, size=(40, 30)))
    mod, ctl = np.arange(5), np.arange(10, 30)
    a = sc.tirosh_score(le, np.zeros(40), mod, ctl)
    b = sc.tirosh_score(le, rng.normal(size=40), mod, ctl)
    assert (a[:15].mean() - a[15:].mean()) == pytest.approx(b[:15].mean() - b[15:].mean())


def ucell_bruteforce(dense, module, max_rank):
    out = []
    for j in range(dense.shape[1]):
        r = rankdata(-dense[:, j], method="average")
        r = np.minimum(r[module], max_rank + 1)
        n = len(module)
        out.append(1 - (r.sum() - n * (n + 1) / 2) / (n * max_rank))
    return np.array(out)


@pytest.mark.parametrize("max_rank", [5, 20, 1000])
def test_ucell_matches_bruteforce_ranks(max_rank):
    rng = np.random.default_rng(2)
    dense = rng.poisson(0.7, size=(60, 25)).astype(float)
    dense[:, 0] = 0
    dense[0, 0] = 1                     # at least one count per nucleus
    module = np.array([0, 3, 7, 11, 40])
    got = sc.ucell_score(sp.csc_matrix(dense), module, max_rank=max_rank)
    assert np.allclose(got, ucell_bruteforce(dense, module, max_rank))


def test_scores_follow_cell_reordering():
    rng = np.random.default_rng(3)
    dense = rng.poisson(1.0, size=(50, 20)).astype(float) + np.eye(50, 20)
    perm = rng.permutation(20)
    m, mp = sp.csc_matrix(dense), sp.csc_matrix(dense[:, perm])
    mod = np.array([1, 2, 3])
    assert np.allclose(sc.ucell_score(m, mod)[perm], sc.ucell_score(mp, mod))
    le, lep = sc.log_cpm(m), sc.log_cpm(mp)
    t, tp = sc.tirosh_score(le, np.zeros(50), mod, np.arange(10, 30)), sc.tirosh_score(lep, np.zeros(50), mod, np.arange(10, 30))
    assert np.allclose(t[perm], tp)
    assert t.mean() == pytest.approx(tp.mean())


def test_higher_module_expression_gives_higher_scores():
    base = np.ones((100, 10))
    up = base.copy()
    up[:5, 5:] = 20
    m = sp.csc_matrix(up)
    mod = np.arange(5)
    u = sc.ucell_score(m, mod, max_rank=50)
    t = sc.tirosh_score(sc.log_cpm(m), np.zeros(100), mod, np.arange(20, 100))
    assert u[5:].mean() > u[:5].mean() and t[5:].mean() > t[:5].mean()


def test_direction_sign_and_ties():
    se = 1.0
    assert inf.direction(2.5, se) == 1          # recurrent - primary > 0 => increase
    assert inf.direction(-2.5, se) == -1
    assert inf.direction(inf.Z * se, se) == 0   # boundary is indeterminate
    assert inf.direction(0.0, se) == 0
    with pytest.raises(ValueError):
        inf.direction(float("nan"), se)


def test_agreement_table_counts_each_patient_once():
    r = [1, 1, -1, 0, 0, 1]
    a = [1, -1, -1, 0, 1, 0]
    t = inf.agreement_table(r, a)
    assert t["n"] == 6
    assert t["opposite_direction"] == 1
    assert t["one_indeterminate"] == 2
    assert t["both_indeterminate"] == 1
    assert t["same_determinate_direction"] == 2
    assert sum(t["table_R_by_A"].values()) == 6
    assert t["any_disagreement"] == 3
    with pytest.raises(ValueError):
        inf.agreement_table([1], [1, 1])


def test_exact_bounds():
    assert inf.one_sided_upper_zero(27) == pytest.approx(0.105019, abs=1e-6)
    lo, hi = inf.clopper_pearson(0, 23)
    assert lo == 0 and hi == pytest.approx(1 - 0.025 ** (1 / 23))
    lo, hi = inf.clopper_pearson(5, 10)
    assert lo == pytest.approx(0.187086, abs=1e-5) and hi == pytest.approx(0.812914, abs=1e-5)


def test_cell_bootstrap_se_close_to_analytic():
    v = np.random.default_rng(4).normal(size=400)
    se = inf.boot_se(v, np.random.default_rng(5), 2000)
    assert se == pytest.approx(v.std(ddof=1) / np.sqrt(len(v)), rel=0.1)


def test_patient_bootstrap_mean_hand_case():
    out = inf.patient_bootstrap_mean(np.array([1.0, 1.0, 1.0]), np.random.default_rng(0), 100)
    assert out["mean"] == 1.0 and out["ci95_percentile"] == [1.0, 1.0] and out["n"] == 3


class FakeSpecimen:
    def __init__(self, counts, tumor, symbols):
        self.counts, self.tumor, self.symbols = sp.csc_matrix(counts), np.asarray(tumor), symbols

    def subset(self, malignant_only):
        return self.counts[:, np.flatnonzero(self.tumor)] if malignant_only else self.counts


def make_fake_cohort(n_pat=8, n_genes=400, n_cells=120, seed=0):
    rng = np.random.default_rng(seed)
    symbols = np.array([f"G{i}" for i in range(n_genes)])
    mes1, mes2, verh = list(symbols[:20]), list(symbols[20:40]), list(symbols[40:60])
    base = rng.gamma(0.5, 1.0, size=n_genes)
    specs, cohort, planted = {}, [], {}
    for p in range(n_pat):
        effect = [3.0, -0.7, 0.0][p % 3]          # up, down, none (log-scale multiplier on MES1)
        planted[p] = np.sign(effect)
        for stage, mult in (("Primary", 0.0), ("Recurrent", effect)):
            lam = np.tile(base[:, None], (1, n_cells)) * rng.uniform(0.5, 1.5, size=(1, n_cells))
            lam[:20] *= np.exp(mult)
            sid = f"S{p}{stage[0]}"
            specs[sid] = FakeSpecimen(rng.poisson(lam) + (np.arange(n_genes)[:, None] == n_genes - 1),
                                      rng.uniform(size=n_cells) < 0.8, symbols)
        cohort.append({"table1_pair": p, "Primary": f"S{p}P", "Recurrent": f"S{p}R"})
    cfg = {"seeds": {"control_genes": 1, "cell_bootstrap": 2, "patient_bootstrap": 3, "matched_depth": 4, "split_half": 5},
           "n_cell_bootstrap": 200, "n_patient_bootstrap": 500, "matched_cell_number": 50, "matched_reps": 20,
           "split_half_reps": 5, "split_half_boot": 100, "tolerance_opposite_fraction": 0.05,
           "independent_capture_libraries": {"S0R": "S0R_b2"}}
    b2 = {"S0R": specs["S0R"]}
    return specs, b2, cohort, {"MESlike1": mes1, "MESlike2": mes2}, verh, cfg, planted


def test_end_to_end_synthetic_recovers_planted_directions():
    from recurrence_reliability.run import analyze
    specs, b2, cohort, neftel, verh, cfg, planted = make_fake_cohort()
    rows, res = analyze(specs, b2, cohort, neftel, verh, cfg)
    assert [r["table1_pair"] for r in rows] == list(range(8))
    assert len({r["primary"] for r in rows} | {r["recurrent"] for r in rows}) == 16   # no specimen twice
    for r in rows:
        if planted[r["table1_pair"]] > 0:
            assert r["dir_R"] == 1 and r["dir_A"] == 1 and r["delta_R"] > 0
        assert np.isfinite(r["delta_R"]) and np.isfinite(r["delta_A"]) and r["se_R"] > 0
    pe = res["primary_endpoint"]
    assert pe["n"] == 8 and sum(pe["table_R_by_A"].values()) == 8
    assert 0 <= pe["opposite_ci95_clopper_pearson"][0] <= pe["opposite_ci95_clopper_pearson"][1] <= 1
    assert res["secondary"]["S6_split_half_false_determinate_rate"]["R"] < 0.25
    s3 = res["secondary"]["S3_library_aware"][0]
    assert s3["R_pooled_libraries"]["dir"] == s3["R_primary_library"]["dir"]


def test_analysis_is_deterministic_given_seeds():
    from recurrence_reliability.run import analyze
    a = analyze(*make_fake_cohort()[:6])
    b = analyze(*make_fake_cohort()[:6])
    assert [r["delta_R"] for r in a[0]] == [r["delta_R"] for r in b[0]]
    assert a[1]["primary_endpoint"] == b[1]["primary_endpoint"]
