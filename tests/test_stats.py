"""Known-answer tests for the patient-level statistics (sctcna.stats)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from statsmodels.stats.multitest import multipletests

from sctcna.stats import bh_adjust, compartment_medians, paired_noise_tests, paired_score_test, rank_biserial

PATIENTS = [f"P{i}" for i in range(1, 7)]
MIN_P_SIX_PAIRS = 2 / 64  # smallest two-sided exact Wilcoxon p for six non-zero paired differences
MIN_P_FIVE_PAIRS = 2 / 32


def test_bh_adjust_matches_statsmodels_and_keeps_nan():
    rng = np.random.default_rng(0)
    p = rng.uniform(size=200) ** 3
    p[[5, 17, 99]] = np.nan
    q = bh_adjust(p)
    ok = np.isfinite(p)
    assert np.isnan(q[~ok]).all()
    np.testing.assert_allclose(q[ok], multipletests(p[ok], method="fdr_bh")[1], rtol=0, atol=1e-15)
    # hand-computed: sorted p 0.01, 0.03, 0.04 -> 0.03, 0.045, 0.04 -> monotone from the top: 0.03, 0.04, 0.04
    np.testing.assert_allclose(bh_adjust([0.01, 0.04, 0.03]), [0.03, 0.04, 0.04], rtol=0, atol=1e-15)
    assert np.isnan(bh_adjust([np.nan, np.nan])).all()
    assert bh_adjust([]).shape == (0,)


def test_rank_biserial_known_values():
    assert rank_biserial([1.0, 2.0, 3.0]) == 1.0
    assert rank_biserial([-0.5, -2.0]) == -1.0
    assert rank_biserial([1.0, -1.0, 2.0, -2.0]) == 0.0
    assert rank_biserial([0.0, 0.0, 1.0]) == 1.0  # zero differences are dropped
    assert rank_biserial([3.0, -1.0, 2.0]) == pytest.approx((3 + 2 - 1) / 6)
    assert np.isnan(rank_biserial([0.0, 0.0]))
    assert np.isnan(rank_biserial([]))
    assert np.isnan(rank_biserial([np.nan]))


def _long_table(dataset: str = "D") -> pd.DataFrame:
    """Per-patient CV for six patients.

    DOUBLE: CV doubled in malignant cells for every patient; SAME: identical CVs; HALF: halved;
    DROP: doubled, but patient P1 has CV 0 in non-malignant cells (pair discarded); FEW: only three patients.
    """
    rows = []
    for i, patient in enumerate(PATIENTS):
        base = 0.5 + 0.1 * i
        genes = {"DOUBLE": (2 * base, base), "SAME": (base, base), "HALF": (base, 2 * base), "DROP": (2 * base, 0.0 if i == 0 else base)}
        if i < 3:
            genes["FEW"] = (3.0, 1.0)
        for gene, (malignant, non_malignant) in genes.items():
            rows.append({"dataset": dataset, "gene": gene, "patient": patient, "compartment": "malignant", "cv": malignant, "mean": 10.0 + i})
            rows.append({"dataset": dataset, "gene": gene, "patient": patient, "compartment": "non_malignant", "cv": non_malignant, "mean": 20.0 + i})
    return pd.DataFrame(rows)


def test_paired_noise_tests_known_answers():
    res = paired_noise_tests(_long_table(), metric="cv", min_patients=5).set_index("gene")
    assert set(res.columns) == {"dataset", "n_pairs", "log2fc_cv_median", "log2fc_cv_mean", "wilcoxon_stat", "p_value", "rank_biserial", "q_value"}
    assert res["n_pairs"].to_dict() == {"DOUBLE": 6, "SAME": 6, "HALF": 6, "DROP": 5, "FEW": 3}
    assert res.loc["DOUBLE", "log2fc_cv_median"] == 1.0
    assert res.loc["DOUBLE", "log2fc_cv_mean"] == 1.0
    assert res.loc["DOUBLE", "rank_biserial"] == 1.0
    assert res.loc["DOUBLE", "wilcoxon_stat"] == 0.0
    assert res.loc["DOUBLE", "p_value"] == MIN_P_SIX_PAIRS
    assert res.loc["HALF", "log2fc_cv_median"] == -1.0
    assert res.loc["HALF", "rank_biserial"] == -1.0
    assert res.loc["HALF", "p_value"] == MIN_P_SIX_PAIRS
    assert res.loc["SAME", "log2fc_cv_median"] == 0.0
    assert res.loc["SAME", "p_value"] == 1.0
    assert np.isnan(res.loc["SAME", "rank_biserial"])
    assert res.loc["DROP", "log2fc_cv_median"] == 1.0
    assert res.loc["DROP", "p_value"] == MIN_P_FIVE_PAIRS
    assert res.loc["FEW", "log2fc_cv_median"] == pytest.approx(np.log2(3.0))
    assert np.isnan(res.loc["FEW", ["wilcoxon_stat", "p_value", "rank_biserial", "q_value"]].to_numpy(dtype=float)).all()
    # BH over the four tested genes: p = 1/32, 1/32, 1/16, 1 -> q = 1/16, 1/16, 1/12, 1
    expected_q = {"DOUBLE": 0.0625, "HALF": 0.0625, "DROP": 0.0625 * 4 / 3, "SAME": 1.0}
    for gene, q in expected_q.items():
        assert res.loc[gene, "q_value"] == pytest.approx(q, abs=1e-15)


def test_q_values_are_adjusted_within_dataset():
    other = _long_table("E")
    other = other[other["gene"] == "DOUBLE"]
    res = paired_noise_tests(pd.concat([_long_table("D"), other], ignore_index=True), "cv", min_patients=5).set_index(["dataset", "gene"])
    assert res.loc[("E", "DOUBLE"), "q_value"] == MIN_P_SIX_PAIRS  # the only test in E
    assert res.loc[("D", "DOUBLE"), "q_value"] == pytest.approx(0.0625, abs=1e-15)
    for dataset in ("D", "E"):
        sub = res.loc[dataset]
        np.testing.assert_allclose(sub["q_value"].to_numpy(dtype=float), bh_adjust(sub["p_value"].to_numpy(dtype=float)), rtol=0, atol=1e-15)


def test_paired_noise_tests_without_both_compartments_is_empty():
    long = _long_table()
    res = paired_noise_tests(long[long["compartment"] == "malignant"], "cv")
    assert res.empty
    assert "p_value" in res.columns and "q_value" in res.columns


def test_compartment_medians_layout_and_values():
    wide = compartment_medians(_long_table())
    assert set(wide.columns) == {"dataset", "gene", "mean_malignant", "mean_non_malignant", "cv_malignant", "cv_non_malignant", "n_patients_malignant", "n_patients_non_malignant"}
    wide = wide.set_index("gene")
    assert wide.loc["DOUBLE", "cv_malignant"] == 1.5  # median of 1.0, 1.2, ..., 2.0
    assert wide.loc["DOUBLE", "cv_non_malignant"] == 0.75
    assert wide.loc["DOUBLE", "mean_malignant"] == 12.5
    assert wide.loc["DOUBLE", "mean_non_malignant"] == 22.5
    assert wide.loc["DOUBLE", "n_patients_malignant"] == 6
    assert wide.loc["FEW", "n_patients_malignant"] == 3
    assert wide.loc["FEW", "cv_malignant"] == 3.0
    assert (wide["dataset"] == "D").all()
    only_cv = compartment_medians(_long_table(), metrics=["cv"])
    assert set(only_cv.columns) == {"dataset", "gene", "cv_malignant", "cv_non_malignant", "n_patients_malignant", "n_patients_non_malignant"}


def test_paired_score_test_shape_and_values():
    scores = pd.DataFrame(
        {
            "patient": PATIENTS * 2,
            "compartment": ["malignant"] * 6 + ["non_malignant"] * 6,
            "depth2_pseudobulk": [1.0, 2.0, 3.0, 4.0, 5.0, 6.0] + [1.0] * 6,
        }
    )
    res = paired_score_test(scores, "depth2_pseudobulk")
    assert set(res) == {"n_pairs", "median_diff", "p_value", "rank_biserial"}
    assert res["n_pairs"] == 6
    assert res["median_diff"] == 2.5  # differences 0, 1, 2, 3, 4, 5
    assert res["p_value"] == MIN_P_FIVE_PAIRS  # five non-zero differences
    assert res["rank_biserial"] == 1.0
    few = paired_score_test(scores, "depth2_pseudobulk", min_patients=7)
    assert few["n_pairs"] == 6 and np.isnan(few["p_value"]) and np.isnan(few["rank_biserial"])
    missing = paired_score_test(scores[scores["compartment"] == "malignant"], "depth2_pseudobulk")
    assert missing["n_pairs"] == 0
    assert np.isnan(missing["median_diff"]) and np.isnan(missing["p_value"]) and np.isnan(missing["rank_biserial"])


def test_one_sample_wilcoxon():
    from sctcna.stats import one_sample_wilcoxon

    res = one_sample_wilcoxon([0.5, 1.0, 0.2, 0.3, 0.8, 0.1], min_n=5)
    assert res["n"] == 6 and res["frac_positive"] == 1.0
    assert res["p_value"] == pytest.approx(2 / 64)
    assert res["rank_biserial"] == 1.0
    few = one_sample_wilcoxon([1.0, np.nan, 2.0], min_n=5)
    assert few["n"] == 2 and np.isnan(few["p_value"])
