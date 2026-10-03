"""Tests for the four-term decomposition reference implementation.

Two layers:

1. Analytic toy tests on tiny hand-checkable datasets, where each term can be
   worked out by hand and the identity must close exactly.
2. An archived regression test against the values reported in the accompanying
   case study, read from the archived patient-level table.

Runs with either `python3 -m unittest discover -s tests` or `pytest tests`.
"""
import math
import os
import sys
import unittest

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "code"))
from composition_decomposition import decompose, _population_cov  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ARCHIVE = os.environ.get(
    "ICOS_PATIENT_STATE_TABLE",
    os.path.join(HERE, "..", "data", "figure2_patient_state_table.csv"),
)

# Values reported in the case study, from analysis_results/figure2_formal_decomposition.csv
REPORTED = {
    "within_state": 0.0036532309026953433,
    "composition": -0.0028159006193097234,
    "interaction": 6.520557245153435e-05,
    "covariance": -0.0004796651992069352,
    "reconstructed_total": 0.0004228706566302191,
    "observed_total": 0.0004228706566302212,
}
REPORTED_COV_BY_STATE = {
    "CM": -0.000576831461694871,
    "Eff/EM": 9.716626248793583e-05,
}


def frame(rows):
    return pd.DataFrame(rows, columns=["patient_id", "group", "state",
                                       "state_fraction", "state_mean"])


class TestAnalyticToy(unittest.TestCase):

    def test_pure_composition_shift(self):
        """State means identical in both groups; only the mixture moves.

        With mu = (0, 1) in both groups and the state-B fraction going from 0.25
        to 0.75, the whole contrast must sit in the composition term:
        within = 0, interaction = 0, covariance = 0 (no within-group spread),
        composition = 0*(0.75-0.25) + 1*(0.75-0.25) = 0.5.
        """
        rows = []
        for p in ("a", "b"):
            rows += [(p, "g0", "A", 0.75, 0.0), (p, "g0", "B", 0.25, 1.0)]
        for p in ("c", "d"):
            rows += [(p, "g1", "A", 0.25, 0.0), (p, "g1", "B", 0.75, 1.0)]
        r = decompose(frame(rows), "g0", "g1")
        self.assertAlmostEqual(r["within_state"], 0.0, places=15)
        self.assertAlmostEqual(r["interaction"], 0.0, places=15)
        self.assertAlmostEqual(r["covariance"], 0.0, places=15)
        self.assertAlmostEqual(r["composition"], 0.5, places=15)
        self.assertAlmostEqual(r["observed_total"], 0.5, places=15)
        self.assertLess(abs(r["closure_error"]), 1e-15)

    def test_pure_within_state_shift(self):
        """Mixture identical in both groups; only the state means move.

        Weights fixed at (0.5, 0.5); means go from (0, 0) to (1, 3).
        within = 0.5*1 + 0.5*3 = 2; composition = interaction = covariance = 0.
        """
        rows = []
        for p in ("a", "b"):
            rows += [(p, "g0", "A", 0.5, 0.0), (p, "g0", "B", 0.5, 0.0)]
        for p in ("c", "d"):
            rows += [(p, "g1", "A", 0.5, 1.0), (p, "g1", "B", 0.5, 3.0)]
        r = decompose(frame(rows), "g0", "g1")
        self.assertAlmostEqual(r["within_state"], 2.0, places=15)
        self.assertAlmostEqual(r["composition"], 0.0, places=15)
        self.assertAlmostEqual(r["interaction"], 0.0, places=15)
        self.assertAlmostEqual(r["covariance"], 0.0, places=15)
        self.assertLess(abs(r["closure_error"]), 1e-15)

    def test_covariance_term_is_not_a_residual(self):
        """A case built so the covariance term alone is non-zero.

        Both groups have the same mean weights (0.5, 0.5) and the same mean state
        means, so within, composition and interaction all vanish. In g1, weight
        and mean are correlated across patients; in g0 they are not. The whole
        contrast is therefore the covariance term, and a three-term
        decomposition would report a total of zero for a contrast that is not
        zero.
        """
        rows = [
            ("p1", "g0", "A", 0.25, 2.0), ("p1", "g0", "B", 0.75, 2.0),
            ("p2", "g0", "A", 0.75, 2.0), ("p2", "g0", "B", 0.25, 2.0),
            ("q1", "g1", "A", 0.25, 1.0), ("q1", "g1", "B", 0.75, 2.0),
            ("q2", "g1", "A", 0.75, 3.0), ("q2", "g1", "B", 0.25, 2.0),
        ]
        r = decompose(frame(rows), "g0", "g1")
        self.assertAlmostEqual(r["within_state"], 0.0, places=15)
        self.assertAlmostEqual(r["composition"], 0.0, places=15)
        self.assertAlmostEqual(r["interaction"], 0.0, places=15)
        self.assertNotAlmostEqual(r["covariance"], 0.0, places=6)
        # covariance in A only. cov_1(A) = mean[(0.25-0.5)(1-2), (0.75-0.5)(3-2)]
        #                                = mean[0.25, 0.25] = 0.25; cov_0(A) = 0
        # because mu is constant in g0. State B has constant mu in both groups.
        self.assertAlmostEqual(r["covariance_by_state"]["A"], 0.25, places=15)
        self.assertAlmostEqual(r["covariance_by_state"]["B"], 0.0, places=15)
        self.assertAlmostEqual(r["observed_total"], 0.25, places=15)
        self.assertLess(abs(r["closure_error"]), 1e-15)

    def test_population_convention_is_required(self):
        """The sample convention (divisor n-1) would not close.

        Same data as the covariance test. Recomputing the covariance term with
        divisor n-1 inflates it by n/(n-1) = 2 for n = 2, so the reconstructed
        total would overshoot the observed contrast by the covariance term itself.
        """
        rows = [
            ("p1", "g0", "A", 0.25, 2.0), ("p1", "g0", "B", 0.75, 2.0),
            ("p2", "g0", "A", 0.75, 2.0), ("p2", "g0", "B", 0.25, 2.0),
            ("q1", "g1", "A", 0.25, 1.0), ("q1", "g1", "B", 0.75, 2.0),
            ("q2", "g1", "A", 0.75, 3.0), ("q2", "g1", "B", 0.25, 2.0),
        ]
        r = decompose(frame(rows), "g0", "g1")
        pop = r["covariance"]
        sample = pop * 2.0 / (2.0 - 1.0)          # n = 2 per group
        wrong_total = (r["within_state"] + r["composition"] + r["interaction"]
                       + sample)
        self.assertGreater(abs(wrong_total - r["observed_total"]), 1e-6)

    def test_population_cov_helper(self):
        a = np.array([1.0, 2.0, 3.0, 4.0])
        b = np.array([2.0, 4.0, 6.0, 8.0])
        self.assertAlmostEqual(_population_cov(a, b), 2.5, places=15)
        self.assertAlmostEqual(_population_cov(a, a), np.var(a), places=15)


class TestGuards(unittest.TestCase):

    def test_missing_state_refuses(self):
        rows = [
            ("p1", "g0", "A", 0.5, 1.0), ("p1", "g0", "B", 0.5, 1.0),
            ("p2", "g0", "A", 1.0, 1.0),                      # no state B
            ("q1", "g1", "A", 0.5, 1.0), ("q1", "g1", "B", 0.5, 1.0),
        ]
        with self.assertRaises(ValueError) as cm:
            decompose(frame(rows), "g0", "g1")
        self.assertIn("missing at least one state", str(cm.exception))

    def test_fractions_must_sum_to_one(self):
        rows = [
            ("p1", "g0", "A", 0.3, 1.0), ("p1", "g0", "B", 0.3, 1.0),
            ("q1", "g1", "A", 0.5, 1.0), ("q1", "g1", "B", 0.5, 1.0),
        ]
        with self.assertRaises(ValueError) as cm:
            decompose(frame(rows), "g0", "g1")
        self.assertIn("do not sum to 1", str(cm.exception))

    def test_subset_of_states_allowed_when_declared(self):
        rows = [
            ("p1", "g0", "A", 0.3, 1.0), ("p1", "g0", "B", 0.3, 1.0),
            ("q1", "g1", "A", 0.5, 2.0), ("q1", "g1", "B", 0.5, 2.0),
        ]
        r = decompose(frame(rows), "g0", "g1", check_fractions=False)
        self.assertLess(abs(r["closure_error"]), 1e-15)

    def test_patient_in_two_groups_refuses(self):
        rows = [
            ("p1", "g0", "A", 0.5, 1.0), ("p1", "g0", "B", 0.5, 1.0),
            ("p1", "g1", "A", 0.5, 1.0), ("p1", "g1", "B", 0.5, 1.0),
            ("q1", "g1", "A", 0.5, 1.0), ("q1", "g1", "B", 0.5, 1.0),
        ]
        with self.assertRaises(ValueError) as cm:
            decompose(frame(rows), "g0", "g1")
        self.assertIn("more than one group", str(cm.exception))

    def test_missing_column_refuses(self):
        df = frame([("p1", "g0", "A", 0.5, 1.0)]).drop(columns=["state_mean"])
        with self.assertRaises(ValueError):
            decompose(df, "g0", "g1")


@unittest.skipUnless(os.path.exists(ARCHIVE),
                     f"archived patient-level table not found at {ARCHIVE}")
class TestArchivedRegression(unittest.TestCase):
    """Reproduce the values reported in the case study from the archived table."""

    @classmethod
    def setUpClass(cls):
        cls.df = pd.read_csv(ARCHIVE)
        cls.res = decompose(cls.df, "CI-NS", "CI-Sep")

    def test_analysis_set(self):
        self.assertEqual(self.res["n_patients"]["CI-NS"], 19)
        self.assertEqual(self.res["n_patients"]["CI-Sep"], 19)
        self.assertEqual(self.res["states"], ("CM", "Eff/EM"))

    def test_reported_terms(self):
        for key, want in REPORTED.items():
            with self.subTest(term=key):
                self.assertTrue(
                    math.isclose(self.res[key], want, rel_tol=0.0, abs_tol=1e-12),
                    f"{key}: got {self.res[key]!r}, reported {want!r}")

    def test_reported_covariance_by_state(self):
        for state, want in REPORTED_COV_BY_STATE.items():
            with self.subTest(state=state):
                self.assertTrue(
                    math.isclose(self.res["covariance_by_state"][state], want,
                                 rel_tol=0.0, abs_tol=1e-12))

    def test_closure_at_machine_precision(self):
        self.assertLess(abs(self.res["closure_error"]), 1e-15)

    def test_three_term_total_is_not_the_contrast(self):
        """The first three terms alone do not equal the observed contrast.

        This is the defect the fourth term exists to fix: their sum is
        0.0009025358558, more than twice the observed patient-weighted contrast
        of 0.0004228706566.
        """
        three = (self.res["within_state"] + self.res["composition"]
                 + self.res["interaction"])
        self.assertGreater(abs(three - self.res["observed_total"]), 1e-6)
        self.assertTrue(
            math.isclose(three, 0.0009025358558371542, rel_tol=0.0, abs_tol=1e-12),
            f"three-term sum is {three!r}")


if __name__ == "__main__":
    unittest.main(verbosity=2)
