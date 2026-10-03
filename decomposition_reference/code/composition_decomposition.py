"""Patient-weighted four-term decomposition of a between-group contrast in a
state-mixed cell population.

This is a reusable reference implementation of the decomposition reported in the
accompanying case study. It is not a new estimator and makes no claim of
methodological priority; it exists so that the decomposition can be recomputed,
checked and reused without reconstructing it from the paper.

The identity
------------
For two groups 0 and 1, patients i, and cell states s, let w_is be the proportion
of state s among patient i's cells in the population of interest and mu_is be
that patient's mean score within state s. Write a bar for the group mean across
patients. The patient-weighted between-group contrast

    Delta = sum_s mean_1(w_s * mu_s) - sum_s mean_0(w_s * mu_s)

decomposes exactly as

    Delta =   sum_s wbar_0s (mubar_1s - mubar_0s)              [within_state]
            + sum_s mubar_0s (wbar_1s - wbar_0s)               [composition]
            + sum_s (wbar_1s - wbar_0s)(mubar_1s - mubar_0s)   [interaction]
            + sum_s {Cov_1(w_s, mu_s) - Cov_0(w_s, mu_s)}      [covariance]

Conventions that the identity depends on
----------------------------------------
1. Patient-weighted, not cell-weighted. Every patient contributes equally to each
   group mean, regardless of how many cells that patient contributed. A
   cell-weighted contrast is a different estimand and will not close against this
   identity.
2. Covariances use the POPULATION convention (divisor n, not n - 1). This is what
   makes the identity close exactly; the sample convention leaves a residual of
   order 1/(n-1).
3. `state_mean` is the arithmetic MEAN of the per-cell scores within that
   patient and state. A median, or any other non-linear summary, breaks the
   identity because the fourth term is derived from the algebra of means.
4. The covariance term is computed DIRECTLY from the patient-level w_is and
   mu_is. It is never obtained by subtracting the first three terms from the
   observed contrast; doing so would make `closure_error` identically zero and
   destroy the only check that the decomposition is correct.
5. Every patient must have a value for every state. If a patient is missing a
   state, the evaluable set is not the same as the set the contrast was computed
   on, and this implementation refuses to proceed rather than silently dropping
   or imputing.

`closure_error` is the diagnostic, not a fitted residual: it is
`reconstructed_total - observed_total` and should be at machine precision.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

__all__ = ["decompose", "DEFAULT_CLOSURE_TOL", "REQUIRED_COLUMNS"]

DEFAULT_CLOSURE_TOL = 1e-12

REQUIRED_COLUMNS = ("patient_id", "group", "state", "state_fraction", "state_mean")


def _population_cov(a: np.ndarray, b: np.ndarray) -> float:
    """Covariance with divisor n (population convention)."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    return float(np.mean((a - a.mean()) * (b - b.mean())))


def decompose(
    data: pd.DataFrame,
    group_0,
    group_1,
    *,
    closure_tol: float = DEFAULT_CLOSURE_TOL,
    check_fractions: bool = True,
) -> dict:
    """Decompose the patient-weighted contrast between two groups.

    Parameters
    ----------
    data
        Long-form patient-level table with one row per (patient, state) and the
        columns in ``REQUIRED_COLUMNS``:

        ``patient_id``      patient or sample identifier; the inferential unit
        ``group``           group label; exactly two are used per call
        ``state``           cell-state label within the population of interest
        ``state_fraction``  w_is, that patient's proportion of cells in that state
        ``state_mean``      mu_is, that patient's MEAN score within that state

        Extra columns are ignored. A patient must appear in exactly one group.
    group_0, group_1
        The reference group and the comparison group. The contrast returned is
        group_1 minus group_0.
    closure_tol
        Maximum absolute closure error tolerated before ``AssertionError``.
    check_fractions
        If True, require each patient's state fractions to sum to 1. Set False
        when the states given are a subset of the population on purpose; the
        decomposition is still exact for the contrast restricted to those states.

    Returns
    -------
    dict with keys
        ``within_state``, ``composition``, ``interaction``, ``covariance``,
        ``reconstructed_total``, ``observed_total``, ``closure_error``,
        ``covariance_by_state`` (dict), ``states`` (tuple),
        ``n_patients`` (dict keyed by group label).

    Raises
    ------
    ValueError
        On missing columns, an unusable group selection, a patient in more than
        one group, a patient missing a state, or duplicated (patient, state) rows.
    AssertionError
        If the identity fails to close within ``closure_tol``. This signals a
        convention violation or a data problem, not a modelling choice.
    """
    missing = [c for c in REQUIRED_COLUMNS if c not in data.columns]
    if missing:
        raise ValueError(f"missing required column(s): {missing}")
    if group_0 == group_1:
        raise ValueError("group_0 and group_1 must differ")

    df = data.loc[data["group"].isin([group_0, group_1]),
                  list(REQUIRED_COLUMNS)].copy()
    for g in (group_0, group_1):
        if not (df["group"] == g).any():
            raise ValueError(f"group {g!r} is not present in the data")

    per_patient_groups = df.groupby("patient_id")["group"].nunique()
    if (per_patient_groups > 1).any():
        bad = per_patient_groups[per_patient_groups > 1].index.tolist()
        raise ValueError(f"patient(s) appear in more than one group: {bad}")

    if df.duplicated(["patient_id", "state"]).any():
        dup = df.loc[df.duplicated(["patient_id", "state"], keep=False), "patient_id"].unique()
        raise ValueError(f"duplicated (patient_id, state) rows for: {list(dup)}")

    states = tuple(sorted(df["state"].unique()))
    if len(states) < 2:
        raise ValueError(f"at least two states are required; found {states}")

    W = df.pivot_table(index=["patient_id", "group"], columns="state",
                       values="state_fraction")
    M = df.pivot_table(index=["patient_id", "group"], columns="state",
                       values="state_mean")
    if W.isna().any().any() or M.isna().any().any():
        incomplete = sorted(set(W[W.isna().any(axis=1)].index.get_level_values(0))
                            | set(M[M.isna().any(axis=1)].index.get_level_values(0)))
        raise ValueError(
            "patient(s) missing at least one state, so the evaluable set is not "
            f"the set the contrast is computed on: {incomplete}")
    if check_fractions:
        tot = W[list(states)].sum(axis=1)
        if not np.allclose(tot.values, 1.0, atol=1e-9):
            off = tot[(tot - 1.0).abs() > 1e-9]
            raise ValueError(
                "state fractions do not sum to 1 for patient(s): "
                f"{off.index.get_level_values(0).tolist()}; pass "
                "check_fractions=False if the states are deliberately a subset")

    W0 = W.xs(group_0, level="group"); M0 = M.xs(group_0, level="group")
    W1 = W.xs(group_1, level="group"); M1 = M.xs(group_1, level="group")
    w0, m0 = W0.mean(), M0.mean()
    w1, m1 = W1.mean(), M1.mean()

    within = float(sum(w0[s] * (m1[s] - m0[s]) for s in states))
    comp = float(sum(m0[s] * (w1[s] - w0[s]) for s in states))
    inter = float(sum((w1[s] - w0[s]) * (m1[s] - m0[s]) for s in states))

    cov_by_state = {
        s: _population_cov(W1[s].values, M1[s].values)
           - _population_cov(W0[s].values, M0[s].values)
        for s in states
    }
    covariance = float(sum(cov_by_state.values()))

    reconstructed = within + comp + inter + covariance
    observed = float(sum(float((W1[s] * M1[s]).mean()) for s in states)
                     - sum(float((W0[s] * M0[s]).mean()) for s in states))
    closure = reconstructed - observed

    assert abs(closure) < closure_tol, (
        f"the identity did not close: closure_error={closure!r} exceeds "
        f"{closure_tol!r}. Check that state_mean is an arithmetic mean, that the "
        "contrast is patient-weighted, and that no patient is missing a state.")

    return {
        "within_state": within,
        "composition": comp,
        "interaction": inter,
        "covariance": covariance,
        "reconstructed_total": reconstructed,
        "observed_total": observed,
        "closure_error": closure,
        "covariance_by_state": cov_by_state,
        "states": states,
        "n_patients": {group_0: int(len(W0)), group_1: int(len(W1))},
    }
