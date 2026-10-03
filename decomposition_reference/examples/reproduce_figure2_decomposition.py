#!/usr/bin/env python3
"""Reproduce the four-term decomposition reported in the case study.

Reads the archived patient-level state table and prints every component, then
compares each against the archived result file the manuscript cites.

    python3 examples/reproduce_figure2_decomposition.py

Optional environment variables:
    ICOS_PATIENT_STATE_TABLE  path to figure2_patient_state_table.csv
    ICOS_REPORTED_CSV         path to figure2_formal_decomposition.csv
"""
import os
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "code"))
from composition_decomposition import decompose  # noqa: E402

TABLE = os.environ.get("ICOS_PATIENT_STATE_TABLE",
                       os.path.join(ROOT, "data", "figure2_patient_state_table.csv"))
REPORTED = os.environ.get("ICOS_REPORTED_CSV",
                          os.path.join(ROOT, "data", "figure2_formal_decomposition.csv"))

# component name in the archived CSV -> key returned by decompose()
MAP = {
    "within_state_at_NS_weights": "within_state",
    "composition_at_NS_levels": "composition",
    "interaction": "interaction",
    "between_patient_covariance": "covariance",
    "reconstructed_total": "reconstructed_total",
    "observed_patient_weighted_total": "observed_total",
    "closure_error": "closure_error",
}

df = pd.read_csv(TABLE)
res = decompose(df, "CI-NS", "CI-Sep")

print(f"input            : {TABLE}")
print(f"states           : {', '.join(res['states'])}")
print(f"patients         : " + ", ".join(f"{g}={n}" for g, n in res["n_patients"].items()))
print()
for key in ("within_state", "composition", "interaction", "covariance",
            "reconstructed_total", "observed_total", "closure_error"):
    print(f"{key:24s} {res[key]!r}")
for state, v in res["covariance_by_state"].items():
    print(f"{'covariance_' + state:24s} {v!r}")

if not os.path.exists(REPORTED):
    print(f"\nreported file not found at {REPORTED}; skipped comparison")
    raise SystemExit(0)

rep = pd.read_csv(REPORTED).set_index("component")["estimate"].to_dict()
print("\ncomparison against the archived reported values")
print(f"{'component':32s} {'recomputed':>24s} {'reported':>24s} {'abs diff':>12s}  ok")
worst = 0.0
for comp, key in MAP.items():
    got, want = res[key], float(rep[comp])
    d = abs(got - want)
    worst = max(worst, d)
    print(f"{comp:32s} {got:24.17g} {want:24.17g} {d:12.3g}  {'yes' if d < 1e-12 else 'NO'}")
for state, key in (("CM", "covariance_CM"), ("Eff/EM", "covariance_Eff/EM")):
    got, want = res["covariance_by_state"][state], float(rep[key])
    d = abs(got - want)
    worst = max(worst, d)
    print(f"{key:32s} {got:24.17g} {want:24.17g} {d:12.3g}  {'yes' if d < 1e-12 else 'NO'}")

print(f"\nlargest absolute difference: {worst:.3g}")
print(f"closure error              : {res['closure_error']:.3g}")
three = res["within_state"] + res["composition"] + res["interaction"]
print(f"first three terms only     : {three!r}")
print(f"observed contrast          : {res['observed_total']!r}")
print("the first three terms alone overstate the contrast by "
      f"{abs(three - res['observed_total']):.3g}, which is what the covariance "
      "term accounts for")
raise SystemExit(0 if worst < 1e-12 else 1)
