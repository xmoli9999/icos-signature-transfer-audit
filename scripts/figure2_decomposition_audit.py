#!/usr/bin/env python3
"""Direct computation of the four-term patient-weighted decomposition (Figure 2D).

The fourth (between-patient covariance) term is computed from the patient-level
w_is and mu_is directly -- it is NOT obtained as observed_total minus the first
three terms. Population convention (divide by n) is used for all covariances,
which is what makes the identity close against the direct patient-weighted
contrast. Acceptance checks are asserted below.
"""
import gzip, os, numpy as np, pandas as pd

B = os.environ.get("ICOS_BASE", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
O = f"{B}/analysis_results"
STATES = ["CM", "Eff/EM"]
FROZEN3 = (0.003653230902695344, -0.0028159006193097243, 6.520557245153402e-05)

sc = pd.read_csv(gzip.open(f"{O}/gate1_cell_scores_calibrated.csv.gz", "rt"))
pat = pd.read_csv(f"{O}/gate1_patient_level_unblinded.csv")[["patient", "group"]]
sc = sc.merge(pat, on="patient", how="left")
sc["state"] = sc["celltype_fine"].map({"CD4 T CM": "CM", "CD4 T Eff/EM": "Eff/EM"})
sc = sc[sc.state.notna()]
sc["group"] = sc["group"].replace({"ICU non-sepsis": "CI-NS", "ICU sepsis": "CI-Sep"})
sc = sc[sc.group.isin(["CI-NS", "CI-Sep"])]

g = (sc.groupby(["patient", "group", "state"])["delta_ucell"]
       .agg(mu="mean", n="size").reset_index())
W = g.pivot_table(index=["patient", "group"], columns="state", values="n", fill_value=0)
M = g.pivot_table(index=["patient", "group"], columns="state", values="mu")
assert (W > 0).all().all() and M.notna().all().all(), "a patient lacks one state; evaluable set would change"
W = W.div(W.sum(axis=1), axis=0)

def split(grp):
    return W.xs(grp, level="group"), M.xs(grp, level="group")
W0, M0 = split("CI-NS"); W1, M1 = split("CI-Sep")
w0, m0, w1, m1 = W0.mean(), M0.mean(), W1.mean(), M1.mean()

within = sum(w0[s] * (m1[s] - m0[s]) for s in STATES)
comp   = sum(m0[s] * (w1[s] - w0[s]) for s in STATES)
inter  = sum((w1[s] - w0[s]) * (m1[s] - m0[s]) for s in STATES)

def pcov(Wg, Mg, s):                      # population covariance, divisor n
    a, b = Wg[s].values, Mg[s].values
    return float(np.mean((a - a.mean()) * (b - b.mean())))
cov_by_state = {s: pcov(W1, M1, s) - pcov(W0, M0, s) for s in STATES}
covariance = sum(cov_by_state.values())

recon = within + comp + inter + covariance
observed = (sum(float((W1[s] * M1[s]).mean()) for s in STATES)
            - sum(float((W0[s] * M0[s]).mean()) for s in STATES))
closure = recon - observed

for name, got, want in zip(("within", "composition", "interaction"), (within, comp, inter), FROZEN3):
    assert abs(got - want) < 1e-12, f"{name} does not reproduce the frozen value: {got} vs {want}"
assert abs(closure) < 1e-12, f"closure error too large: {closure}"

rows = [("within_state_at_NS_weights", within),
        ("composition_at_NS_levels", comp),
        ("interaction", inter),
        ("between_patient_covariance", covariance),
        ("reconstructed_total", recon),
        ("observed_patient_weighted_total", observed),
        ("closure_error", closure),
        ("covariance_CM", cov_by_state["CM"]),
        ("covariance_Eff/EM", cov_by_state["Eff/EM"]),
        ("n_patients_CI_NS", float(len(W0))),
        ("n_patients_CI_Sep", float(len(W1)))]
pd.DataFrame(rows, columns=["component", "estimate"]).to_csv(
    f"{O}/figure2_formal_decomposition.csv", index=False)
for k, v in rows: print(f"{k:34s} {v!r}")
