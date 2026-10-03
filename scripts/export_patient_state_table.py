#!/usr/bin/env python3
"""Export the patient-level state table that the four-term decomposition consumes.

This is an export of an intermediate that figure2_decomposition_audit.py already
computes from archived inputs. It performs no new analysis: the grouping, the
state mapping and the group relabelling are byte-for-byte the same operations.
Output is the long-form table used by the reference implementation and its
regression test, so that neither needs the cell-level score archive.
"""
import gzip, os, pandas as pd

B = os.environ.get("ICOS_BASE", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
O = f"{B}/analysis_results"

sc = pd.read_csv(gzip.open(f"{O}/gate1_cell_scores_calibrated.csv.gz", "rt"))
pat = pd.read_csv(f"{O}/gate1_patient_level_unblinded.csv")[["patient", "group"]]
sc = sc.merge(pat, on="patient", how="left")
sc["state"] = sc["celltype_fine"].map({"CD4 T CM": "CM", "CD4 T Eff/EM": "Eff/EM"})
sc = sc[sc.state.notna()]
sc["group"] = sc["group"].replace({"ICU non-sepsis": "CI-NS", "ICU sepsis": "CI-Sep"})
sc = sc[sc.group.isin(["CI-NS", "CI-Sep"])]

g = (sc.groupby(["patient", "group", "state"])["delta_ucell"]
       .agg(state_mean="mean", n_cells="size").reset_index())
g["state_fraction"] = g["n_cells"] / g.groupby("patient")["n_cells"].transform("sum")
g = g.rename(columns={"patient": "patient_id"})
g = g[["patient_id", "group", "state", "state_fraction", "state_mean", "n_cells"]]
g = g.sort_values(["group", "patient_id", "state"]).reset_index(drop=True)

# integrity: every patient has both states, and fractions sum to 1 per patient
piv = g.pivot_table(index="patient_id", columns="state", values="state_fraction")
assert piv.notna().all().all(), "a patient lacks one state"
assert (piv.sum(axis=1).sub(1.0).abs() < 1e-12).all(), "state fractions do not sum to 1"

out = f"{O}/figure2_patient_state_table.csv"
g.to_csv(out, index=False)
print(f"wrote {out}: {len(g)} rows, {g.patient_id.nunique()} patients")
print(g.groupby("group").patient_id.nunique().to_string())
