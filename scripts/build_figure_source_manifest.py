#!/usr/bin/env python3
"""Build FIGURE_SOURCE_MANIFEST.tsv: every figure panel -> its machine source.

provenance_class:
  machine_read      the generating script reads this file at run time; the plotted
                    or annotated value is never typed into the script
  literal_traceable the value is hard-coded in the script, but the same value is
                    present verbatim in the named archived file
  literal_derived   the value is hard-coded in the script and is obtained from the
                    named archived file by the arithmetic stated in `note`
  text_only         the panel carries no data value
Checksums are computed here, never typed.
"""
import hashlib, os, sys

B = os.environ.get("ICOS_BASE", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

ROWS = [
 # figure, panel, source_file, generating_script, key_variables, provenance_class, note
 ("Figure 1","A","metadata/ortholog_attrition.csv","manuscript/make_fig1.py",
  "stage; n","literal_derived",
  "155 and 107 are rows of this file; the displayed -48 is 155 minus 107"),
 ("Figure 1","B","gate1_archive_2026-09-15/provenance/analysis_plan_frozen.md","manuscript/make_fig1.py",
  "ambient retention; sign concordance; Sham1-excluded overlap; leave-one-rat-out AUROC","literal_derived",
  "retention 0.974 and sign 0.983 appear verbatim (ambient exclusion line); overlap 0.884 appears as 88.4%; "
  "the plotted 0.896 is the mean of the three held-out AUROCs 0.928/0.872/0.889 recorded for G0B-2, "
  "not a single AUROC -- the panel label does not currently say this"),
 ("Figure 1","B","metadata/gate1_SAP.md","manuscript/make_fig1.py",
  "prespecified thresholds 0.60 and 0.70","literal_traceable",
  "the two orange threshold ticks; both floors are stated in the frozen plan and in Methods"),

 ("Figure 2","A","analysis_results/gate1_gate_test_results.json","manuscript/make_figures.py",
  "primary.beta/ci_lo/ci_hi; decomp_CD4 T CM; decomp_CD4 T Eff/EM","machine_read",
  "read via _fig_common.py as res[...]"),
 ("Figure 2","B","analysis_results/gate1_cell_scores_calibrated.csv.gz","manuscript/make_figures.py",
  "patient; celltype_fine; delta_ucell","machine_read",
  "per-patient state means mCM and mEM computed in _fig_common.py"),
 ("Figure 2","B","analysis_results/gate1_cell_scores.csv.gz","manuscript/make_figures.py",
  "cell_index; group; pop_primary_memory","machine_read",
  "supplies the group label and the primary-memory mask"),
 ("Figure 2","C","analysis_results/gate1_cell_scores_calibrated.csv.gz","manuscript/make_figures.py",
  "per-patient Eff/EM fraction","machine_read",
  "the plotted points are computed, not typed"),
 ("Figure 2","C","analysis_results/figure2_state_proportions.csv","manuscript/make_figures.py",
  "Eff/EM_prop","literal_traceable",
  "the annotation '0.607 -> 0.420' is hard-coded; the group means of Eff/EM_prop in this file are "
  "0.606651 (CI-NS) and 0.420310 (CI-Sep)"),
 ("Figure 2","D","analysis_results/figure2_formal_decomposition.csv","manuscript/make_figures.py",
  "within_state_at_NS_weights; composition_at_NS_levels; interaction; between_patient_covariance; "
  "observed_patient_weighted_total","machine_read",
  "all five bars and the footnote total are read from this file"),
 ("Figure 2","D","scripts/figure2_decomposition_audit.py","manuscript/make_figures.py",
  "-- generator of the source table --","machine_read",
  "computes the four terms directly from patient-level w_is and mu_is; asserts closure below 1e-12"),

 ("Figure 3","A","analysis_results/gate1_cell_scores_calibrated.csv.gz","manuscript/make_fig34.py",
  "nFeature_RNA; ucell_sig","machine_read","cell-level scatter, all donors"),
 ("Figure 3","B","analysis_results/gate1_cell_scores_calibrated.csv.gz","manuscript/make_fig34.py",
  "nFeature_RNA; delta_ucell","machine_read","cell-level scatter, all donors"),
 ("Figure 3","A,B,C","analysis_results/figure3_depth_correlations.csv","manuscript/make_fig34.py",
  "cell_raw; cell_calibrated; donor_raw; donor_calibrated; within_group_raw; within_group_calibrated",
  "machine_read","every rho printed on the panels is read from this file"),
 ("Figure 3","C","analysis_results/gate1_patient_level_unblinded.csv","manuscript/make_fig34.py",
  "med_nfeat; mean_raw; mean_delta; group","machine_read","donor-level scatter, 47 donors"),
 ("Figure 3","D","analysis_results/gate1_calibration_diagnostics.json","manuscript/make_fig34.py",
  "signature and matched-control detection frequency and mean expression","machine_read",
  "control-matching bars"),
 ("Figure 3","A,B,C","scripts/figure3_depth_correlation_audit.py","manuscript/make_fig34.py",
  "-- generator of the source table --","machine_read",
  "recomputes all six rho values and asserts each against the previously hard-coded value"),

 ("Figure 4","A","analysis_results/gate1_patient_level_unblinded.csv","manuscript/make_fig34.py",
  "group; mean_delta","machine_read","per-patient points and group means"),
 ("Figure 4","A","analysis_results/gate1_patient_level_unblinded.csv","manuscript/make_fig34.py",
  "group counts 9 / 19 / 19","literal_derived",
  "the three x-axis tick labels are hard-coded; they are the group value counts of this file"),
 ("Figure 4","B","analysis_results/gate1_gate_test_results.json","manuscript/make_fig34.py",
  "primary; sens1_unadjusted_delta; sens2_raw_adjusted; sens3_aucell_adjusted","machine_read",
  "all four estimates and intervals"),
 ("Figure 4","C","metadata/gate2_preregistration_DRAFT.md","manuscript/make_fig34.py",
  "repositories searched; closest supportive cohorts; cohorts meeting all five criteria","literal_traceable",
  "the bar values 8 / 2 / 0 are hard-coded; the frozen plan enumerates exactly eight repositories "
  "(GEO, Broad Single Cell Portal, cellxgene Discover, dbGaP, EGA, GSA-Human/HRA, "
  "ArrayExpress/BioStudies, Zenodo) and states that none satisfies all five criteria"),
 ("Figure 4","D","metadata/gate2b_preregistration_FROZEN_v1.0.md","manuscript/make_fig34.py",
  "Kwok and Reyes qualifying patients per arm","literal_traceable",
  "the fractions 14/26, 3/8, 4/7, 7/7 are hard-coded; the frozen Gate 2b plan and Supplementary "
  "Table S3 carry the same counts"),

 ("Figure S1","--","(none)","manuscript/make_figS1.py",
  "-- no data value --","text_only",
  "records the frozen model specification and its not-estimable status; contains no estimate"),
]

def sha(path):
    p = os.path.join(B, path)
    if not os.path.exists(p):
        return "FILE_NOT_FOUND"
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()

out = [("figure", "panel", "source_file", "generating_script", "key_variables",
        "provenance_class", "source_sha256", "generating_script_sha256", "note")]
missing = []
for fig, panel, src, script, keys, cls, note in ROWS:
    s1 = "n/a" if src == "(none)" else sha(src)
    s2 = sha(script)
    if "FILE_NOT_FOUND" in (s1, s2):
        missing.append((fig, panel, src, script))
    out.append((fig, panel, src, script, keys, cls, s1, s2, note))

dest = os.path.join(B, "FIGURE_SOURCE_MANIFEST.tsv")
with open(dest, "w", encoding="utf-8") as f:
    for r in out:
        f.write("\t".join(r) + "\n")

print(f"wrote {dest}: {len(out) - 1} rows")
from collections import Counter
print(Counter(r[5] for r in out[1:]))
if missing:
    print("MISSING FILES:", missing, file=sys.stderr)
    sys.exit(1)
