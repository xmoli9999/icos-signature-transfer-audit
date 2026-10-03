# Patient-weighted four-term decomposition — reference implementation

**What this is.** A reusable reference implementation of the decomposition
reported in the accompanying case study, so that the decomposition can be
recomputed, checked and reused without reconstructing it from the paper.

**What this is not.** It is not a new estimator, a new decomposition method, or a
software package with a claim of methodological priority. Variance and mean
decompositions of this form are standard algebra; what is provided here is a
tested implementation with the conventions the identity depends on written down
and enforced, together with a regression test against the archived values the
paper reports.

## The identity

For two groups 0 and 1, patients *i*, and cell states *s*, let *w<sub>is</sub>* be
the proportion of state *s* among patient *i*'s cells in the population of
interest, and *μ<sub>is</sub>* that patient's mean score within state *s*. A bar
denotes the group mean across patients. The patient-weighted between-group
contrast

    Δ = Σ_s mean_1(w_s·μ_s) − Σ_s mean_0(w_s·μ_s)

decomposes exactly as

    Δ =   Σ_s w̄_0s (μ̄_1s − μ̄_0s)                    within_state
        + Σ_s μ̄_0s (w̄_1s − w̄_0s)                    composition
        + Σ_s (w̄_1s − w̄_0s)(μ̄_1s − μ̄_0s)            interaction
        + Σ_s { Cov_1(w_s, μ_s) − Cov_0(w_s, μ_s) }   covariance

The fourth term is required because the group mean of a per-patient product is
not the product of the group means. Dropping it does not give an approximation
with a small error: in the case study the first three terms sum to 0.0009025,
more than twice the observed patient-weighted contrast of 0.0004229.

## Conventions the identity depends on

These are enforced in code and covered by tests. Violating any of them breaks
closure.

1. **Patient-weighted, not cell-weighted.** Each patient contributes equally to
   each group mean regardless of cell count. A cell-weighted contrast is a
   different estimand and will not close against this identity.
2. **Population covariance convention (divisor *n*, not *n* − 1).** The sample
   convention leaves a residual of order 1/(*n* − 1). A test demonstrates the
   failure explicitly.
3. **`state_mean` is the arithmetic mean** of the per-cell scores within that
   patient and state. A median or any other non-linear summary breaks the
   identity, because the fourth term comes from the algebra of means.
4. **The covariance term is computed directly** from the patient-level
   *w<sub>is</sub>* and *μ<sub>is</sub>*. It is never obtained by subtracting the
   first three terms from the observed contrast. Obtaining it by subtraction
   would force `closure_error` to zero and destroy the only check that the
   decomposition is right.
5. **Every patient must have every state.** If a patient is missing a state, the
   evaluable set is not the set the contrast was computed on. The implementation
   raises rather than dropping or imputing.

`closure_error` is `reconstructed_total − observed_total`. It is a diagnostic,
not a fitted residual, and should sit at machine precision. The default
tolerance is 1e-12.

## Layout

    code/composition_decomposition.py            the implementation
    tests/test_composition_decomposition.py      analytic + archived regression tests
    examples/reproduce_figure2_decomposition.py  reproduces the reported values
    data/figure2_patient_state_table.csv         archived patient-level input
    data/figure2_formal_decomposition.csv        archived reported values

## Input format

Long form, one row per (patient, state):

| column | meaning |
|---|---|
| `patient_id` | patient or sample identifier; the inferential unit |
| `group` | group label; exactly two are used per call |
| `state` | cell-state label within the population of interest |
| `state_fraction` | *w<sub>is</sub>*, that patient's proportion of cells in that state |
| `state_mean` | *μ<sub>is</sub>*, that patient's **mean** score within that state |

Extra columns are ignored. A patient must appear in exactly one group.

## Use

```python
import pandas as pd
from composition_decomposition import decompose

df = pd.read_csv("data/figure2_patient_state_table.csv")
res = decompose(df, "CI-NS", "CI-Sep")
print(res["composition"], res["covariance"], res["closure_error"])
```

Returned keys: `within_state`, `composition`, `interaction`, `covariance`,
`reconstructed_total`, `observed_total`, `closure_error`,
`covariance_by_state`, `states`, `n_patients`.

## Running the tests

No test framework is required beyond the standard library:

    python3 -m unittest discover -s tests -t .

`pytest tests` also works if pytest is installed. Fifteen tests run in two
layers:

- **Analytic toy tests** on tiny hand-checkable datasets: a pure composition
  shift, a pure within-state shift, a case where the covariance term carries the
  entire contrast, an explicit demonstration that the sample covariance
  convention fails to close, and the covariance helper itself. Guard tests cover
  a patient missing a state, fractions that do not sum to one, a deliberate
  state subset, a patient appearing in two groups, and a missing column.
- **Archived regression test** against the values the case study reports,
  recomputed from `data/figure2_patient_state_table.csv`: all four terms, both
  per-state covariances, the analysis set (19 and 19 patients, states CM and
  Eff/EM), closure below 1e-15, and an assertion that the first three terms alone
  do not equal the observed contrast.

The regression test uses an absolute tolerance of 1e-12. The recomputed values
agree with the archived ones to about 7e-17, not bit-for-bit: the reference
implementation sums from the long-form table while the original audit script
pivots on cell counts, so floating-point summation order differs. That
difference is fourteen orders of magnitude below the smallest reported term.

## Provenance of the bundled data

`data/figure2_patient_state_table.csv` is produced by
`../scripts/export_patient_state_table.py` from the archived calibrated
cell-level scores and the archived patient-level group assignments. That script
performs no new analysis: its grouping, state mapping and group relabelling are
the same operations as in `../scripts/figure2_decomposition_audit.py`. The export
exists so that neither the tests nor the example needs the cell-level archive.

`data/figure2_formal_decomposition.csv` is a copy of the archived result file the
manuscript cites.

## Licence

MIT, as for the rest of the project code.
