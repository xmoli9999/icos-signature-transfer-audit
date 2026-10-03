# Analysis archive

Code, frozen analysis plans, archived results and figure provenance for the
manuscript *Cell-state mixture and detection depth can change what a transferred
single-cell signature reports*.

Licence: MIT (see `LICENSE`). Written for deposition on Zenodo and mirroring on
GitHub, as required by the journal's open-source and data-availability policy.

Repository: https://github.com/xmoli9999/icos-signature-transfer-audit  
Archived version of record: https://doi.org/10.5281/zenodo.23118493

## What this is, and what it is not

This study is a **secondary analysis of previously published, publicly available
data**. No primary data are redistributed here. Every source dataset is obtained
from its own repository:

| Dataset | Role | Where to get it |
|---|---|---|
| GSE285325 | rat lymph discovery; source of the frozen 107-gene signature | NCBI GEO |
| GSE290679 | human Gate 1 cohort | NCBI GEO |
| Kwok et al. whole blood | supportive Gate 2b cohort | Zenodo record 7924238 (v1.1) |
| Reyes et al. | supportive Gate 2b cohort | Single Cell Portal, SCP548 |
| Hao et al. PBMC multimodal reference | reference mapping | CELLxGENE, nygc multimodal PBMC collection |

What is here is everything needed to follow and re-run the analysis once those
inputs are in place, plus the frozen plans that fix what was decided before the
data were seen.

## Layout

The tree deliberately mirrors the working project directory, so the archived
scripts run unmodified — they resolve paths from `ICOS_BASE` and expect these
folder names.

```
scripts/                     analysis and audit scripts
manuscript/                  figure-generation scripts
decomposition_reference/     reference implementation of the four-term
                             decomposition, with its own README and tests
metadata/                    frozen analysis plans, the 107-gene signature,
                             ortholog attrition, label mappings
analysis_results/            archived intermediate and final result objects
figures/                     Figures 1-4 and S1 as PDF, SVG and 600 dpi PNG
provenance/                  figure-source manifest, the Kwok deposition check,
                             and the Gate 1 frozen archive
CHECKSUMS.sha256             sha256 for every file in this archive
```

## Requirements

No installation step is needed: the archive is run in place, from its own root.

Python 3.10 or later. The audit and figure scripts import `numpy`, `pandas`,
`scipy`, `statsmodels` and `matplotlib`; the versions these results were produced
and re-verified under are Python 3.10.12, numpy 2.2.6, pandas 2.3.3 and scipy
1.15.3. The reference implementation in `decomposition_reference/` needs only
`numpy` and `pandas`, and its test suite needs nothing beyond the standard
library.

```
pip install numpy pandas scipy statsmodels matplotlib
```

Every command below is run from the root of this archive.

## Reproducing the reported numbers

The four-term decomposition and the six detection-depth correlations are the two
places where a reader is most likely to want to check the arithmetic. Both have
standalone audit scripts that recompute from archived inputs and assert against the
reported values:

```
export ICOS_BASE=/path/to/this/archive
python3 scripts/figure2_decomposition_audit.py
python3 scripts/figure3_depth_correlation_audit.py
```

Both read `ICOS_BASE` from the environment; set it to the root of this archive. Both were run from inside this archive, with no other inputs, before it was deposited: they recompute the values and their assertions pass.

The decomposition also has a dependency-free test suite:

```
cd decomposition_reference && python3 -m unittest discover -s tests -t .
```

Fifteen tests, no third-party test framework required; all pass from inside this archive. They cover analytic toy
cases, the input guards, and a regression test against the archived values.

## Figure provenance

`provenance/FIGURE_SOURCE_MANIFEST.tsv` maps every figure panel to the file it
reads and the script that draws it, with checksums, and classifies each row:

- `machine_read` — the script reads the value at run time (14 rows)
- `literal_traceable` — hard-coded in the script, same value present in the named
  archived file (4 rows)
- `literal_derived` — hard-coded, obtained from the named file by the arithmetic
  stated in the row (3 rows)
- `text_only` — the panel carries no data value (1 row)

Not every figure value is machine-read, and the manifest exists to say which are
not rather than to imply they all are.

Compare PNG checksums across runs. PDF and SVG exports embed a creation timestamp
and change checksum even when the figure itself has not.

## Fixed stochastic parameters

| Parameter | Value |
|---|---|
| Control-set construction seed | 20260915 |
| Matched control sets | 100 |
| Detection bins / expression bins | 20 / 20 |
| Control gene universe | 28,167 genes |
| UCell MAXRANK | 1500 |
| Figure 3A,B plotting subsample | 25,000 cells, `random_state = 3` |
| Closure tolerance, four-term identity | 1e-12 (observed closure error -2.11e-18) |

## Deliberately not included

- `metadata/scp_meta_updated.txt` — third-party cohort metadata from the Single
  Cell Portal. Obtain it from SCP548 rather than from here.
- Files marked `.OLD_DO_NOT_USE`, `_SUPERSEDED` or `.bak_*` in the working project
  directory. They are superseded working copies; including them would invite someone
  to run the wrong one. The frozen plans that do govern the analysis are in
  `metadata/` and `provenance/gate1_archive_2026-09-15/`.

## A note on the frozen plans

The plans and amendments in `metadata/` were each fixed in writing before the
analysis that used them, and the blinded amendments were finalized before any group
comparison. The archive itself, however, was created after Gate 1 was unblinded. Its
deposition timestamp therefore locks the completed Gate 1 analysis and everything
after it; it is not evidence of prospective registration of Gate 1. This is stated
the same way in the manuscript.
