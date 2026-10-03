# Kwok deposition provenance check

Run 2026-09-30. Closes the loop from the files actually on disk to Zenodo record
7924238, and records the relationship between records 7723202 and 7924238.

## Result

**PASS, 4/4 files.** Byte sizes and md5 checksums of the four local files agree
exactly with (a) the md5 table recorded in Appendix F of the frozen Gate 2b plan
and (b) the checksums published on Zenodo record 7924238. Per-file detail is in
`KWOK_ZENODO_PROVENANCE_CHECK.tsv`.

Local files, at `Desktop/Kwok2023_whole_blood/source/`:

| File | Bytes | md5 |
|---|---:|---|
| `..._BD-Rhapsody_whole-blood_RNA-counts.mtx` | 3,494,797,339 | `719ef6a602474e78c2331e39be9bee02` |
| `..._BD-Rhapsody_whole-blood_RNA-counts_barcodes.tsv` | 9,270,174 | `9812e0d17d1a25c2b65006edfc068e86` |
| `..._BD-Rhapsody_whole-blood_RNA-counts_features.tsv` | 177,721 | `f927042500a6b4e9f3e4344b87815626` |
| `..._BD-Rhapsody_whole-blood_cell-metadata.tsv` | 151,531,112 | `3c71e9e7947bd8a5b962b4d4a7ba4664` |

md5 values were computed from the files on disk in this check; they were not
copied from the frozen plan and then compared against themselves.

## Record relationship, as observed on Zenodo

- Record **7723202** (DOI 10.5281/zenodo.7723202) carries a "newer version
  available" notice that resolves to record **7924238**
  (DOI 10.5281/zenodo.7924238). Both belong to the deposition for Kwok et al.
  2023, *Neutrophils and emergency granulopoiesis drive immune suppression and an
  extreme response endotype during sepsis*. They are the same dataset's version
  chain, not two cohorts.
- Record 7723202 lists **7 files** and does **not** include the whole-blood
  barcodes or features lists. Of the four files this study uses, only two —
  `RNA-counts.mtx` and `cell-metadata.tsv` — are present there, and both carry the
  **same md5** as in 7924238.
- Record 7924238 lists all four files this study uses, each with a md5 matching the
  local copy. **Two of the four files are obtainable only from 7924238.**

Two descriptive discrepancies, recorded and not resolved because neither bears on
cohort identity, eligibility or any result:

1. The 7723202 page currently displays its version as **1.1**, whereas Appendix E
   of the frozen Gate 2b plan recorded that record as Version 1.0. The file list
   observed on 7723202 today — 7 files, no barcodes, no features — matches what the
   frozen plan described, so the substance the plan recorded is unchanged; only the
   displayed version label differs.
2. Appendix E recorded 7924238 as containing 15 files; the record page read today
   reports 12. The four files this study uses are present with matching checksums
   either way.

## Disposition

- The submitted manuscript's Data availability section keeps **Zenodo record
  7924238**. It is the record that actually holds all four files used.
- The frozen Gate 2 and Gate 2b plans are **not** edited. Their reference to
  7723202 is the true record of what was cited at plan fixation, and Appendix E of
  the Gate 2b plan already documents the correction to 7924238 as a factual
  correction that does not alter any analysis rule.
- A clarifying sentence is added to the Supplementary Methods transportability
  section so that a reader of the manuscript alone can see why two record numbers
  appear in the archive.

## Note on an earlier misstatement in this project's working notes

An earlier pass flagged the 7723202 / 7924238 pair as an unresolved discrepancy
needing an author decision. That flag was wrong: Appendix E of the frozen Gate 2b
plan had already identified and documented the correction, and Appendix F had
already recorded the four md5 checksums. The flag came from searching the
archive for the record number rather than reading the appendices, and the entry in
`SUBMISSION_CONTENT_BASELINE_2026-09-30.md` has been corrected accordingly. What
this check adds that the frozen plan did not already contain is the independent
recomputation of the four md5 values from the files on disk today, and the
confirmation against the live Zenodo checksums.
