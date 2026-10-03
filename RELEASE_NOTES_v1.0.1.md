# v1.0.1 — 2026-10-03

Bibliographic and submission-metadata correction. No analysis definition, input
file, result value, figure or figure source table changed; every checksum in
`CHECKSUMS.sha256` for those files is identical to v1.0.0.

- `README.md`: the three commands under "Reproducing the reported numbers"
  pointed at `code/analysis/` and `code/decomposition_reference/`, paths from an
  earlier layout that was abandoned because it broke the scripts' own path
  resolution. They are corrected to `scripts/` and `decomposition_reference/`,
  which is what the Layout section of the same file already described. All three
  were re-run from the archive root after the correction and pass.
- `README.md`: a Requirements section was added, naming the Python and package
  versions the results were produced and re-verified under, to meet the
  journal's requirement that a code deposition carry installation instructions
  alongside its manual and usage example.
- `README.md`, `CITATION.cff`, `.zenodo.json`: the archive DOI is updated from
  10.5281/zenodo.23116737 (v1.0.0) to 10.5281/zenodo.23118493 (v1.0.1), and the
  recorded version is 1.0.1.
- `CHECKSUMS.sha256`: regenerated. v1.0.0 listed a macOS `.DS_Store` file that is
  excluded by `.gitignore` and so never reached the public repository, which made
  the manifest one entry longer than the published tree.

The accompanying manuscript cites this version DOI. The concept DOI
10.5281/zenodo.23116736 always resolves to the newest version.
