# Project status — discovery analysis

## Decision

The discovery project is not blocked by Cell Ranger re-quantification. The working analysis uses the six-sample shared feature space (21,622 genes), defines CD4 T cells with `Cd3e+ / Cd8a- / Cd8b-`, and proceeds with the frozen G0A/G0B-2/G0C/P0 results.

## Statistical limitation

After the frozen QC/availability audit, animal-level inference has two usable animals per group. This is recorded as a limitation; no additional biological interpretation is added until the manuscript-level evidence is assembled.

## Frozen signature audit

The 107-gene signature is stable without Sham1: overlap 88.4%, with 100% directional concordance. This supports continuing the discovery line without waiting for a Cell Ranger rerun.

## Separate manuscript repair

Figure 1C can be corrected from the existing matrices using shared combination markers (`Cd3d/Cd3e`, `Cd8a/Cd8b`, `Cd40lg`, `Il7r`, `Themis`, `Lef1`). A unified Cell Ranger rerun is an optional manuscript/GEO repair and external validation, not a prerequisite for the current project.
