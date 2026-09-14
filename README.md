# TransferRes 1.0.0
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22736901.svg)](https://doi.org/10.5281/zenodo.22736901)
**TransferRes — Transferable Perturbational Resolution**

TransferRes is a desktop analysis tool for estimating the hierarchy level at
which perturbational response patterns transfer across independent biological
units using leakage-aware held-out prediction.

## Status

Version **1.0.0** is the frozen first stable release. The core estimator is
frozen. Subsequent changes should be versioned and must not silently alter the
v1.0.0 statistical contract.

## Primary statistical contract

- The discrete transferable resolution is selected by **minimum raw held-out SSE**.
- Biological-unit holdout defines the transfer target.
- Nested feature selection is the recommended leakage-aware default.
- Partial-pooling gain is a diagnostic; it does not replace the discrete selection rule.
- Biological-unit and perturbation bootstraps are distinct and are labeled by resampling unit.
- Structural identifiability describes whether finer resolution is evaluable under the observed sampling geometry.
- Signal-injection calibration describes detectability, not a biological effect-size estimate.
- Zero transferable gain does not imply absence of finer biological structure.

## Windows quick start

For a built Windows release, open `TransferRes.exe`. For this source-and-builder
package, run `BUILD_FINAL_WINDOWS_RELEASE.bat` on Windows to create the portable
EXE release.

Inside the GUI:

1. Click **Load Sample Data** for the bundled example.
2. Keep the **Recommended** profile.
3. Validate the metadata mapping and hierarchy.
4. Run the analysis.
5. Review the **Results** tab and generated offline HTML report.

For new datasets, use **Generate Input Template** and map the biological unit,
perturbation, optional weight/depth column, and hierarchy from coarse to fine.

## Inputs

TransferRes uses:
- a metadata CSV containing biological-unit, perturbation and hierarchy labels;
- a response matrix (`.npy`, `.npz`, or supported CSV representation);
- an output folder.

See `QUICK_START.md` and the generated template for the expected layout.

## Outputs

Each analysis creates a timestamped result folder containing machine-readable
tables, diagnostic plots, logs, an analysis manifest, and an offline HTML report.

## Frozen acceptance gates

The v1.0.0 release lineage passed:
- Windows standalone build and packaged-EXE self-test;
- sample regression and feature-count overflow tests;
- paths containing spaces and non-ASCII characters;
- GUI launch, scrolling, Results page, HTML tables and plots;
- real-data manuscript reproduction for Whole-brain Glut Broad→Class:
  **PASS=14, FAIL=0, SKIP=0**.

The real-data gate reproduced NestedFS Top500/Top1000 fingerprints, final
geometry, and common-taxonomy permutation fingerprints.

## Real-data provenance

The manuscript-final Broad→Class regression uses held-out-source Top500/Top1000
feature lists frozen from the Cluster-analysis population. Re-ranking features
inside the Broad→Class subset is leakage-free but defines a different analysis
and is therefore not the v1.0.0 manuscript reproduction target.

The dedicated audit package is distributed separately as
`TransferRes_v1_0_0_REAL_DATA_REPRO_AUDIT_v3_1_IMPORT_HOTFIX.zip`.

## Citation and license

TransferRes v1.0.0 is archived on Zenodo:

**Liu, J.-M. (2026). TransferRes: Transferable Perturbational Resolution (Version 1.0.0). Zenodo. https://doi.org/10.5281/zenodo.22736901**

Source code, release assets, and documentation are available from the GitHub
repository: https://github.com/jinming3553-creator/TransferRes

See `CITATION.cff` for machine-readable citation metadata and `LICENSE` for
licensing information.
