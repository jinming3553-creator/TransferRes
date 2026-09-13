# TransferRes 1.0.0 release checklist

## Statistical contract
- [x] Minimum raw held-out SSE remains the discrete selection rule.
- [x] Partial pooling remains diagnostic only.
- [x] Nested feature selection default preserved.
- [x] Biological-unit and perturbation bootstrap remain distinct.
- [x] Common-taxonomy permutation contract preserved.
- [x] Structural-identifiability and calibration interpretation boundaries frozen.
- [x] Frozen manuscript feature-universe provenance documented.

## Desktop / packaging
- [x] Bundled sample data and template workflow.
- [x] Input validation.
- [x] Core / Recommended / Full profiles.
- [x] Progress reporting and cooperative cancellation.
- [x] Error logs and analysis manifest.
- [x] Results tab, PNG plots and offline HTML report.
- [x] Windows standalone EXE build PASS.
- [x] Packaged EXE self-test PASS.
- [x] Paths with spaces PASS.
- [x] Non-ASCII paths PASS.
- [x] Scrollable GUI PASS.
- [x] HTML tables and plots PASS.

## Reproducibility
- [x] Sample regression PASS.
- [x] Real-data Whole-brain Glut Broad→Class: PASS=14, FAIL=0, SKIP=0.
- [x] NestedFS Top500 and Top1000 frozen fingerprints reproduced.
- [x] Final-geometry fingerprints reproduced.
- [x] Common-taxonomy permutation fingerprints reproduced.
- [x] Frozen feature-selection provenance bundled in dedicated audit package.

## Public-release hygiene
- [x] LICENSE present.
- [x] CITATION.cff present.
- [x] VERSION.txt frozen at 1.0.0.
- [x] Release notes updated.
- [x] SHA256 checksums generated.
- [ ] Add public repository URL after repository creation.
- [ ] Add archival DOI after Zenodo/deposition.
- [ ] Optional: code-sign Windows executable if a signing certificate is available.
