# Changelog

## 1.0.0

First stable release.

### Stable features
- Point-and-click Windows GUI.
- Sample data and input-template workflow.
- Core / Recommended / Full analysis profiles.
- Leakage-aware nested feature selection.
- Biological-unit holdout resolution selection.
- Partial-pooling diagnostics.
- Biological-unit bootstrap.
- Structural identifiability diagnostics.
- Optional metric sensitivity, feature-count sensitivity, common-taxonomy
  permutation, and signal-injection calibration.
- Software-assisted interpretation with explicit claim boundaries.
- HTML report with semantic tables, plots, and selected-resolution highlighting.
- Standalone Windows EXE build and acceptance harness.

### Release freeze
The primary estimator and inference contract are unchanged from the accepted
v0.7.3 release candidate.

## 0.7.0-rc1

Packaging release candidate.

### Added
- Standalone Windows EXE acceptance self-test.
- Automated build + acceptance workflow.
- Portable release packaging.
- SHA256 checksum generation.
- Version manifest.
- CITATION.cff.
- Reproducibility documentation.
- Release manifest and release checklist.
- Conservative software-assisted interpretation.
- Input template generator.
- Core / Recommended / Full run profiles.
- Progress reporting, cooperative cancellation, and automatic error logs.

### Statistical contract
No change to the frozen estimator.

- Discrete resolution selection remains minimum raw held-out SSE.
- Partial-pooling gain remains diagnostic only.
- Nested feature selection remains the default.
- Biological-unit bootstrap and perturbation bootstrap remain distinct.
- Common-taxonomy permutation retains one mapping per biological-unit × parent stratum, shared across perturbations.
- Zero transferable gain is not interpreted as absence of biological structure.
