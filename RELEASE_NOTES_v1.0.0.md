# TransferRes 1.0.0 release notes

TransferRes 1.0.0 is the first stable release and the release-engineering freeze
of the accepted estimator lineage.

## Frozen method

The selected transferable perturbational resolution is the tested hierarchy
level with minimum raw held-out SSE across biological-unit holdouts. Nested
feature selection is the recommended leakage-aware default. Partial pooling,
bootstrap uncertainty, structural identifiability, permutation and
signal-injection calibration remain separate diagnostics with conservative
interpretation boundaries.

## Final acceptance

### Windows / desktop acceptance
PASS:
- standalone PyInstaller build;
- packaged EXE self-test;
- frozen sample regression;
- feature-count overflow handling;
- paths containing spaces;
- non-ASCII paths;
- GUI launch and scrollable Setup & Run page;
- Recommended sample workflow;
- Results page;
- offline HTML tables and plots.

### Real-data manuscript reproduction
Target: **Whole-brain Glut Broad→Class**

Final gate: **PASS=14, FAIL=0, SKIP=0**

Passed fingerprints:
- NestedFS Top500: raw incremental gain, partial-pooling gain, mean weight;
- NestedFS Top1000: raw incremental gain, partial-pooling gain, mean weight;
- final geometry: raw incremental gain, partial-pooling gain, mean weight;
- common-taxonomy permutation: greater-tail p, two-sided p, null median,
  null q025 and null q975.

The final audit reproduced the manuscript-final fingerprint.

## Provenance correction resolved before freeze

During release auditing, a six-target NestedFS discrepancy was traced to feature
universe provenance rather than estimator mathematics. The frozen manuscript
analysis used training-only held-out-source feature lists defined on the frozen
Cluster-analysis population. A later audit implementation had re-ranked
features inside the Broad→Class subset, which is a different analysis.

Differential audits independently confirmed response construction,
coarse/fine predictors, A/B/C sufficient-statistic formulas and biological-unit
membership before the exact frozen feature universe was restored.

No tolerance was relaxed to obtain the final PASS.

## Freeze rule

The v1.0.0 estimator and interpretation contract are frozen. Any future change
that can alter numerical output or inferential meaning requires a new version.
