# Reproducibility contract — TransferRes 1.0.0

This document records the numerical and interpretive boundaries frozen for
TransferRes 1.0.0.

## Selection rule

The discrete transferable resolution is the tested hierarchy level with minimum
raw held-out SSE across biological-unit holdouts.

Partial-pooling gain is a diagnostic and must not be substituted for the
discrete selection rule.

## Leakage control

When feature selection is requested, feature ranking must be performed using
training information appropriate to the frozen analysis provenance. The
manuscript-final Whole-brain Glut Broad→Class target uses held-out-source
Top500/Top1000 feature lists frozen from the Cluster-analysis population.

## Uncertainty and diagnostics

Biological-unit bootstrap and perturbation bootstrap answer different questions
and must be labeled by their resampling unit. Structural identifiability
describes evaluability under the sampling geometry. Signal-injection fractions
are calibration parameters, not biological effect-size estimates.

## v1.0.0 release gates

The Windows acceptance suite passed standalone build, packaged self-test, sample
regression, overflow handling, path robustness, GUI workflow and HTML rendering.

The real-data manuscript gate passed **14/14** frozen numerical targets with
**FAIL=0, SKIP=0**, including NestedFS Top500/Top1000, final geometry and
common-taxonomy permutation.

## Change control

Do not silently modify estimator behavior, feature-universe provenance,
resampling definitions, selection rules, or interpretation language under the
v1.0.0 version number. Any output-changing change requires a new version and a
fresh acceptance record.
