# sci-Plex manuscript reproduction

These scripts reproduce the sci-Plex analyses used in the TransferRes manuscript.

## Scripts

1. `build_response_matrix.py` — constructs the 24-h A549/K562/MCF7 response matrices and the single global variance ranking used as nested 500/1000/2000-gene feature sets.
2. `reproduce_resolution_analysis.py` — reproduces reciprocal-replicate held-out SSE and primary resolution selection.
3. `reproduce_featurecount_bootstrap.py` — reproduces 500/1000/2000-gene robustness and the 20,000-sample drug-level bootstrap.
4. `reproduce_out_of_panel.py` — reproduces the 300 random 94/94 drug-split analysis, unseen-drug prediction, response-map reproducibility and MOA recovery used in Fig. 5 and Extended Data Fig. 6.

## Analysis contract

Candidate resolutions are `Shared`, `Cell line`, `Cell line × 2 dose bins`, and `Cell line × exact dose`. The two dose bins are 10/100 and 1000/10000. Resolution selection uses `n_match`-weighted training centroids and raw held-out squared Euclidean error. The drug bootstrap and 300-split analysis both use NumPy `default_rng(20260913)` according to the manuscript Methods.

The out-of-panel script intentionally distinguishes the downstream split-level response-map construction (L2-normalize exact-dose vectors first, then unweighted context averaging) from the raw-SSE prediction layer.
