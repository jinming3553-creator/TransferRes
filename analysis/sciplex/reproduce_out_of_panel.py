#!/usr/bin/env python3
"""
Reproduce the sci-Plex 300-split out-of-panel validation used in Fig. 5 and Extended Data Fig. 6.

This script implements the frozen analysis contract used for the manuscript. When the optional frozen source-data CSV is supplied, all 300 split-level outputs are validated to floating-point precision.

Inputs
------
SciPlex_response_matrix_meta.csv
SciPlex_response_matrix_1000genes.npz
SciPlex_drug_annotations.csv
(optional) frozen source-data CSV exported from the manuscript workbook

Analysis contract
---------------------------
RNG:
    rng = np.random.default_rng(20260913)
    for split in range(300):
        perm = rng.permutation(sorted_drug_names)
        train = perm[:94]
        unseen = perm[94:]

Prediction / selection:
    Four resolutions:
      Shared
      Cell line
      Cell line × 2 dose bins (10/100 vs 1000/10000)
      Cell line × exact dose

    For each drug and held-out replicate, predictions are n_match-weighted context
    means from the opposite replicate. Raw squared Euclidean errors are summed.
    Training resolution = argmin raw held-out SSE.
    Unseen gain vs Shared = 1 - SSE_selected / SSE_Shared.

Downstream response-map representation (used for split-level reproducibility/MOA):
    1. L2-normalize each exact-dose 1000-gene response vector FIRST.
    2. Within each resolution context, take the UNWEIGHTED arithmetic mean of the
       normalized exact-dose vectors.
    3. For every drug pair, compute cosine similarity within each context and then
       take the UNWEIGHTED mean across contexts. This yields a pairwise drug
       response-similarity map for each replicate and resolution.

Reproducibility:
    Restrict the pairwise maps to the unseen drugs and compute Spearman correlation
    between rep1 and rep2 pairwise similarities.

MOA AUROC:
    Annotation layers: pathway, pathway_level_1, pathway_level_2, target.
    Within each unseen panel and annotation layer, retain only annotation classes
    represented by >=2 unseen drugs.
    Positive pair = same annotation class.
    Score = pairwise response-map cosine similarity defined above.
    Compute AUROC separately for each annotation layer and replicate, then average
    the 8 AUROCs (4 layers × 2 replicates).

Note:
    The separate whole-panel Fig.5 benchmark summaries used n_match-weighted raw
    context means before cosine similarity. The 300-split downstream validation
    used the normalized-first/unweighted map construction above; this distinction
    is required to reproduce the frozen split-level values exactly.
"""

import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score

LEVELS = [
    "Shared",
    "Cell line",
    "Cell line × 2 dose bins",
    "Cell line × exact dose",
]
ANNOTATIONS = ["pathway", "pathway_level_1", "pathway_level_2", "target"]
SEED = 20260913
N_SPLITS = 300

def context_cols(level):
    return {
        "Shared": [],
        "Cell line": ["cell_line"],
        "Cell line × 2 dose bins": ["cell_line", "dose_bin"],
        "Cell line × exact dose": ["cell_line", "dose"],
    }[level]

def load_inputs(meta_path, npz_path, annotation_path):
    meta = pd.read_csv(meta_path).copy()
    z = np.load(npz_path)
    R = np.asarray(z["response"], dtype=np.float64)
    unit_id = np.asarray(z["unit_id"])
    if len(meta) != R.shape[0] or not np.array_equal(unit_id, meta["unit_id"].to_numpy()):
        raise ValueError("Response matrix and metadata are not aligned.")
    meta["dose_bin"] = np.where(meta["dose"].astype(float) <= 100.0, "low", "high")
    ann = pd.read_csv(annotation_path)
    return meta, R, ann

def weighted_mean(X, w):
    return np.average(X, axis=0, weights=np.asarray(w, dtype=float))

def prediction_drug_sse(meta, R, drug, level):
    cols = context_cols(level)
    idx = np.flatnonzero(meta["perturbation"].to_numpy() == drug)
    total = 0.0
    for held in ("rep1", "rep2"):
        te = idx[meta.iloc[idx]["replicate"].to_numpy() == held]
        tr = idx[meta.iloc[idx]["replicate"].to_numpy() != held]

        if not cols:
            pred = {(): weighted_mean(R[tr], meta.iloc[tr]["n_match"])}
            keys = [()] * len(te)
        else:
            pred = {}
            for key, g in meta.iloc[tr].groupby(cols, sort=True):
                key = key if isinstance(key, tuple) else (key,)
                ii = g.index.to_numpy()
                pred[key] = weighted_mean(R[ii], g["n_match"])
            keys = []
            for _, row in meta.iloc[te].iterrows():
                keys.append(tuple(row[c] for c in cols))

        for ii, key in zip(te, keys):
            d = R[ii] - pred[key]
            total += float(d @ d)
    return total

def build_drug_sse(meta, R, drugs):
    out = pd.DataFrame(index=drugs, columns=LEVELS, dtype=float)
    for drug in drugs:
        for level in LEVELS:
            out.loc[drug, level] = prediction_drug_sse(meta, R, drug, level)
    return out

def normalized_context_map(meta, R, drug, rep, level):
    """Normalize exact vectors first, then unweighted-average within contexts."""
    cols = context_cols(level)
    g = meta[(meta["perturbation"] == drug) & (meta["replicate"] == rep)]
    groups = [g] if not cols else [x for _, x in g.groupby(cols, sort=True)]
    ctx = []
    for gg in groups:
        X = R[gg.index.to_numpy()]
        X = X / np.linalg.norm(X, axis=1, keepdims=True)
        ctx.append(X.mean(axis=0))
    return np.vstack(ctx)

def build_pair_maps(meta, R, drugs):
    pairs = [(drugs[i], drugs[j]) for i in range(len(drugs)) for j in range(i+1, len(drugs))]
    pair_index = {p:i for i,p in enumerate(pairs)}
    maps = {}

    cache = {}
    for level in LEVELS:
        for rep in ("rep1", "rep2"):
            cache[(level, rep)] = {
                d: normalized_context_map(meta, R, d, rep, level) for d in drugs
            }
            scores = np.empty(len(pairs), dtype=float)
            for k, (a,b) in enumerate(pairs):
                A = cache[(level, rep)][a]
                B = cache[(level, rep)][b]
                cos = np.sum(A*B, axis=1) / (
                    np.linalg.norm(A, axis=1) * np.linalg.norm(B, axis=1)
                )
                scores[k] = cos.mean()
            maps[(level, rep)] = scores
    return pairs, pair_index, maps

def unseen_pair_indices(unseen, pair_index):
    u = sorted(unseen)
    return np.asarray(
        [pair_index[(u[i],u[j])] for i in range(len(u)) for j in range(i+1,len(u))],
        dtype=int
    )

def reproducibility(unseen, level, pair_index, maps):
    ii = unseen_pair_indices(unseen, pair_index)
    return float(spearmanr(maps[(level,"rep1")][ii], maps[(level,"rep2")][ii]).statistic)

def moa_auc(unseen, level, ann, pair_index, maps):
    unseen = set(unseen)
    aa0 = ann[ann["perturbation"].isin(unseen)].copy()
    vals = []

    for rep in ("rep1","rep2"):
        score_map = maps[(level,rep)]
        for label in ANNOTATIONS:
            counts = aa0[label].value_counts()
            valid_labels = counts[counts >= 2].index
            aa = aa0[aa0[label].isin(valid_labels)].copy()
            aa = aa.sort_values("perturbation")
            ds = aa["perturbation"].tolist()
            lab = dict(zip(aa["perturbation"], aa[label].astype(str)))

            scores, y = [], []
            for i,a in enumerate(ds):
                for b in ds[i+1:]:
                    key = tuple(sorted((a,b)))
                    scores.append(score_map[pair_index[key]])
                    y.append(lab[a] == lab[b])
            vals.append(roc_auc_score(y, scores))

    return float(np.mean(vals))

def run(meta, R, ann):
    drugs = sorted(meta["perturbation"].unique())
    if len(drugs) != 188:
        raise ValueError(f"Expected 188 drugs, found {len(drugs)}")

    drug_sse = build_drug_sse(meta, R, drugs)
    _, pair_index, maps = build_pair_maps(meta, R, drugs)

    rng = np.random.default_rng(SEED)
    rows = []

    for split in range(N_SPLITS):
        perm = rng.permutation(np.asarray(drugs))
        train = perm[:94]
        unseen = perm[94:]

        train_sse = drug_sse.loc[train].sum(axis=0)
        selected = train_sse.idxmin()

        test_sse = drug_sse.loc[unseen].sum(axis=0)
        pred_gain = 1.0 - test_sse[selected] / test_sse["Shared"]

        r_sel = reproducibility(unseen, selected, pair_index, maps)
        r_coarse = reproducibility(unseen, "Shared", pair_index, maps)
        r_finest = reproducibility(unseen, "Cell line × exact dose", pair_index, maps)

        m_sel = moa_auc(unseen, selected, ann, pair_index, maps)
        m_coarse = moa_auc(unseen, "Shared", ann, pair_index, maps)
        m_finest = moa_auc(unseen, "Cell line × exact dose", ann, pair_index, maps)

        rows.append({
            "split": split,
            "selected_level": selected,
            "test_prediction_gain_vs_shared": float(pred_gain),
            "test_repro_selected": r_sel,
            "test_repro_coarse": r_coarse,
            "test_repro_finest": r_finest,
            "test_MOA_selected": m_sel,
            "test_MOA_coarse": m_coarse,
            "test_MOA_finest": m_finest,
            "repro_adv_vs_coarse": r_sel-r_coarse,
            "repro_adv_vs_finest": r_sel-r_finest,
            "MOA_adv_vs_coarse": m_sel-m_coarse,
            "MOA_adv_vs_finest": m_sel-m_finest,
        })
    return pd.DataFrame(rows)

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--meta", required=True)
    p.add_argument("--npz", required=True)
    p.add_argument("--annotations", required=True)
    p.add_argument("--out", default="SciPlex_drug_split_validation_RECONSTRUCTED.csv")
    p.add_argument("--frozen-csv", default=None,
                   help="Optional 300-row frozen CSV with the same 13 columns for exact verification.")
    args = p.parse_args()

    meta, R, ann = load_inputs(args.meta, args.npz, args.annotations)
    out = run(meta, R, ann)
    out.to_csv(args.out, index=False)

    print(out["selected_level"].value_counts())
    print("Median unseen prediction gain:", out["test_prediction_gain_vs_shared"].median())

    if args.frozen_csv:
        frozen = pd.read_csv(args.frozen_csv)
        if "index" in frozen.columns:
            frozen = frozen.drop(columns=["index"])
        if list(frozen.columns) != list(out.columns):
            raise ValueError("Frozen columns do not match reconstructed columns.")
        numeric = [c for c in out.columns if c not in ("selected_level",)]
        maxerr = max(
            np.max(np.abs(out[c].to_numpy(float)-frozen[c].to_numpy(float)))
            for c in numeric
        )
        same_levels = np.array_equal(out["selected_level"].to_numpy(),
                                     frozen["selected_level"].to_numpy())
        print("Max numeric absolute error:", maxerr)
        print("Selected levels identical:", same_levels)
        if maxerr > 1e-12 or not same_levels:
            raise RuntimeError("Frozen 300-row validation failed.")
        print("EXACT FROZEN VALIDATION PASSED.")

if __name__ == "__main__":
    main()
