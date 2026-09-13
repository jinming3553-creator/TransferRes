from __future__ import annotations

import numpy as np
import pandas as pd

from .features import select_features_training_only
from .types import ResolutionResult
from .validation import validate_input


def _weighted_mean(X: np.ndarray, w: np.ndarray) -> np.ndarray:
    w = np.asarray(w, dtype=float)
    if np.any(w < 0):
        raise ValueError("Weights must be non-negative.")
    s = w.sum()
    if s <= 0:
        raise ValueError("Weights must contain at least one positive value.")
    return np.average(X, axis=0, weights=w)


def _weights(meta, mask, weight_col):
    if weight_col is None:
        return np.ones(int(mask.sum()), dtype=float)
    return meta.loc[mask, weight_col].to_numpy(float)


def _raw_level_sse(
    X,
    meta,
    train_mask,
    test_mask,
    level_col,
    weight_col,
    feature_idx,
):
    Xf = X[:, feature_idx]
    sse = 0.0
    n_eval = 0

    for idx in np.flatnonzero(test_mask):
        row = meta.iloc[idx]
        mask = (
            train_mask
            & meta["perturbation"].eq(row["perturbation"]).to_numpy()
            & meta[level_col].eq(row[level_col]).to_numpy()
        )
        if mask.sum() == 0:
            continue

        mu = _weighted_mean(Xf[mask], _weights(meta, mask, weight_col))
        e = Xf[idx] - mu
        sse += float(e @ e)
        n_eval += 1

    return sse, n_eval


def _transition_fold_stats(
    X,
    meta,
    train_mask,
    test_mask,
    parent_col,
    child_col,
    weight_col,
    feature_idx,
):
    Xf = X[:, feature_idx]

    A=B=C=raw_gain=0.0
    n_eval=0

    for idx in np.flatnonzero(test_mask):
        row=meta.iloc[idx]
        same_pert=meta["perturbation"].eq(row["perturbation"]).to_numpy()

        coarse_mask=(
            train_mask & same_pert
            & meta[parent_col].eq(row[parent_col]).to_numpy()
        )
        fine_mask=(
            coarse_mask
            & meta[child_col].eq(row[child_col]).to_numpy()
        )

        if coarse_mask.sum()==0 or fine_mask.sum()==0:
            continue

        mu_c=_weighted_mean(Xf[coarse_mask],_weights(meta,coarse_mask,weight_col))
        mu_f=_weighted_mean(Xf[fine_mask],_weights(meta,fine_mask,weight_col))
        y=Xf[idx]

        a=y-mu_c
        c=mu_f-mu_c
        Ai=float(a@a)
        Bi=float(a@c)
        Ci=float(c@c)

        A+=Ai; B+=Bi; C+=Ci
        raw_gain += Ai - float((y-mu_f)@(y-mu_f))
        n_eval+=1

    return A,B,C,raw_gain,n_eval


def fit_resolution(
    response,
    meta: pd.DataFrame,
    hierarchy: list[str],
    *,
    weight_col: str | None=None,
    n_features: int | None=None,
    feature_method: str="variance",
    nested_feature_selection: bool=True,
) -> ResolutionResult:
    """
    Estimate transferable perturbational resolution.

    Two distinct outputs are intentionally separated:

    1. Discrete resolution selection:
       L* is chosen by RAW leave-one-biological-unit-out SSE at each
       candidate hierarchy level.

    2. Partial-pooling diagnostic:
       For each adjacent parent->child transition, a training-only
       shrinkage weight estimates whether the child deviation contains
       a transferable component. These shrunken gains are NOT accumulated
       to choose L*.

    This distinction is essential: a fully specified fine model may
    overfit (negative raw transfer) while a small, heavily shrunk fine
    component remains reproducibly transferable.
    """
    X=np.asarray(response,dtype=float)
    meta=meta.reset_index(drop=True).copy()
    validate_input(meta,X,hierarchy,weight_col)

    units=list(pd.unique(meta["biological_unit"]))
    preprocessing_note="nested" if nested_feature_selection else "transductive"

    # --------------------------------------------------------
    # 1) Raw level-wise held-out SSE used ONLY for L* selection.
    # --------------------------------------------------------
    raw_fold_rows=[]
    feature_cache={}

    for held in units:
        train_mask=meta["biological_unit"].ne(held).to_numpy()
        test_mask=~train_mask

        if nested_feature_selection:
            feat_idx=select_features_training_only(
                X[train_mask], n_features, feature_method
            )
        else:
            feat_idx=select_features_training_only(
                X, n_features, feature_method
            )
        feature_cache[held]=feat_idx

        for lev in hierarchy:
            sse,n_eval=_raw_level_sse(
                X,meta,train_mask,test_mask,lev,weight_col,feat_idx
            )
            raw_fold_rows.append({
                "heldout_unit":held,
                "level":lev,
                "sse":sse,
                "n_eval":n_eval,
                "n_features_used":len(feat_idx),
            })

    raw_fold=pd.DataFrame(raw_fold_rows)
    level_rows=[]
    coarse_sse=raw_fold.loc[raw_fold.level==hierarchy[0],"sse"].sum()

    previous_sse=None
    for lev in hierarchy:
        ss=raw_fold.loc[raw_fold.level==lev,"sse"].sum()
        gain_vs_coarse=np.nan if coarse_sse<=0 else (coarse_sse-ss)/coarse_sse
        if previous_sse is None or previous_sse<=0:
            inc=np.nan
        else:
            inc=(previous_sse-ss)/previous_sse
        level_rows.append({
            "level":lev,
            "heldout_sse":ss,
            "gain_vs_coarsest":gain_vs_coarse,
            "raw_incremental_gain_vs_parent":inc,
            "n_biological_units":len(units),
        })
        previous_sse=ss

    level_table=pd.DataFrame(level_rows)
    selected_level=level_table.loc[level_table.heldout_sse.idxmin(),"level"]

    # --------------------------------------------------------
    # 2) Adjacent transition partial-pooling diagnostics.
    # --------------------------------------------------------
    transition_rows=[]
    transition_folds=[]

    for parent_col,child_col in zip(hierarchy[:-1],hierarchy[1:]):
        suff=[]
        for held in units:
            train_mask=meta["biological_unit"].ne(held).to_numpy()
            test_mask=~train_mask
            feat_idx=feature_cache[held]

            A,B,C,raw_gain,n_eval=_transition_fold_stats(
                X,meta,train_mask,test_mask,parent_col,child_col,
                weight_col,feat_idx
            )
            suff.append({
                "heldout_unit":held,
                "A":A,"B":B,"C":C,
                "raw_gain":raw_gain,
                "n_eval":n_eval,
                "n_features_used":len(feat_idx),
            })

        fs=pd.DataFrame(suff)
        Btot=fs.B.sum()
        Ctot=fs.C.sum()

        ws=[]; gs=[]
        for _,r in fs.iterrows():
            train_B=Btot-r.B
            train_C=Ctot-r.C
            w=0.0 if train_C<=0 else float(np.clip(train_B/train_C,0,1))
            g=2*w*r.B-w*w*r.C
            ws.append(w); gs.append(g)

        fs["weight"]=ws
        fs["partial_pooling_gain_absolute"]=gs
        fs["parent_level"]=parent_col
        fs["child_level"]=child_col

        Aall=fs.A.sum()
        pp=np.nan if Aall<=0 else fs.partial_pooling_gain_absolute.sum()/Aall
        raw=np.nan if Aall<=0 else fs.raw_gain.sum()/Aall

        transition_rows.append({
            "parent_level":parent_col,
            "child_level":child_col,
            "raw_incremental_gain":raw,
            "partial_pooling_gain":pp,
            "mean_weight":fs.weight.mean(),
            "n_biological_units":len(units),
            "n_evaluable_rows":int(fs.n_eval.sum()),
        })
        transition_folds.append(fs)

    transition_table=pd.DataFrame(transition_rows)
    transition_fold_table=(
        pd.concat(transition_folds,ignore_index=True)
        if transition_folds else pd.DataFrame()
    )

    # Combine raw-level and transition fold diagnostics with a table tag.
    rf=raw_fold.copy()
    rf["table_type"]="raw_level"
    tf=transition_fold_table.copy()
    if len(tf):
        tf["table_type"]="transition"
    fold_table=pd.concat([rf,tf],ignore_index=True,sort=False)

    return ResolutionResult(
        selected_level=selected_level,
        level_order=hierarchy,
        level_table=level_table,
        transition_table=transition_table,
        fold_table=fold_table,
        metadata={
            "n_biological_units":len(units),
            "n_features_requested":n_features,
            "feature_method":feature_method,
            "feature_selection":preprocessing_note,
            "primary_metric":"euclidean_sse",
            "resolution_selection_rule":"minimum raw held-out SSE",
            "partial_pooling_role":"diagnostic only; not used to choose selected_level",
            "conditional_on_matched_control":True,
        },
    )
