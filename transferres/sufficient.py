from __future__ import annotations

import numpy as np
import pandas as pd


def summarize_transition_sufficient_stats(
    table: pd.DataFrame,
    *,
    biological_unit_col: str,
    A_col: str,
    B_col: str,
    C_col: str,
    raw_gain_col: str | None = None,
) -> dict:
    """
    Reconstruct a nested partial-pooling transition from fold-level A/B/C
    sufficient statistics.

    Parameters
    ----------
    table
        One or more rows per held-out biological unit. Rows belonging to the
        same biological unit are summed before the nested shrinkage calculation.
    biological_unit_col
        Column identifying the held-out biological unit.
    A_col, B_col, C_col
        Quadratic sufficient statistics:
            A = ||y - mu_parent||^2
            B = (y - mu_parent)'(mu_child - mu_parent)
            C = ||mu_child - mu_parent||^2
    raw_gain_col
        Optional column containing raw parent-SSE minus raw child-SSE.
        If omitted, raw gain is reconstructed exactly as 2B-C.

    Returns
    -------
    dict
        raw_incremental_gain
        partial_pooling_gain
        mean_weight
        min_weight
        max_weight
        n_biological_units
        A_total
        fold_table

    Notes
    -----
    For each held-out biological unit u, the shrinkage weight is learned only
    from OTHER biological units:

        w_u = clip((sum B - B_u) / (sum C - C_u), 0, 1)

    The partial-pooling gain for u is:

        2 w_u B_u - w_u^2 C_u

    The returned raw and partial-pooling gains are distinct estimands.
    """
    required={biological_unit_col,A_col,B_col,C_col}
    if raw_gain_col is not None:
        required.add(raw_gain_col)
    missing=required-set(table.columns)
    if missing:
        raise ValueError(f"Missing sufficient-statistic columns: {sorted(missing)}")

    t=table.copy()

    if raw_gain_col is None:
        t["_raw_gain_reconstructed"]=2.0*t[B_col]-t[C_col]
        rg="_raw_gain_reconstructed"
    else:
        rg=raw_gain_col

    fold=(t.groupby(biological_unit_col,sort=False)
          .agg(
              A=(A_col,"sum"),
              B=(B_col,"sum"),
              C=(C_col,"sum"),
              raw_gain=(rg,"sum"),
          ).reset_index())

    if len(fold)<2:
        raise ValueError("At least two biological units are required.")

    Atot=float(fold.A.sum())
    Btot=float(fold.B.sum())
    Ctot=float(fold.C.sum())

    ws=[]
    pg=[]
    for r in fold.itertuples(index=False):
        train_B=Btot-float(r.B)
        train_C=Ctot-float(r.C)
        w=0.0 if train_C<=0 else float(np.clip(train_B/train_C,0.0,1.0))
        gain=2.0*w*float(r.B)-w*w*float(r.C)
        ws.append(w)
        pg.append(gain)

    fold["weight"]=ws
    fold["partial_pooling_gain_absolute"]=pg

    raw=np.nan if Atot<=0 else float(fold.raw_gain.sum()/Atot)
    partial=np.nan if Atot<=0 else float(fold.partial_pooling_gain_absolute.sum()/Atot)

    return {
        "raw_incremental_gain":raw,
        "partial_pooling_gain":partial,
        "mean_weight":float(np.mean(ws)),
        "min_weight":float(np.min(ws)),
        "max_weight":float(np.max(ws)),
        "n_biological_units":int(len(fold)),
        "A_total":Atot,
        "fold_table":fold,
    }
