from __future__ import annotations
import numpy as np
import pandas as pd


def common_taxonomy_relabel(
    meta: pd.DataFrame,
    parent_col: str,
    child_col: str,
    rng: np.random.Generator,
) -> pd.Series:
    """
    One relabeling per biological_unit × parent stratum,
    shared across all perturbations.

    This is intentionally NOT a perturbation-wise permutation.
    """
    out = meta[child_col].astype(object).copy()

    for (_, _), g in meta.groupby(
        ["biological_unit", parent_col], sort=False, dropna=False
    ):
        labels = pd.Index(g[child_col].drop_duplicates())
        if len(labels) <= 1:
            continue

        permuted = rng.permutation(labels.to_numpy())
        mapping = dict(zip(labels.tolist(), permuted.tolist()))
        out.loc[g.index] = g[child_col].map(mapping)

    return out
