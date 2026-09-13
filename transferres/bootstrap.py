from __future__ import annotations

from dataclasses import dataclass
import numpy as np
import pandas as pd

from .types import ResolutionResult


@dataclass
class BootstrapResult:
    raw_level_table: pd.DataFrame
    transition_table: pd.DataFrame
    selection_frequency: pd.DataFrame
    metadata: dict


def _quantile_interval(x, alpha):
    x=np.asarray(x,dtype=float)
    return (
        float(np.nanquantile(x, alpha/2)),
        float(np.nanquantile(x, 1-alpha/2)),
    )


def bootstrap_biological_units(
    result: ResolutionResult,
    *,
    n_boot: int = 2000,
    seed: int | None = 1,
    alpha: float = 0.05,
) -> BootstrapResult:
    """
    Biological-unit bootstrap for a fitted transferable-resolution result.

    The bootstrap unit is the same independent biological unit used in the
    outer cross-validation loop (e.g. mouse/source, donor, replicate).

    Raw resolution uncertainty
    --------------------------
    Held-out SSE contributions are resampled by biological unit and summed.
    The selected level in each bootstrap replicate is the level with minimum
    resampled raw held-out SSE.

    Partial-pooling transition uncertainty
    --------------------------------------
    Fold-level A/B/C sufficient statistics are resampled by biological unit.
    For each bootstrap pseudo-unit j, the shrinkage weight is re-estimated
    from all OTHER bootstrap pseudo-units:

        w_j = clip((sum B - B_j) / (sum C - C_j), 0, 1)

    This preserves the training-only definition of the transition estimator
    inside each bootstrap replicate.

    Important
    ---------
    This is a biological-unit bootstrap interval only when the resampled
    units are genuinely independent biological replicates.

    With very few biological units, especially n < 5, percentile intervals
    should be interpreted as stability diagnostics rather than strong
    population-level confidence intervals.
    """
    if n_boot <= 0:
        raise ValueError("n_boot must be positive.")
    if not (0 < alpha < 1):
        raise ValueError("alpha must lie between 0 and 1.")

    ft=result.fold_table.copy()

    raw=ft[ft["table_type"].eq("raw_level")].copy()
    tr=ft[ft["table_type"].eq("transition")].copy()

    if raw.empty:
        raise ValueError("Result does not contain raw-level fold data.")

    units=list(pd.unique(raw["heldout_unit"]))
    n_units=len(units)
    if n_units < 2:
        raise ValueError("At least two biological units are required.")

    unit_to_i={u:i for i,u in enumerate(units)}
    levels=result.level_order

    # matrix: unit x level raw SSE
    raw_mat=np.full((n_units,len(levels)),np.nan,dtype=float)
    for li,lev in enumerate(levels):
        z=(raw[raw["level"].eq(lev)]
           .groupby("heldout_unit",sort=False)["sse"].sum())
        for u,v in z.items():
            raw_mat[unit_to_i[u],li]=float(v)

    if np.isnan(raw_mat).any():
        raise ValueError("Missing raw SSE contribution for at least one unit × level.")

    rng=np.random.default_rng(seed)
    boot_idx=rng.integers(0,n_units,size=(n_boot,n_units))

    # ---------- raw levels ----------
    boot_sse=raw_mat[boot_idx,:].sum(axis=1)
    base=boot_sse[:,[0]]
    boot_gain=(base-boot_sse)/base
    best=np.argmin(boot_sse,axis=1)

    raw_rows=[]
    for li,lev in enumerate(levels):
        lo,hi=_quantile_interval(boot_gain[:,li],alpha)
        raw_rows.append({
            "level":lev,
            "point_gain_vs_coarsest":float(
                result.level_table.loc[
                    result.level_table["level"].eq(lev),
                    "gain_vs_coarsest"
                ].iloc[0]
            ),
            "bootstrap_mean_gain_vs_coarsest":float(np.mean(boot_gain[:,li])),
            "bootstrap_median_gain_vs_coarsest":float(np.median(boot_gain[:,li])),
            "interval_low":lo,
            "interval_high":hi,
        })
    raw_table=pd.DataFrame(raw_rows)

    sf=pd.DataFrame({
        "level":levels,
        "selection_frequency":[float(np.mean(best==li)) for li in range(len(levels))]
    })

    # ---------- transitions ----------
    transition_rows=[]
    if not tr.empty:
        for _,point in result.transition_table.iterrows():
            parent=point["parent_level"]
            child=point["child_level"]
            z=tr[
                tr["parent_level"].eq(parent)
                & tr["child_level"].eq(child)
            ].copy()

            # Aggregate if a unit appears in multiple transition rows.
            z=(z.groupby("heldout_unit",sort=False)
               .agg(A=("A","sum"),B=("B","sum"),C=("C","sum"),
                    raw_gain=("raw_gain","sum"))
               .reindex(units))

            if z[["A","B","C","raw_gain"]].isna().any().any():
                raise ValueError(
                    f"Missing transition sufficient statistics for {parent}->{child}."
                )

            M=z[["A","B","C","raw_gain"]].to_numpy(float)
            sampled=M[boot_idx]  # boot x pseudo-unit x 4

            A=sampled[:,:,0]
            B=sampled[:,:,1]
            C=sampled[:,:,2]
            RG=sampled[:,:,3]

            At=A.sum(axis=1)
            Bt=B.sum(axis=1)
            Ct=C.sum(axis=1)

            # Re-estimate training-only weight for every bootstrap pseudo-fold.
            trainB=Bt[:,None]-B
            trainC=Ct[:,None]-C
            W=np.zeros_like(trainB)
            ok=trainC>0
            W[ok]=np.clip(trainB[ok]/trainC[ok],0,1)

            partial_abs=(2*W*B-W*W*C).sum(axis=1)
            partial=np.divide(
                partial_abs,At,
                out=np.full_like(partial_abs,np.nan),
                where=At>0
            )
            raw_gain=np.divide(
                RG.sum(axis=1),At,
                out=np.full_like(At,np.nan),
                where=At>0
            )

            plo,phi=_quantile_interval(partial,alpha)
            rlo,rhi=_quantile_interval(raw_gain,alpha)

            transition_rows.append({
                "parent_level":parent,
                "child_level":child,
                "point_raw_incremental_gain":float(point["raw_incremental_gain"]),
                "raw_bootstrap_mean":float(np.nanmean(raw_gain)),
                "raw_interval_low":rlo,
                "raw_interval_high":rhi,
                "point_partial_pooling_gain":float(point["partial_pooling_gain"]),
                "partial_bootstrap_mean":float(np.nanmean(partial)),
                "partial_interval_low":plo,
                "partial_interval_high":phi,
                "point_mean_weight":float(point["mean_weight"]),
                "bootstrap_mean_weight":float(np.nanmean(W)),
            })

    transition_table=pd.DataFrame(transition_rows)

    return BootstrapResult(
        raw_level_table=raw_table,
        transition_table=transition_table,
        selection_frequency=sf,
        metadata={
            "bootstrap_unit":"biological_unit",
            "n_biological_units":n_units,
            "n_boot":int(n_boot),
            "seed":seed,
            "alpha":float(alpha),
            "interval_type":"percentile biological-unit bootstrap",
            "low_unit_count_warning":bool(n_units < 5),
            "population_ci_language_recommended":bool(n_units >= 5),
        },
    )
