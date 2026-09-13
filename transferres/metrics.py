from __future__ import annotations

from dataclasses import dataclass
import numpy as np
import pandas as pd

from .features import select_features_training_only
from .validation import validate_input


@dataclass
class MetricSensitivityResult:
    score_table: pd.DataFrame
    bootstrap_table: pd.DataFrame
    primary_contrast_table: pd.DataFrame
    interpretation_table: pd.DataFrame
    row_scores: pd.DataFrame
    metadata: dict


def _weighted_mean(X, w):
    w=np.asarray(w,dtype=float)
    if np.any(w<0):
        raise ValueError("Weights must be non-negative.")
    if w.sum()<=0:
        raise ValueError("Weights must contain at least one positive value.")
    return np.average(X,axis=0,weights=w)


def _score_euclidean_normalized(y,pred):
    den=float(y@y)
    if den<=0:
        return np.nan
    e=y-pred
    return 1.0-float(e@e)/den


def _score_cosine(y,pred):
    den=float(np.linalg.norm(y)*np.linalg.norm(pred))
    if den<=0:
        return np.nan
    return float(y@pred)/den


def _score_pearson(y,pred):
    sy=float(np.std(y))
    sp=float(np.std(pred))
    if sy<=0 or sp<=0:
        return np.nan
    return float(np.corrcoef(y,pred)[0,1])


_METRICS={
    "Euclidean":_score_euclidean_normalized,
    "Cosine":_score_cosine,
    "Pearson":_score_pearson,
}


def diagnose_metric_sensitivity(
    response,
    meta: pd.DataFrame,
    hierarchy: list[str],
    *,
    weight_col: str | None = None,
    n_features: int | None = None,
    feature_method: str = "variance",
    nested_feature_selection: bool = True,
    primary_level: str | None = None,
    n_boot: int = 5000,
    seed: int | None = 1,
    alpha: float = 0.05,
) -> MetricSensitivityResult:
    """
    Compare raw held-out resolution rankings under Euclidean, cosine and Pearson metrics.

    This is a sensitivity diagnostic. It does NOT replace the package primary
    resolution selector, which remains raw held-out SSE.

    Scoring
    -------
    For each held-out biological-unit row and hierarchy level, the predictor is
    estimated from training biological units with the same perturbation and
    level label.

    Metrics:
      Euclidean = 1 - ||y-pred||^2 / ||y||^2
      Cosine    = cosine(y,pred)
      Pearson   = corr(y,pred)

    Scores are averaged within perturbation first, so each perturbation receives
    equal weight in the reported metric means.

    Perturbation bootstrap
    ----------------------
    Bootstrap intervals resample perturbations, not biological units. Therefore
    these intervals describe robustness to perturbation-panel composition and
    must be labeled perturbation-bootstrap intervals, not population CIs over
    biological replicates.

    Interpretation
    --------------
    For each metric:
      - exact_agreement:
          metric-best level == primary SSE-selected level.
      - statistically_unresolved:
          a different metric-best level has a primary-minus-best bootstrap
          interval that includes zero.
      - metric_sensitive_shift:
          a different metric-best level is favored and the primary-minus-best
          interval excludes zero.

    No claim is made that all metrics must choose exactly the same resolution.
    """
    X=np.asarray(response,dtype=float)
    m=meta.reset_index(drop=True).copy()
    validate_input(m,X,hierarchy,weight_col)

    if n_boot<=0:
        raise ValueError("n_boot must be positive.")
    if not (0<alpha<1):
        raise ValueError("alpha must lie between 0 and 1.")

    if primary_level is not None and primary_level not in hierarchy:
        raise ValueError("primary_level must be one of the hierarchy columns.")

    units=list(pd.unique(m["biological_unit"]))
    row_records=[]

    for held in units:
        train_mask=m["biological_unit"].ne(held).to_numpy()
        test_idx=np.flatnonzero(~train_mask)

        if nested_feature_selection:
            feat=select_features_training_only(
                X[train_mask],n_features,feature_method
            )
        else:
            feat=select_features_training_only(
                X,n_features,feature_method
            )

        Xf=X[:,feat]

        for idx in test_idx:
            r=m.iloc[idx]
            y=Xf[idx]
            same_pert=m["perturbation"].eq(r["perturbation"]).to_numpy()

            for lev in hierarchy:
                q=(
                    train_mask
                    & same_pert
                    & m[lev].eq(r[lev]).to_numpy()
                )
                if q.sum()==0:
                    continue

                if weight_col is None:
                    w=np.ones(q.sum(),dtype=float)
                else:
                    w=m.loc[q,weight_col].to_numpy(float)

                pred=_weighted_mean(Xf[q],w)

                for metric,fn in _METRICS.items():
                    row_records.append({
                        "heldout_unit":held,
                        "perturbation":r["perturbation"],
                        "level":lev,
                        "metric":metric,
                        "score":fn(y,pred),
                    })

    row_scores=pd.DataFrame(row_records)

    # Equal weight per perturbation.
    pert=(
        row_scores.groupby(["perturbation","metric","level"],as_index=False)
        .agg(score=("score","mean"))
    )

    score_table=(
        pert.groupby(["metric","level"],as_index=False)
        .agg(mean_score=("score","mean"))
    )
    score_table["rank"]=(
        score_table.groupby("metric")["mean_score"]
        .rank(ascending=False,method="min")
    )

    # If caller did not supply the primary SSE-selected level, use the
    # Euclidean-normalized best only as a diagnostic reference. In manuscript
    # workflows callers should pass fit_resolution(...).selected_level.
    if primary_level is None:
        primary_level=(
            score_table[score_table.metric.eq("Euclidean")]
            .sort_values(["rank","level"])
            .iloc[0]["level"]
        )
        primary_level_source="Euclidean diagnostic best"
    else:
        primary_level_source="external raw-SSE selected level"

    perturbations=list(pd.unique(pert["perturbation"]))
    rng=np.random.default_rng(seed)
    boot_idx=rng.integers(0,len(perturbations),size=(n_boot,len(perturbations)))

    # Matrices per metric: perturbation x hierarchy level.
    bootstrap_rows=[]
    contrast_rows=[]
    interp_rows=[]

    for metric in _METRICS:
        z=(
            pert[pert.metric.eq(metric)]
            .pivot(index="perturbation",columns="level",values="score")
            .reindex(index=perturbations,columns=hierarchy)
        )

        if z.isna().any().any():
            raise ValueError(
                f"Missing score for at least one perturbation × level under {metric}."
            )

        M=z.to_numpy(float)
        BM=np.mean(M[boot_idx,:],axis=1)  # n_boot x levels
        wins=np.argmax(BM,axis=1)

        point=M.mean(axis=0)
        best_i=int(np.argmax(point))
        metric_best=hierarchy[best_i]
        primary_i=hierarchy.index(primary_level)

        for j,lev in enumerate(hierarchy):
            lo=float(np.quantile(BM[:,j],alpha/2))
            hi=float(np.quantile(BM[:,j],1-alpha/2))
            bootstrap_rows.append({
                "metric":metric,
                "level":lev,
                "mean_score":float(point[j]),
                "interval_low":lo,
                "interval_high":hi,
                "P_best":float(np.mean(wins==j)),
            })

        d=BM[:,primary_i]-BM[:,best_i]
        dlo=float(np.quantile(d,alpha/2))
        dhi=float(np.quantile(d,1-alpha/2))
        dpoint=float(point[primary_i]-point[best_i])

        if metric_best==primary_level:
            status="exact_agreement"
        elif dlo<=0<=dhi:
            status="statistically_unresolved"
        else:
            status="metric_sensitive_shift"

        contrast_rows.append({
            "metric":metric,
            "primary_level":primary_level,
            "metric_best_level":metric_best,
            "primary_minus_metric_best":dpoint,
            "interval_low":dlo,
            "interval_high":dhi,
            "P_primary_better":float(np.mean(d>0)),
        })

        interp_rows.append({
            "metric":metric,
            "primary_level":primary_level,
            "metric_best_level":metric_best,
            "status":status,
            "exact_level_agreement":bool(metric_best==primary_level),
            "contrast_interval_includes_zero":bool(dlo<=0<=dhi),
        })

    return MetricSensitivityResult(
        score_table=score_table,
        bootstrap_table=pd.DataFrame(bootstrap_rows),
        primary_contrast_table=pd.DataFrame(contrast_rows),
        interpretation_table=pd.DataFrame(interp_rows),
        row_scores=row_scores,
        metadata={
            "diagnostic_type":"response-space metric sensitivity",
            "primary_resolution_rule_unchanged":"minimum raw held-out SSE",
            "primary_level":primary_level,
            "primary_level_source":primary_level_source,
            "metrics":["Euclidean","Cosine","Pearson"],
            "bootstrap_unit":"perturbation",
            "interval_type":"perturbation-bootstrap interval",
            "n_boot":int(n_boot),
            "seed":seed,
            "alpha":float(alpha),
            "not_a_biological_unit_population_CI":True,
            "does_not_require_exact_metric_agreement":True,
        },
    )
