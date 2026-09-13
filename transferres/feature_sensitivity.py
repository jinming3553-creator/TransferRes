from __future__ import annotations

from dataclasses import dataclass
import numpy as np
import pandas as pd

from .core import fit_resolution


@dataclass
class FeatureCountSensitivityResult:
    level_table: pd.DataFrame
    selection_table: pd.DataFrame
    contrast_table: pd.DataFrame
    interpretation_table: pd.DataFrame
    metadata: dict


def diagnose_feature_count_sensitivity(
    response,
    meta: pd.DataFrame,
    hierarchy: list[str],
    *,
    feature_counts,
    weight_col: str | None = None,
    feature_method: str = "variance",
    nested_feature_selection: bool = True,
    reference_level: str | None = None,
    comparator_level: str | None = None,
) -> FeatureCountSensitivityResult:
    """
    Evaluate whether transferable-resolution conclusions are stable to feature count.

    Each requested feature count is fit independently using the same frozen
    estimator. If nested_feature_selection=True, features are selected within
    each training biological-unit fold for every feature count.

    Parameters
    ----------
    response, meta, hierarchy
        Same inputs as fit_resolution.
    feature_counts
        Positive integers, e.g. [500, 1000, 2000].
        Counts greater than the number of available features are not allowed.
    weight_col
        Optional reliability/depth weight.
    feature_method
        Training-fold feature ranking method.
    nested_feature_selection
        Should remain True for leakage-free sensitivity analysis.
    reference_level
        Optional level whose selected-vs-reference stability should be tracked.
        Usually the primary resolution selected in the manuscript.
    comparator_level
        Optional specific neighboring alternative, e.g. coarse dose grouping.

    Returns
    -------
    FeatureCountSensitivityResult

    Interpretation
    --------------
    The module classifies feature-count robustness as:

      exact_selection_stability
          all feature counts select the same level;

      local_plateau_stability
          selection varies, but only between reference_level and comparator_level;

      feature_count_sensitive
          at least one feature count selects outside that local plateau.

    This diagnostic does not alter the primary feature count used by the main
    analysis.
    """
    X=np.asarray(response,dtype=float)
    p=X.shape[1]

    counts=sorted({int(x) for x in feature_counts})
    if not counts:
        raise ValueError("feature_counts must contain at least one positive integer.")
    if any(x<=0 for x in counts):
        raise ValueError("feature_counts must be positive.")
    if any(x>p for x in counts):
        raise ValueError(
            f"feature_counts cannot exceed the available feature dimension ({p})."
        )

    if reference_level is not None and reference_level not in hierarchy:
        raise ValueError("reference_level must be in hierarchy.")
    if comparator_level is not None and comparator_level not in hierarchy:
        raise ValueError("comparator_level must be in hierarchy.")

    level_rows=[]
    sel_rows=[]
    contrast_rows=[]

    for nf in counts:
        fit=fit_resolution(
            X,
            meta,
            hierarchy=hierarchy,
            weight_col=weight_col,
            n_features=nf,
            feature_method=feature_method,
            nested_feature_selection=nested_feature_selection,
        )

        for _,r in fit.level_table.iterrows():
            level_rows.append({
                "n_features":nf,
                "level":r["level"],
                "heldout_sse":float(r["heldout_sse"]),
                "gain_vs_coarsest":float(r["gain_vs_coarsest"]),
                "raw_incremental_gain_vs_parent":(
                    np.nan
                    if pd.isna(r["raw_incremental_gain_vs_parent"])
                    else float(r["raw_incremental_gain_vs_parent"])
                ),
            })

        sel_rows.append({
            "n_features":nf,
            "selected_level":fit.selected_level,
        })

        if reference_level is not None:
            tab=fit.level_table.set_index("level")
            ref_sse=float(tab.loc[reference_level,"heldout_sse"])

            if comparator_level is not None:
                cmp_sse=float(tab.loc[comparator_level,"heldout_sse"])
                denom=cmp_sse if cmp_sse>0 else np.nan
                contrast=(cmp_sse-ref_sse)/denom if np.isfinite(denom) else np.nan
                contrast_rows.append({
                    "n_features":nf,
                    "reference_level":reference_level,
                    "comparator_level":comparator_level,
                    "reference_gain_vs_comparator":contrast,
                    "reference_better":bool(ref_sse<cmp_sse),
                })

    level_table=pd.DataFrame(level_rows)
    selection_table=pd.DataFrame(sel_rows)
    contrast_table=pd.DataFrame(contrast_rows)

    selected_unique=list(pd.unique(selection_table["selected_level"]))

    if len(selected_unique)==1:
        status="exact_selection_stability"
    else:
        plateau=set(x for x in [reference_level,comparator_level] if x is not None)
        if plateau and set(selected_unique).issubset(plateau):
            status="local_plateau_stability"
        else:
            status="feature_count_sensitive"

    # Consistency of the full ranking across feature counts:
    # Spearman correlation between level gain vectors.
    ranking_rows=[]
    for i,a in enumerate(counts):
        va=(level_table[level_table.n_features.eq(a)]
            .set_index("level")
            .reindex(hierarchy)["gain_vs_coarsest"])
        for b in counts[i+1:]:
            vb=(level_table[level_table.n_features.eq(b)]
                .set_index("level")
                .reindex(hierarchy)["gain_vs_coarsest"])
            rho=float(va.corr(vb,method="spearman"))
            ranking_rows.append({
                "feature_count_a":a,
                "feature_count_b":b,
                "ranking_spearman":rho,
            })
    ranking_table=pd.DataFrame(ranking_rows)

    interpretation=pd.DataFrame([{
        "status":status,
        "n_feature_counts":len(counts),
        "all_counts_select_same_level":bool(len(selected_unique)==1),
        "selected_levels":" | ".join(map(str,selected_unique)),
        "minimum_pairwise_ranking_spearman":(
            np.nan if ranking_table.empty
            else float(ranking_table["ranking_spearman"].min())
        ),
    }])

    return FeatureCountSensitivityResult(
        level_table=level_table,
        selection_table=selection_table,
        contrast_table=contrast_table,
        interpretation_table=interpretation,
        metadata={
            "diagnostic_type":"feature-count sensitivity",
            "feature_counts":counts,
            "feature_method":feature_method,
            "nested_feature_selection":bool(nested_feature_selection),
            "primary_estimator_unchanged":True,
            "reference_level":reference_level,
            "comparator_level":comparator_level,
            "ranking_table":ranking_table,
        },
    )


def summarize_feature_count_results(
    level_results: pd.DataFrame,
    *,
    hierarchy: list[str],
    feature_count_col: str="n_features",
    level_col: str="resolution",
    gain_col: str="gain_vs_shared",
    reference_level: str | None=None,
    comparator_level: str | None=None,
) -> FeatureCountSensitivityResult:
    """
    Summarize already-computed feature-count sensitivity results.

    This helper exists for manuscript regression testing when fold-specific
    response matrices were produced externally and each feature count has
    already been evaluated with the frozen estimator.

    It does not re-fit models and therefore must only be used on results that
    were generated under the documented estimator contract.
    """
    req={feature_count_col,level_col,gain_col}
    missing=req-set(level_results.columns)
    if missing:
        raise ValueError(f"Missing columns: {sorted(missing)}")

    z=level_results.copy()
    counts=sorted(pd.unique(z[feature_count_col]).astype(int))

    rows=[]
    sel=[]
    ctr=[]

    for nf in counts:
        g=z[z[feature_count_col].eq(nf)].copy()
        g=g.set_index(level_col).reindex(hierarchy)

        if g[gain_col].isna().any():
            raise ValueError(f"Missing hierarchy level at feature count {nf}.")

        for lev,val in g[gain_col].items():
            rows.append({
                "n_features":nf,
                "level":lev,
                "gain_vs_coarsest":float(val),
            })

        best=str(g[gain_col].idxmax())
        sel.append({"n_features":nf,"selected_level":best})

        if reference_level is not None and comparator_level is not None:
            diff=float(
                g.loc[reference_level,gain_col]
                - g.loc[comparator_level,gain_col]
            )
            ctr.append({
                "n_features":nf,
                "reference_level":reference_level,
                "comparator_level":comparator_level,
                "reference_gain_minus_comparator_gain":diff,
                "reference_better":bool(diff>0),
            })

    level_table=pd.DataFrame(rows)
    selection_table=pd.DataFrame(sel)
    contrast_table=pd.DataFrame(ctr)

    unique=list(pd.unique(selection_table["selected_level"]))

    if len(unique)==1:
        status="exact_selection_stability"
    else:
        plateau=set(x for x in [reference_level,comparator_level] if x is not None)
        status=(
            "local_plateau_stability"
            if plateau and set(unique).issubset(plateau)
            else "feature_count_sensitive"
        )

    ranking_rows=[]
    for i,a in enumerate(counts):
        va=(level_table[level_table.n_features.eq(a)]
            .set_index("level").reindex(hierarchy)["gain_vs_coarsest"])
        for b in counts[i+1:]:
            vb=(level_table[level_table.n_features.eq(b)]
                .set_index("level").reindex(hierarchy)["gain_vs_coarsest"])
            ranking_rows.append({
                "feature_count_a":a,
                "feature_count_b":b,
                "ranking_spearman":float(va.corr(vb,method="spearman")),
            })
    ranking_table=pd.DataFrame(ranking_rows)

    interpretation=pd.DataFrame([{
        "status":status,
        "n_feature_counts":len(counts),
        "all_counts_select_same_level":bool(len(unique)==1),
        "selected_levels":" | ".join(map(str,unique)),
        "minimum_pairwise_ranking_spearman":(
            np.nan if ranking_table.empty
            else float(ranking_table.ranking_spearman.min())
        ),
    }])

    return FeatureCountSensitivityResult(
        level_table=level_table,
        selection_table=selection_table,
        contrast_table=contrast_table,
        interpretation_table=interpretation,
        metadata={
            "diagnostic_type":"feature-count sensitivity summary",
            "feature_counts":counts,
            "primary_estimator_unchanged":True,
            "reference_level":reference_level,
            "comparator_level":comparator_level,
            "ranking_table":ranking_table,
            "precomputed_results":True,
        },
    )
