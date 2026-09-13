from __future__ import annotations

from dataclasses import dataclass
import numpy as np
import pandas as pd


@dataclass
class IdentifiabilityResult:
    parent_level: str
    child_level: str
    summary_table: pd.DataFrame
    support_frontier: pd.DataFrame
    perturbation_table: pd.DataFrame
    stratum_table: pd.DataFrame
    metadata: dict


def diagnose_identifiability(
    meta: pd.DataFrame,
    *,
    parent_col: str,
    child_col: str,
    biological_unit_col: str = "biological_unit",
    perturbation_col: str = "perturbation",
    max_training_sources: int = 5,
) -> IdentifiabilityResult:
    """
    Quantify structural support for an adjacent hierarchy transition.

    A fine-response stratum is defined as:

        perturbation × parent × child

    For each stratum, the function counts how many independent biological
    units contribute observations.

    The leave-one-unit-out support frontier then asks:

        after holding out one biological unit, what fraction of fine strata
        (or observed rows) retain at least k TRAINING biological units?

    Therefore support at k training units requires:

        total biological units in the stratum >= k + 1

    This is a structural identifiability diagnostic, NOT a hypothesis test.

    It must not be interpreted as evidence that a finer taxonomy has or lacks
    biological meaning. Signal detectability also depends on effect size,
    noise, dimensionality, and estimator behavior; those require separate
    calibration.

    Returns
    -------
    IdentifiabilityResult
        summary_table:
            global counts and replication summaries.
        support_frontier:
            support as a function of required training biological units k.
        perturbation_table:
            perturbation-specific support summaries.
        stratum_table:
            one row per perturbation × parent × child stratum.
    """
    if max_training_sources < 1:
        raise ValueError("max_training_sources must be at least 1.")

    required={
        parent_col,child_col,biological_unit_col,perturbation_col
    }
    missing=required-set(meta.columns)
    if missing:
        raise ValueError(f"Missing metadata columns: {sorted(missing)}")

    m=meta.reset_index(drop=True).copy()

    if m[biological_unit_col].isna().any():
        raise ValueError("biological-unit column contains missing values.")
    if m[perturbation_col].isna().any():
        raise ValueError("perturbation column contains missing values.")

    # One row can represent an aggregated analysis unit. If accidental
    # duplicates exist, count biological units, not rows, for replication.
    strata_cols=[perturbation_col,parent_col,child_col]

    strata=(
        m.groupby(strata_cols,dropna=False)
        .agg(
            n_biological_units=(biological_unit_col,"nunique"),
            n_rows=(biological_unit_col,"size"),
        )
        .reset_index()
    )

    # Attach stratum support back to observed rows for row-weighted coverage.
    row_support=m.merge(
        strata[strata_cols+["n_biological_units"]],
        on=strata_cols,
        how="left",
        validate="many_to_one",
    )

    # Global summary.
    counts=strata["n_biological_units"].to_numpy(int)

    summary=pd.DataFrame([{
        "n_biological_units_total":int(m[biological_unit_col].nunique()),
        "n_perturbations":int(m[perturbation_col].nunique()),
        "n_parent_labels":int(m[parent_col].nunique(dropna=False)),
        "n_child_labels":int(m[child_col].nunique(dropna=False)),
        "n_fine_strata":int(len(strata)),
        "median_sources_per_stratum":float(np.median(counts)) if len(counts) else np.nan,
        "q25_sources_per_stratum":float(np.quantile(counts,.25)) if len(counts) else np.nan,
        "q75_sources_per_stratum":float(np.quantile(counts,.75)) if len(counts) else np.nan,
        "min_sources_per_stratum":int(np.min(counts)) if len(counts) else 0,
        "max_sources_per_stratum":int(np.max(counts)) if len(counts) else 0,
        "structurally_unevaluable_strata_LOO":int(np.sum(counts < 2)),
        "structurally_unevaluable_strata_fraction_LOO":float(np.mean(counts < 2)) if len(counts) else np.nan,
    }])

    # Support frontier.
    frontier_rows=[]
    for k in range(1,max_training_sources+1):
        threshold=k+1
        strata_supported=strata["n_biological_units"]>=threshold
        rows_supported=row_support["n_biological_units"]>=threshold

        # Perturbation has at least one fine stratum meeting the support level.
        per_pert=(
            strata.assign(_supported=strata_supported)
            .groupby(perturbation_col,dropna=False)["_supported"]
            .any()
        )

        frontier_rows.append({
            "required_training_biological_units_after_holdout":k,
            "required_total_biological_units_per_stratum":threshold,
            "supported_strata_count":int(strata_supported.sum()),
            "supported_strata_fraction":float(strata_supported.mean()) if len(strata) else np.nan,
            "supported_observed_row_fraction":float(rows_supported.mean()) if len(row_support) else np.nan,
            "perturbations_with_any_supported_child_count":int(per_pert.sum()),
            "perturbations_with_any_supported_child_fraction":float(per_pert.mean()) if len(per_pert) else np.nan,
        })

    frontier=pd.DataFrame(frontier_rows)

    # Perturbation-specific summaries.
    pert_rows=[]
    for pert,g in strata.groupby(perturbation_col,dropna=False):
        rec={
            perturbation_col:pert,
            "n_fine_strata":int(len(g)),
            "median_sources_per_stratum":float(g["n_biological_units"].median()),
            "max_sources_per_stratum":int(g["n_biological_units"].max()),
        }
        for k in range(1,max_training_sources+1):
            threshold=k+1
            rec[f"strata_fraction_train_ge_{k}"]=float(
                (g["n_biological_units"]>=threshold).mean()
            )
            rec[f"any_child_train_ge_{k}"]=bool(
                (g["n_biological_units"]>=threshold).any()
            )
        pert_rows.append(rec)

    perturbation_table=pd.DataFrame(pert_rows)

    return IdentifiabilityResult(
        parent_level=parent_col,
        child_level=child_col,
        summary_table=summary,
        support_frontier=frontier,
        perturbation_table=perturbation_table,
        stratum_table=strata,
        metadata={
            "diagnostic_type":"structural replication support frontier",
            "stratum_definition":"perturbation × parent × child",
            "holdout_definition":"leave one biological unit out",
            "not_a_hypothesis_test":True,
            "does_not_establish_absence_of_biology":True,
            "requires_positive_control_calibration_for_signal_detectability":True,
        },
    )
