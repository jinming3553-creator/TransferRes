from __future__ import annotations

from dataclasses import dataclass
import numpy as np
import pandas as pd


@dataclass
class SignalCalibrationResult:
    calibration_table: pd.DataFrame
    metadata: dict


def calibrate_signal_detectability(
    table: pd.DataFrame,
    *,
    biological_unit_col: str,
    A_col: str,
    B_col: str,
    C_col: str,
    fractions,
    n_boot: int = 0,
    seed: int | None = 1,
    alpha: float = 0.05,
) -> SignalCalibrationResult:
    """
    Calibrate fine-level signal detectability on the observed sampling geometry.

    A/B/C geometry
    --------------
    For each held-out biological unit:

        A = ||y - mu_parent||^2
        B = (y - mu_parent)'(mu_child - mu_parent)
        C = ||mu_child - mu_parent||^2

    A fraction f of the observed child-direction vector is injected into the
    held-out response:

        (y - mu_parent)' = (y - mu_parent) + f (mu_child - mu_parent)

    This changes the sufficient statistics exactly to:

        A_f = A + 2 f B + f^2 C
        B_f = B + f C
        C_f = C

    The biological-unit structure, hierarchy coverage, and sampling geometry
    are otherwise unchanged.

    For every f, the training-only shrinkage weight for each held-out unit is
    re-estimated from all OTHER biological units.

    Parameters
    ----------
    table
        One or more A/B/C rows per held-out biological unit.
    biological_unit_col
        Held-out biological unit identifier.
    A_col, B_col, C_col
        Quadratic sufficient-statistic columns.
    fractions
        Iterable of non-negative injected signal fractions.
    n_boot
        Optional biological-unit bootstrap replicates. Zero disables bootstrap.
    seed
        RNG seed for bootstrap.
    alpha
        Percentile interval alpha.

    Interpretation
    --------------
    This is a positive-control / detectability calibration, not a biological
    effect-size estimate. The injected fraction scales the empirically observed
    child-direction geometry.

    A near-zero observed gain should only be interpreted as informative when
    signals of scientifically relevant magnitude are recoverable under the
    observed sampling geometry.
    """
    required={biological_unit_col,A_col,B_col,C_col}
    missing=required-set(table.columns)
    if missing:
        raise ValueError(f"Missing columns: {sorted(missing)}")

    fr=np.asarray(list(fractions),dtype=float)
    if fr.ndim != 1 or len(fr)==0:
        raise ValueError("fractions must be a non-empty one-dimensional iterable.")
    if np.any(~np.isfinite(fr)) or np.any(fr < 0):
        raise ValueError("fractions must contain finite non-negative values.")
    if n_boot < 0:
        raise ValueError("n_boot must be >= 0.")
    if not (0 < alpha < 1):
        raise ValueError("alpha must lie between 0 and 1.")

    # Aggregate to one A/B/C row per biological unit.
    s=(table.groupby(biological_unit_col,sort=False)
       .agg(A=(A_col,"sum"),B=(B_col,"sum"),C=(C_col,"sum"))
       .reset_index())

    n_units=len(s)
    if n_units < 2:
        raise ValueError("At least two biological units are required.")

    A0=s["A"].to_numpy(float)
    B0=s["B"].to_numpy(float)
    C0=s["C"].to_numpy(float)

    def evaluate(A,B,C):
        At=A.sum()
        Bt=B.sum()
        Ct=C.sum()

        ws=np.zeros(len(A),dtype=float)
        gains=np.zeros(len(A),dtype=float)

        for i in range(len(A)):
            trainB=Bt-B[i]
            trainC=Ct-C[i]
            w=0.0 if trainC<=0 else float(np.clip(trainB/trainC,0.0,1.0))
            ws[i]=w
            gains[i]=2*w*B[i]-w*w*C[i]

        relative=np.nan if At<=0 else float(gains.sum()/At)

        return {
            "mean_weight":float(ws.mean()),
            "relative_gain":relative,
            "positive_gain_unit_fraction":float(np.mean(gains>0)),
        }

    rows=[]
    rng=np.random.default_rng(seed)

    # Use the same bootstrap index matrix for every injected fraction so
    # calibration curves are paired across f.
    if n_boot>0:
        boot_idx=rng.integers(0,n_units,size=(n_boot,n_units))
    else:
        boot_idx=None

    for f in fr:
        A=A0 + 2.0*f*B0 + (f*f)*C0
        B=B0 + f*C0
        C=C0.copy()

        point=evaluate(A,B,C)

        rec={
            "injected_signal_fraction":float(f),
            "mean_weight":point["mean_weight"],
            "relative_gain":point["relative_gain"],
            "positive_gain_unit_fraction":point["positive_gain_unit_fraction"],
        }

        if n_boot>0:
            bg=np.empty(n_boot,dtype=float)
            bw=np.empty(n_boot,dtype=float)

            for b in range(n_boot):
                ix=boot_idx[b]
                ev=evaluate(A[ix],B[ix],C[ix])
                bg[b]=ev["relative_gain"]
                bw[b]=ev["mean_weight"]

            rec.update({
                "bootstrap_mean_gain":float(np.nanmean(bg)),
                "gain_interval_low":float(np.nanquantile(bg,alpha/2)),
                "gain_interval_high":float(np.nanquantile(bg,1-alpha/2)),
                "bootstrap_mean_weight":float(np.nanmean(bw)),
            })

        rows.append(rec)

    cal=pd.DataFrame(rows).sort_values("injected_signal_fraction").reset_index(drop=True)

    # Descriptive threshold: first tested fraction whose bootstrap percentile
    # lower bound is > 0. This is grid-dependent and is NOT a formal detection limit.
    if n_boot>0:
        hit=cal[cal["gain_interval_low"]>0]
        first_bootstrap_positive=(
            None if hit.empty else float(hit.iloc[0]["injected_signal_fraction"])
        )
    else:
        first_bootstrap_positive=None

    return SignalCalibrationResult(
        calibration_table=cal,
        metadata={
            "calibration_type":"data-matched child-direction signal injection",
            "injection_formula":"A_f=A+2fB+f^2C; B_f=B+fC; C_f=C",
            "n_biological_units":int(n_units),
            "n_boot":int(n_boot),
            "seed":seed,
            "alpha":float(alpha),
            "first_tested_fraction_with_bootstrap_lower_bound_above_zero":first_bootstrap_positive,
            "threshold_is_grid_dependent":True,
            "not_a_biological_effect_size_estimate":True,
            "not_a_formal_power_calculation":True,
        },
    )
