
from __future__ import annotations
import pandas as pd


def build_interpretation(fit, extra: dict) -> list[str]:
    """
    Conservative, software-generated interpretation statements.

    These statements intentionally avoid:
    - claiming that zero transferable gain means no biology;
    - calling partial-pooling gain the selected resolution;
    - treating perturbation bootstrap as a biological-population CI;
    - claiming exact cross-metric agreement when only an intermediate plateau is robust.
    """
    lines=[]

    lines.append(
        f"The raw held-out SSE criterion selected '{fit.selected_level}' "
        "as the transferable perturbational resolution among the tested levels."
    )

    if len(fit.transition_table):
        for _,r in fit.transition_table.iterrows():
            raw=float(r["raw_incremental_gain"])
            pp=float(r["partial_pooling_gain"])
            w=float(r["mean_weight"])
            parent=r["parent_level"]
            child=r["child_level"]

            if raw < 0 and pp > 0:
                lines.append(
                    f"{parent}→{child}: the full finer model reduced held-out transfer "
                    f"(raw gain {100*raw:.3f}%), while a small shrunken finer-level "
                    f"component remained transferable (partial gain {100*pp:.3f}%; "
                    f"mean weight {w:.3f})."
                )
            elif raw < 0 and abs(pp) < 1e-4:
                lines.append(
                    f"{parent}→{child}: the full finer model reduced held-out transfer "
                    f"(raw gain {100*raw:.3f}%), and the shrunken transferable component "
                    "was approximately zero at the tested resolution."
                )
            elif raw > 0:
                lines.append(
                    f"{parent}→{child}: the finer level improved raw held-out prediction "
                    f"by {100*raw:.3f}% relative to its parent."
                )
            else:
                lines.append(
                    f"{parent}→{child}: no positive raw held-out gain was observed. "
                    "This does not imply that the finer annotation lacks biological structure."
                )

    if extra.get("bootstrap") is not None:
        b=extra["bootstrap"]
        if b.metadata.get("low_unit_count_warning"):
            lines.append(
                "Biological-unit count is low; bootstrap intervals should be treated as "
                "stability diagnostics rather than strong population-level confidence intervals."
            )
        else:
            lines.append(
                "Biological-unit bootstrap was used to quantify uncertainty across independent "
                "biological units."
            )

    if extra.get("ident") is not None:
        i=extra["ident"]
        s=i.summary_table.iloc[0]
        frac=float(s["structurally_unevaluable_strata_fraction_LOO"])
        if frac > 0:
            lines.append(
                f"Structural support is incomplete: {100*frac:.1f}% of fine strata are "
                "unevaluable under leave-one-biological-unit-out because they occur in only one unit. "
                "Near-zero fine-level gain should therefore be interpreted together with this coverage limit."
            )

    if extra.get("metric") is not None:
        tab=extra["metric"].interpretation_table
        unresolved=tab[tab["status"].eq("statistically_unresolved")]
        shifted=tab[tab["status"].eq("metric_sensitive_shift")]
        if len(shifted):
            lines.append(
                "At least one alternative response-space metric favored a different resolution "
                "with a perturbation-bootstrap contrast excluding zero; the exact optimum is metric-sensitive."
            )
        elif len(unresolved):
            lines.append(
                "Alternative metrics changed the point-estimate optimum but did not resolve the "
                "difference from the primary level; the broader resolution regime is more robust than the exact boundary."
            )
        else:
            lines.append(
                "The tested response-space metrics agreed on the same best resolution."
            )

    if extra.get("feature_sensitivity") is not None:
        status=extra["feature_sensitivity"].interpretation_table.iloc[0]["status"]
        if status=="exact_selection_stability":
            lines.append(
                "The selected resolution was unchanged across the tested feature counts."
            )
        elif status=="local_plateau_stability":
            lines.append(
                "Feature count changed the exact optimum only within the prespecified local resolution plateau."
            )
        else:
            lines.append(
                "The selected resolution was sensitive to feature count; interpretation should be qualified accordingly."
            )

    if extra.get("permutation") is not None:
        p=extra["permutation"]
        if p.p_partial_greater < 0.05:
            lines.append(
                f"The shrunken finer-level component exceeded the common-taxonomy permutation null "
                f"(one-sided P={p.p_partial_greater:.4g}). This supports shared taxonomic alignment "
                "across perturbations, not necessarily a large full-model predictive gain."
            )
        else:
            lines.append(
                f"The shrunken finer-level component did not exceed the common-taxonomy permutation null "
                f"(one-sided P={p.p_partial_greater:.4g})."
            )

    if extra.get("calibration") is not None:
        lines.append(
            "Signal-injection calibration reports detectability under the observed sampling geometry; "
            "injected fractions are calibration parameters, not biological effect-size estimates."
        )

    return lines
