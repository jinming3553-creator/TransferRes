from __future__ import annotations

from dataclasses import dataclass
import numpy as np
import pandas as pd

from .core import fit_resolution
from .permutation import common_taxonomy_relabel
from .fast_permutation import (
    prepare_fast_common_taxonomy,
    permute_child_codes,
    sufficient_from_codes,
    nested_partial_gain,
)


@dataclass
class PermutationResult:
    parent_level: str
    child_level: str
    observed_raw_incremental_gain: float
    observed_partial_pooling_gain: float
    null_table: pd.DataFrame
    p_raw_greater: float
    p_raw_two_sided: float
    p_partial_greater: float
    p_partial_two_sided: float
    metadata: dict



def common_taxonomy_permutation_test(
    response,
    meta: pd.DataFrame,
    *,
    parent_col: str,
    child_col: str,
    weight_col: str | None = None,
    n_features: int | None = None,
    feature_method: str = "variance",
    nested_feature_selection: bool = True,
    n_permutations: int = 1000,
    seed: int | None = 1,
    engine: str = "auto",
) -> PermutationResult:
    """
    Hierarchy-preserving common-taxonomy permutation test.

    `engine="auto"` uses a compiled sufficient-statistics fast path when
    `n_features is None`. This path is mathematically equivalent to the full
    response-vector estimator but avoids re-fitting the entire model for every
    permutation.

    When fold-specific feature selection is requested (`n_features` not None),
    the function falls back to the reference estimator because each permutation
    must respect the frozen feature-selection contract.
    """
    if n_permutations <= 0:
        raise ValueError("n_permutations must be positive.")
    if parent_col == child_col:
        raise ValueError("parent_col and child_col must differ.")
    if engine not in {"auto","fast","reference"}:
        raise ValueError("engine must be 'auto', 'fast', or 'reference'.")

    X=np.asarray(response,dtype=float)
    m=meta.reset_index(drop=True).copy()

    observed=fit_resolution(
        X,m,
        hierarchy=[parent_col,child_col],
        weight_col=weight_col,
        n_features=n_features,
        feature_method=feature_method,
        nested_feature_selection=nested_feature_selection,
    )
    tr=observed.transition_table.iloc[0]
    obs_raw=float(tr["raw_incremental_gain"])
    obs_partial=float(tr["partial_pooling_gain"])

    use_fast=(engine=="fast") or (engine=="auto" and n_features is None)
    if engine=="fast" and n_features is not None:
        raise ValueError(
            "Fast permutation engine currently requires n_features=None. "
            "Use engine='reference' for fold-specific feature selection."
        )

    rng=np.random.default_rng(seed)
    rows=[]

    if use_fast:
        prepared=prepare_fast_common_taxonomy(
            X,m,parent_col,child_col,weight_col
        )

        # Warm compiled kernel using observed codes, then verify exact agreement
        # with the reference estimator before any permutation is accepted.
        A,B,C=sufficient_from_codes(prepared,prepared["child"])
        fast_raw,fast_partial,fast_w=nested_partial_gain(A,B,C)

        if not (
            np.isclose(fast_raw,obs_raw,rtol=1e-9,atol=1e-10)
            and np.isclose(fast_partial,obs_partial,rtol=1e-9,atol=1e-10)
        ):
            raise RuntimeError(
                "Fast permutation engine failed equivalence check against "
                "the reference estimator."
            )

        for i in range(n_permutations):
            mapped=permute_child_codes(prepared,rng)
            A,B,C=sufficient_from_codes(prepared,mapped)
            raw,partial,w=nested_partial_gain(A,B,C)
            rows.append({
                "permutation":i,
                "raw_incremental_gain":raw,
                "partial_pooling_gain":partial,
                "mean_weight":w,
            })
        engine_used="fast_sufficient_statistics"

    else:
        for i in range(n_permutations):
            mp=m.copy()
            mp[child_col]=common_taxonomy_relabel(
                m,parent_col=parent_col,child_col=child_col,rng=rng
            ).to_numpy()
            fit=fit_resolution(
                X,mp,
                hierarchy=[parent_col,child_col],
                weight_col=weight_col,
                n_features=n_features,
                feature_method=feature_method,
                nested_feature_selection=nested_feature_selection,
            )
            r=fit.transition_table.iloc[0]
            rows.append({
                "permutation":i,
                "raw_incremental_gain":float(r["raw_incremental_gain"]),
                "partial_pooling_gain":float(r["partial_pooling_gain"]),
                "mean_weight":float(r["mean_weight"]),
            })
        engine_used="reference_refit"

    null=pd.DataFrame(rows)
    nr=null["raw_incremental_gain"].to_numpy(float)
    npart=null["partial_pooling_gain"].to_numpy(float)

    p_raw_g=(1+int(np.sum(nr>=obs_raw)))/(n_permutations+1)
    p_raw_2=(1+int(np.sum(np.abs(nr)>=abs(obs_raw))))/(n_permutations+1)
    p_part_g=(1+int(np.sum(npart>=obs_partial)))/(n_permutations+1)
    p_part_2=(1+int(np.sum(np.abs(npart)>=abs(obs_partial))))/(n_permutations+1)

    return PermutationResult(
        parent_level=parent_col,
        child_level=child_col,
        observed_raw_incremental_gain=obs_raw,
        observed_partial_pooling_gain=obs_partial,
        null_table=null,
        p_raw_greater=float(p_raw_g),
        p_raw_two_sided=float(p_raw_2),
        p_partial_greater=float(p_part_g),
        p_partial_two_sided=float(p_part_2),
        metadata={
            "null_type":"common taxonomy relabeling",
            "mapping_scope":"one mapping per biological_unit × parent stratum, shared across perturbations",
            "n_permutations":int(n_permutations),
            "seed":seed,
            "plus_one_correction":True,
            "nested_feature_selection":bool(nested_feature_selection),
            "engine":engine_used,
        },
    )
