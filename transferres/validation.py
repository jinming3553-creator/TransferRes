from __future__ import annotations
import pandas as pd


REQUIRED_BASE_COLUMNS = {"biological_unit", "perturbation"}


def validate_input(
    meta: pd.DataFrame,
    response,
    hierarchy: list[str],
    weight_col: str | None = None,
) -> None:
    if len(meta) != len(response):
        raise ValueError("meta and response must contain the same number of rows.")

    missing = REQUIRED_BASE_COLUMNS - set(meta.columns)
    if missing:
        raise ValueError(f"Missing required metadata columns: {sorted(missing)}")

    if len(hierarchy) < 1:
        raise ValueError("hierarchy must contain at least one biological resolution column.")

    missing_h = [c for c in hierarchy if c not in meta.columns]
    if missing_h:
        raise ValueError(f"Hierarchy columns not found in metadata: {missing_h}")

    if weight_col is not None and weight_col not in meta.columns:
        raise ValueError(f"weight_col '{weight_col}' not found in metadata.")

    if meta["biological_unit"].isna().any():
        raise ValueError("biological_unit contains missing values.")

    if meta["perturbation"].isna().any():
        raise ValueError("perturbation contains missing values.")

    if meta["biological_unit"].nunique() < 2:
        raise ValueError(
            "At least two independent biological units are required. "
            "Cells, wells, or technical partitions must not be substituted "
            "for biological replicates unless they are genuinely independent."
        )

    if meta.duplicated(
        subset=["biological_unit", "perturbation"] + hierarchy
    ).any():
        # Duplicates are not forbidden, but raw cell-level data should not silently
        # enter a unit-level estimator.
        raise ValueError(
            "Duplicate biological_unit × perturbation × hierarchy rows detected. "
            "Aggregate raw cells to one response vector per analysis unit before fitting."
        )
