from dataclasses import dataclass
from typing import Any
import pandas as pd


@dataclass
class TransitionResult:
    parent_level: str
    child_level: str
    raw_incremental_gain: float
    partial_pooling_gain: float
    mean_weight: float
    fold_table: pd.DataFrame


@dataclass
class ResolutionResult:
    selected_level: str
    level_order: list[str]
    level_table: pd.DataFrame
    transition_table: pd.DataFrame
    fold_table: pd.DataFrame
    metadata: dict[str, Any]
