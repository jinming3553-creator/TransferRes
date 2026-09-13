from __future__ import annotations
import numpy as np


def select_features_training_only(
    X_train: np.ndarray,
    n_features: int | None,
    method: str = "variance",
) -> np.ndarray:
    p = X_train.shape[1]

    if n_features is None or n_features >= p:
        return np.arange(p, dtype=int)

    if n_features <= 0:
        raise ValueError("n_features must be positive or None.")

    if method == "variance":
        score = np.var(X_train, axis=0, ddof=1)
    elif method == "detection":
        score = np.sum(X_train != 0, axis=0)
    else:
        raise ValueError("feature_method must be 'variance' or 'detection'.")

    # Stable deterministic ranking.
    return np.argsort(-score, kind="stable")[:n_features]
