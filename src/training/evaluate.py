"""
Evaluation utilities for the paralog retention GAT.

Functions here operate on already-computed predictions (e.g. the contents
of results/predicted_pi.csv) and produce summary metrics. They are pure
in the sense documented for src/data/features.py: numpy/pandas in,
numpy/pandas out, no I/O of their own.
"""

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score


def compute_auc_roc(y_true, y_pred):
    """ROC AUC for binary regulatory inheritance prediction.

    Returns NaN if y_true has only one class (cannot compute ROC).
    """
    y_true = np.asarray(y_true)
    y_pred = np.asarray(y_pred)
    if len(np.unique(y_true)) < 2:
        return float("nan")
    return float(roc_auc_score(y_true, y_pred))


def compute_family_correlation(predictions_df):
    """Pearson correlation between predicted pi per family and observed
    fraction of retained pairs per family.

    Per-family pi is the mean of `predicted_pi` across the family's
    paralog pairs. Per-family observed retention is the mean of
    `true_label` (0/1) across the same pairs. The correlation tells you
    how well the model's family-level summary tracks the empirical
    retention rate.

    Parameters
    ----------
    predictions_df : pandas.DataFrame
        Output of src/training/train.py, with columns 'family_id',
        'predicted_pi', 'true_label'.

    Returns
    -------
    dict
        {'pearson_r': float, 'p_value': float, 'n_families': int}
    """
    per_family = (
        predictions_df.groupby("family_id")
        .agg(
            predicted_pi=("predicted_pi", "mean"),
            observed_retention=("true_label", "mean"),
        )
        .reset_index()
    )
    if len(per_family) < 3:
        return {"pearson_r": float("nan"), "p_value": float("nan"), "n_families": len(per_family)}
    r, p = stats.pearsonr(
        per_family["predicted_pi"], per_family["observed_retention"]
    )
    return {
        "pearson_r": float(r),
        "p_value": float(p),
        "n_families": int(len(per_family)),
    }


def per_fold_auc(predictions_df):
    """Recompute AUC per fold from the predictions DataFrame.

    Useful for verifying the values reported by train.py and for
    ablation runs that read from disk rather than re-train.
    """
    rows = []
    for fold, group in predictions_df.groupby("fold"):
        rows.append(
            {
                "fold": int(fold),
                "n_pairs": int(len(group)),
                "auc_roc": compute_auc_roc(
                    group["true_label"], group["predicted_pi"]
                ),
            }
        )
    return pd.DataFrame(rows).sort_values("fold").reset_index(drop=True)
