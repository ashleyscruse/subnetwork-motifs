"""
Build the motif catalog: compute Z-scores, apply multiple testing correction,
and annotate significant motifs.

The catalog is the primary scientific output of this project.
"""

import pandas as pd
from statsmodels.stats.multitest import multipletests

from src.significance.formulas import z_score


def compute_catalog(motif_counts, n, m, pi_values):
    """
    For each enumerated motif, compute Z-score using exact formulas
    with GAT-estimated pi values.

    Args:
        motif_counts: list of dicts from enumerate.py
        n: total genes
        m: total families
        pi_values: dict of family_id -> predicted pi

    Returns:
        DataFrame with columns: families, topology, observed, expected, variance, z_score
    """
    raise NotImplementedError


def apply_multiple_testing(catalog_df, method="fdr_bh", alpha=0.05):
    """
    Apply multiple testing correction to the catalog Z-scores.

    Args:
        catalog_df: DataFrame with z_score column
        method: correction method (fdr_bh = Benjamini-Hochberg)
        alpha: significance level

    Returns:
        DataFrame with added columns: p_value, p_adjusted, significant
    """
    raise NotImplementedError


def annotate_orthologs(catalog_df, ortholog_mapping):
    """
    Map yeast gene families in significant motifs to human orthologs
    using Ensembl Compara.
    """
    raise NotImplementedError
