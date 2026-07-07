"""
Baseline 2: Sequence-similarity pi.
pi_i is proportional to mean paralog sequence identity within each family.
"""


def predict_pi(families, similarity_scores):
    """
    For each family, pi = mean pairwise sequence identity among paralogs.
    Families with higher sequence conservation get higher pi.
    """
    raise NotImplementedError
