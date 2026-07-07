"""
Baseline 1: Uniform pi (pi_i = 0.5 for all families).
Uninformed null model.
"""


def predict_pi(families):
    """Returns pi = 0.5 for every family."""
    return {fam: 0.5 for fam in families}
