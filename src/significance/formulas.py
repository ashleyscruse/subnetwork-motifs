"""
Exact closed-form expressions for subnetwork motif expected count and variance
under the gene duplication model.

From: Scruse, Arnold & Robinson (2026), "Counting Subnetworks Under Gene
Duplication in Genetic Regulatory Networks," Bulletin of Mathematical Biology.

Full Duplication Model (pi_i = 1 for all families):

    E(M) = C(n+k-1, m+k-1) / C(n-1, m-1)

    E(M^2) = sum_{i=0}^{k} 2^i * C(k,i) * C(n+k-1, m-1+k+i) / C(n-1, m-1)

    Var(M) = E(M^2) - E(M)^2

Where:
    n = total number of genes
    m = number of gene families
    k = motif size (number of families in the subnetwork motif)
    C(a,b) = binomial coefficient "a choose b"

The Partial Duplication Model (0 <= pi_i <= 1) generalizes these expressions.
See the published paper for the full derivation.
"""

from math import comb


def expected_count(n, m, k):
    """
    Expected number of subnetwork motif instances under the Full Duplication model.

    Args:
        n: total number of genes
        m: number of gene families
        k: motif size
    """
    return comb(n + k - 1, m + k - 1) / comb(n - 1, m - 1)


def second_moment(n, m, k):
    """
    Second moment E(M^2) under the Full Duplication model.
    """
    total = 0
    for i in range(k + 1):
        total += (2 ** i) * comb(k, i) * comb(n + k - 1, m - 1 + k + i)
    return total / comb(n - 1, m - 1)


def variance(n, m, k):
    """Variance of subnetwork motif count under the Full Duplication model."""
    em = expected_count(n, m, k)
    em2 = second_moment(n, m, k)
    return em2 - em ** 2


def z_score(observed, n, m, k):
    """
    Z-score for a subnetwork motif.

    Args:
        observed: observed count of this motif in the network
        n, m, k: parameters for the duplication model

    Returns:
        Z-score. |Z| > 1.96 suggests significance at the 5% level.
    """
    em = expected_count(n, m, k)
    var = variance(n, m, k)
    if var <= 0:
        return 0.0
    return (observed - em) / (var ** 0.5)


def partial_duplication_expected_count(n, m, k, pi_values):
    """
    Expected count under the Partial Duplication model.
    pi_values is a list of pi_i for the k families in the motif.

    TODO: Implement from the published paper's partial duplication expressions.
    """
    raise NotImplementedError


def partial_duplication_variance(n, m, k, pi_values):
    """
    Variance under the Partial Duplication model.

    TODO: Implement from the published paper's partial duplication expressions.
    """
    raise NotImplementedError
