"""
Baseline 4: Logistic regression on the same hand-crafted family features as the random forest.
Linear baseline to test whether the GAT's advantage comes from nonlinearity or graph structure.
"""

from sklearn.linear_model import LogisticRegression


def train_and_predict(features, labels, cv_splits):
    """Train LR with same CV splits as GAT. Returns predictions and AUC-ROC."""
    raise NotImplementedError
