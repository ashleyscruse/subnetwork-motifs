"""
Baseline 3: Random forest on hand-crafted family features.

Features per family:
    - Family size (number of genes)
    - Mean expression correlation between paralogs
    - Mean regulatory target count
    - Mean sequence conservation
"""

from sklearn.ensemble import RandomForestClassifier


def extract_family_features(families, graph_data):
    """Extract hand-crafted features for each family."""
    raise NotImplementedError


def train_and_predict(features, labels, cv_splits):
    """Train RF with same CV splits as GAT. Returns predictions and AUC-ROC."""
    raise NotImplementedError
