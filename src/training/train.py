"""
5-fold cross-validation training loop for the paralog retention GAT.

Loads the consolidated graph artifact (data/processed/yeast_graph.pt),
trains five models using the precomputed `family_fold` field for splits,
and writes per-pair retention predictions to results/predicted_pi.csv.

Loss
----
Per-pair binary cross-entropy with positive-class reweighting. The
positive class (regulation retained) outweighs the negative class roughly
88/12 in the Wapinski data; pos_weight = n_neg / n_pos rebalances the
gradient so the rare class is not ignored.

Splits
------
Cross-validation folds are computed at the family level inside the data
artifact (see docs/data.md section 3.4). Both members of any paralog pair
share a family, so taking the fold of either member identifies the pair's
fold unambiguously. Pairs with fold == -1 (genes in families with no
labeled pairs) cannot occur in this artifact but are checked for safety.

Outputs
-------
results/predicted_pi.csv with columns:
    gene1_idx, gene2_idx, gene1, gene2, family_id,
    true_label, predicted_pi, fold

The per-family pi values that feed Phase 6 significance testing are the
mean of `predicted_pi` within each `family_id` group.

Usage
-----
    python -m src.training.train
    python -m src.training.train --device cuda

Defaults are deliberately fixed; pass through `HParams` to override or
wire up Optuna in a follow-up.
"""

import argparse
from dataclasses import dataclass, asdict
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from sklearn.metrics import roc_auc_score
from torch import nn

from src.models.gat import RegNetGAT

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_GRAPH = REPO_ROOT / "data" / "processed" / "yeast_graph.pt"
DEFAULT_OUTPUT = REPO_ROOT / "results" / "predicted_pi.csv"

NUM_FOLDS = 5
RANDOM_SEED = 0


@dataclass
class HParams:
    """Fixed defaults for the first end-to-end run.

    All values are research starting points, not tuned. Hyperparameter
    search (Optuna) is a follow-up that will sample over these.
    """

    family_embed_dim: int = 16
    hidden_dim: int = 64
    num_heads: int = 4
    num_layers: int = 2
    dropout: float = 0.5
    learning_rate: float = 1e-3
    weight_decay: float = 1e-5
    max_epochs: int = 200
    patience: int = 20


def set_seed(seed):
    """Seed Python and torch RNGs so each fold is reproducible.

    Cross-fold variation comes from the fold split itself plus the
    per-fold seed offset; within a fold, model init and dropout are
    deterministic given the seed.
    """
    torch.manual_seed(seed)
    np.random.seed(seed)


def auc_safe(probs, labels):
    """ROC AUC with a guard for single-class validation folds."""
    if len(np.unique(labels)) < 2:
        return float("nan")
    return float(roc_auc_score(labels, probs))


def train_one_fold(data, fold, hp, device, verbose=True):
    """Train and evaluate a single CV fold.

    Returns
    -------
    best_val_auc : float
    best_epoch : int
    best_logits : torch.Tensor of shape (num_pairs,)
        Logits at the epoch that achieved best validation AUC. Used to
        record predictions for the fold's validation pairs.
    """
    set_seed(RANDOM_SEED + fold)

    # Pair-level fold derived from per-node family_fold. Both members of
    # a pair share a family (and therefore a fold) by construction.
    pair_fold = data.family_fold[data.paralog_edges[0]]
    val_mask = pair_fold == fold
    train_mask = (pair_fold != fold) & (pair_fold != -1)

    # Class-imbalance reweighting: pos_weight is applied to the positive
    # term in BCEWithLogitsLoss so the rarer class contributes more.
    train_labels = data.paralog_labels[train_mask]
    n_pos = (train_labels == 1).sum().item()
    n_neg = (train_labels == 0).sum().item()
    pos_weight = torch.tensor([n_neg / max(n_pos, 1)], device=device)

    num_families = int(data.family_id.max().item()) + 1
    model = RegNetGAT(
        num_continuous_features=data.x.shape[1],
        num_families=num_families,
        family_embed_dim=hp.family_embed_dim,
        hidden_dim=hp.hidden_dim,
        num_heads=hp.num_heads,
        num_layers=hp.num_layers,
        dropout=hp.dropout,
    ).to(device)

    optimizer = torch.optim.Adam(
        model.parameters(),
        lr=hp.learning_rate,
        weight_decay=hp.weight_decay,
    )
    criterion = nn.BCEWithLogitsLoss(pos_weight=pos_weight)

    train_targets = data.paralog_labels[train_mask].float()
    val_labels_np = data.paralog_labels[val_mask].cpu().numpy()

    best_val_auc = -float("inf")
    best_epoch = 0
    best_logits = None
    epochs_since_best = 0

    for epoch in range(1, hp.max_epochs + 1):
        # Train step (full-graph batch; the network is small enough that
        # mini-batching does not help).
        model.train()
        optimizer.zero_grad()
        logits = model(
            data.x, data.family_id, data.edge_index, data.paralog_edges
        )
        loss = criterion(logits[train_mask], train_targets)
        loss.backward()
        optimizer.step()

        # Eval step on the held-out fold.
        model.eval()
        with torch.no_grad():
            logits = model(
                data.x, data.family_id, data.edge_index, data.paralog_edges
            )
            val_probs = torch.sigmoid(logits[val_mask]).cpu().numpy()
            val_auc = auc_safe(val_probs, val_labels_np)

        if val_auc > best_val_auc:
            best_val_auc = val_auc
            best_epoch = epoch
            best_logits = logits.detach().clone()
            epochs_since_best = 0
        else:
            epochs_since_best += 1

        if verbose and epoch % 20 == 0:
            print(
                f"    epoch {epoch:3d}  train_loss={loss.item():.4f}  "
                f"val_auc={val_auc:.4f}  (best {best_val_auc:.4f} @ ep {best_epoch})"
            )

        if epochs_since_best >= hp.patience:
            if verbose:
                print(
                    f"    early stop at epoch {epoch} "
                    f"({hp.patience} epochs without val improvement)"
                )
            break

    return best_val_auc, best_epoch, best_logits


def cross_validate(data, hp, device):
    """Run NUM_FOLDS independent training runs; record per-fold metrics."""
    pair_fold_np = data.family_fold[data.paralog_edges[0]].cpu().numpy()
    pair_predictions = np.full(data.paralog_edges.shape[1], np.nan)
    fold_aucs = []
    fold_epochs = []

    for fold in range(NUM_FOLDS):
        n_train = ((pair_fold_np != fold) & (pair_fold_np != -1)).sum()
        n_val = (pair_fold_np == fold).sum()
        print(f"\n=== Fold {fold} ===  train pairs: {n_train}, val pairs: {n_val}")
        best_auc, best_epoch, best_logits = train_one_fold(
            data, fold, hp, device
        )
        fold_aucs.append(best_auc)
        fold_epochs.append(best_epoch)

        # Record this fold's val-set predictions. Each pair appears in the
        # validation set of exactly one fold, so pair_predictions ends up
        # fully populated.
        val_mask_np = pair_fold_np == fold
        val_probs = torch.sigmoid(best_logits.cpu()).numpy()
        pair_predictions[val_mask_np] = val_probs[val_mask_np]
        print(
            f"  fold {fold} done: best val AUC = {best_auc:.4f} (epoch {best_epoch})"
        )

    return fold_aucs, fold_epochs, pair_predictions, pair_fold_np


def write_predictions(data, pair_predictions, pair_fold_np, output_path):
    p_edges = data.paralog_edges.cpu().numpy()
    df = pd.DataFrame(
        {
            "gene1_idx": p_edges[0],
            "gene2_idx": p_edges[1],
            "gene1": [data.gene_names[i] for i in p_edges[0]],
            "gene2": [data.gene_names[i] for i in p_edges[1]],
            "family_id": data.family_id[p_edges[0]].cpu().numpy(),
            "true_label": data.paralog_labels.cpu().numpy(),
            "predicted_pi": pair_predictions,
            "fold": pair_fold_np,
        }
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, index=False)
    return df


def main():
    parser = argparse.ArgumentParser(
        description="5-fold CV training for the paralog retention GAT"
    )
    parser.add_argument("--graph", type=Path, default=DEFAULT_GRAPH)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--device",
        default="cuda" if torch.cuda.is_available() else "cpu",
        help="torch device (cpu or cuda).",
    )
    args = parser.parse_args()

    print(f"Loading {args.graph}")
    data = torch.load(args.graph, weights_only=False)
    print(
        f"  N={data.num_nodes:,}  E={data.edge_index.shape[1]:,}  "
        f"P={data.paralog_edges.shape[1]:,}"
    )

    hp = HParams()
    print(f"Device: {args.device}")
    print(f"Hyperparameters: {asdict(hp)}")

    data = data.to(args.device)
    fold_aucs, fold_epochs, pair_predictions, pair_fold_np = cross_validate(
        data, hp, args.device
    )

    print("\n========== Cross-validation summary ==========")
    for i, (auc, ep) in enumerate(zip(fold_aucs, fold_epochs)):
        print(f"  fold {i}: AUC {auc:.4f} (best ep {ep})")
    valid_aucs = [a for a in fold_aucs if not np.isnan(a)]
    if valid_aucs:
        print(
            f"Mean ± std: {np.mean(valid_aucs):.4f} ± "
            f"{np.std(valid_aucs):.4f}  (n={len(valid_aucs)})"
        )

    df = write_predictions(data, pair_predictions, pair_fold_np, args.output)
    print(f"\nWrote {args.output}")

    # Per-family pi summary (mean predicted retention prob within each family).
    family_pi = df.groupby("family_id")["predicted_pi"].mean()
    print(
        f"Per-family pi range: [{family_pi.min():.3f}, {family_pi.max():.3f}], "
        f"median {family_pi.median():.3f}, n_families={len(family_pi):,}"
    )


if __name__ == "__main__":
    main()
