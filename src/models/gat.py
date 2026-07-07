"""
Graph attention network for predicting paralog regulatory inheritance.

The model produces per-node embeddings via GAT message passing on the
regulatory network, then predicts a per-pair retention probability by
combining the embeddings of the two paralogs symmetrically. Per-family
regulatory inheritance probability (pi) is derived at inference time as
the mean of per-pair probabilities within each family.

Inputs are read from the consolidated artifact described in
docs/data.md section 3.5:
    x              - continuous per-node features, shape (N, F)
    family_id      - categorical family ID per node, shape (N,)
    edge_index     - regulatory edges, shape (2, E)
    paralog_edges  - paralog pairs as edges, shape (2, P)

Output: per-pair logits, shape (P,). Apply sigmoid to obtain retention
probabilities; use BCEWithLogitsLoss for training to keep numerics stable.
"""

import torch
from torch import nn
from torch_geometric.nn import GATConv


class RegNetGAT(nn.Module):
    """GAT that predicts paralog-pair regulatory retention probability.

    Architecture
    ------------
    1. Family ID is embedded and concatenated with continuous node features.
    2. An input projection lifts to `hidden_dim` so all GAT layers operate
       in the same space and residual connections are well-defined.
    3. A stack of `num_layers` GATConv layers (multi-head attention with
       output concatenation back to `hidden_dim`), each followed by
       LayerNorm, ELU, dropout, and a residual add.
    4. Per-pair prediction: the embeddings of the two paralog nodes are
       averaged (symmetric pooling, so swapping (g1, g2) gives the same
       output) and passed through a small MLP that returns a single logit
       per pair.

    The mean-pooling step is the place where pair symmetry is enforced.
    Concatenation would make output depend on input order, which is not
    a property paralog pairs have.
    """

    def __init__(
        self,
        num_continuous_features,
        num_families,
        family_embed_dim=16,
        hidden_dim=64,
        num_heads=4,
        num_layers=2,
        dropout=0.5,
    ):
        super().__init__()
        if hidden_dim % num_heads != 0:
            raise ValueError(
                f"hidden_dim ({hidden_dim}) must be divisible by num_heads "
                f"({num_heads}) so head outputs concatenate cleanly back to "
                f"hidden_dim."
            )
        head_dim = hidden_dim // num_heads

        # Categorical family ID -> learned embedding, then concatenate with
        # the continuous feature vector. The embedding lets the model treat
        # genes in the same family as related without forcing a particular
        # similarity structure.
        self.family_embedding = nn.Embedding(num_families, family_embed_dim)
        in_dim = num_continuous_features + family_embed_dim
        self.input_proj = nn.Linear(in_dim, hidden_dim)

        # Stack of GAT layers, all operating in hidden_dim space so
        # residuals are simple element-wise adds.
        self.convs = nn.ModuleList(
            [
                GATConv(hidden_dim, head_dim, heads=num_heads, dropout=dropout)
                for _ in range(num_layers)
            ]
        )
        self.norms = nn.ModuleList(
            [nn.LayerNorm(hidden_dim) for _ in range(num_layers)]
        )

        self.dropout = nn.Dropout(dropout)
        self.activation = nn.ELU()

        # Pair classifier: takes a symmetric pair embedding, returns one logit.
        self.pair_mlp = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim),
            nn.ELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, 1),
        )

    def encode_nodes(self, x, family_id, edge_index):
        """Compute per-node embeddings via GAT message passing.

        Returns
        -------
        torch.Tensor of shape (num_nodes, hidden_dim)
        """
        family_emb = self.family_embedding(family_id)
        h = torch.cat([x, family_emb], dim=-1)
        h = self.activation(self.input_proj(h))

        for conv, norm in zip(self.convs, self.norms):
            h_new = conv(h, edge_index)
            h_new = norm(h_new)
            h_new = self.activation(h_new)
            h_new = self.dropout(h_new)
            # Residual add. Both tensors have shape (N, hidden_dim) by
            # construction, so this is just element-wise addition.
            h = h + h_new
        return h

    def predict_pair_logits(self, node_embeddings, paralog_edges):
        """Per-pair logits from per-node embeddings.

        Symmetric pooling (mean of the two paralog embeddings) ensures the
        output for pair (g1, g2) equals the output for (g2, g1).
        """
        h1 = node_embeddings[paralog_edges[0]]
        h2 = node_embeddings[paralog_edges[1]]
        pair_repr = (h1 + h2) / 2.0
        return self.pair_mlp(pair_repr).squeeze(-1)

    def forward(self, x, family_id, edge_index, paralog_edges):
        """End-to-end: per-pair retention logits.

        Returns
        -------
        torch.Tensor of shape (num_pairs,), unactivated logits.
        Apply torch.sigmoid for probabilities or pass to
        BCEWithLogitsLoss for training.
        """
        h = self.encode_nodes(x, family_id, edge_index)
        return self.predict_pair_logits(h, paralog_edges)

    @torch.no_grad()
    def predict_family_pi(self, x, family_id, edge_index, paralog_edges):
        """Aggregate per-pair retention probability into per-family pi.

        For each family that contains at least one paralog pair, pi is the
        mean of the model's predicted retention probability across the
        pairs in that family. This is the per-family pi value that feeds
        the closed-form motif significance framework (Scruse et al. 2026).

        Returns
        -------
        dict[int, float]
            family_id -> predicted pi
        """
        logits = self.forward(x, family_id, edge_index, paralog_edges)
        probs = torch.sigmoid(logits)
        # Both members of any paralog pair share a family by construction
        # (families were derived as connected components of the paralog
        # graph), so taking family_id[paralog_edges[0]] is equivalent to
        # taking it for paralog_edges[1].
        pair_families = family_id[paralog_edges[0]]

        unique_families, inverse = torch.unique(
            pair_families, return_inverse=True
        )
        sums = torch.zeros(
            len(unique_families), dtype=probs.dtype, device=probs.device
        )
        counts = torch.zeros_like(sums)
        sums.scatter_add_(0, inverse, probs)
        counts.scatter_add_(0, inverse, torch.ones_like(probs))
        means = sums / counts
        return dict(zip(unique_families.cpu().tolist(), means.cpu().tolist()))
