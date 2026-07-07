# Contributing

## Getting Started

### 1. Clone the repo
```bash
git clone https://github.com/[USERNAME]/subnetwork-motifs.git
cd subnetwork-motifs
```

### 2. Set up your environment
```bash
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

### 3. Read the background
- Read the published paper: [Scruse, Arnold & Robinson 2026](https://link.springer.com/article/10.1007/s11538-025-01592-1)
- Read `docs/plan.md` for the project plan and timeline
- Read `data/README.md` for dataset descriptions

### 4. Understand the key concept

The core idea: gene regulatory networks have patterns (motifs) that might be biologically significant or might just be artifacts of gene duplication. Our published math can test which is which, but it needs one parameter per gene family: **pi_i** (the probability that regulatory relationships are inherited during duplication). We train a graph attention network to estimate pi_i from the network data.

## Project Structure

```
src/
├── data/          # Data loading and graph construction
├── models/        # GAT model
├── baselines/     # Four comparison methods
├── training/      # Training loop and evaluation
└── significance/  # Exact formulas, motif enumeration, catalog
```

## Workflow

1. Download and process data (see `data/README.md`)
2. Build the graph: `python -m src.data.build_graph`
3. Train the GAT: `python -m src.training.train --data data/processed/yeast_graph.pt`
4. Run baselines for comparison
5. Run significance testing: `python -m src.significance.enumerate`

On TACC, use the SLURM scripts in `hpc/`.

## Code Style

- Python 3.10+
- Use type hints where practical
- Keep functions short and single-purpose
- Every module has a docstring explaining what it does
- Notebooks are for exploration; production code goes in `src/`

## Git Workflow

- Work on a branch, not main
- Write descriptive commit messages
- Don't commit large data files (they're gitignored)
- Don't commit results (except summary CSVs and figures for the paper)
