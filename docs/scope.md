---
title: Your Scope
layout: default
---

<div class="page-nav">
  <a href="./">Getting Started</a>
  <a href="scope.html" class="current">Your Scope</a>
</div>

# Your Scope

## What you own

Four baseline estimators of the regulatory-inheritance probability π_i (Phase 4 in `docs/plan.md`). These plug into the seven-way cross-validation the paper promises. The four files to write:

- `src/baselines/uniform.py`: π_i = 0.5 for all families. A trivial floor; also computes downstream Z-scores as a sanity baseline.
- `src/baselines/sequence_similarity.py`: π_i proportional to mean paralog percent identity within each family, using `data/families/sequence_similarity.csv`.
- `src/baselines/random_forest.py`: random forest classifier on hand-crafted family features (family size, mean expression correlation, regulatory target count, sequence conservation).
- `src/baselines/logistic_regression.py`: same features as the random forest, logistic regression classifier.

## Rules

- Each baseline is a function that takes the loaded `data/processed/yeast_graph.pt` and returns a dict `{family_id: pi_value}`.
- Use the same 5-fold family-level cross-validation as the GAT (`data.family_fold`).
- Report AUC-ROC on the paralog-pair task, matching the output format of `src/training/train.py`.
- Save per-baseline predictions to `results/baseline_pi_<method>.csv` using the same column set as `results/predicted_pi.csv`.
- Do not modify `src/models/gat.py` or `src/training/train.py`. If you need shared utilities, add a new file or ask Ashley first.

## Where to start

Start with `src/baselines/uniform.py`. It is the shortest and lets you get end-to-end wiring right before you tackle the more involved baselines.

Read first, in this order:

1. `docs/plan.md`: phased project plan. You are Phase 4.
2. `docs/data.md`: data documentation. Explains what each field on `data/processed/yeast_graph.pt` means.
3. `src/training/train.py`: mirror its cross-validation loop structure and output format.
4. `src/data/features.py`: pure feature functions you will reuse for the random-forest and logistic-regression baselines.

## Success measure

By the end of the summer, produce a single results notebook (in `notebooks/`) that compares your four baselines and the GAT on:

1. AUC-ROC on the Wapinski paralog-pair task (5-fold CV).
2. Mean squared error against synthetic true-π profiles.
3. Downstream Z-scores of the top motifs.

## Non-overlap with the other collaborator

Another collaborator on the paper owns two additional estimators (an evidence-code based method and a count-based method) built around a seven-family Gene Ontology framework and their own web application. Your work uses the 599-family Wapinski framework in this repo. Same math, different granularity, non-overlapping scope. Do not build the GO-family framework or the evidence-code estimator.

## Workflow reminders

- Work on the `delva` branch, never on `main`.
- Commit often with short, clear messages.
- Open a pull request from `delva` into `main` when you want work reviewed and merged.
- If you use an AI assistant, point it at `docs/plan.md` and `docs/data.md` first, and keep it on your branch.
