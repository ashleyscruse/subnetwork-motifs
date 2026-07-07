---
type: project
status: active
start: 2026-04-01
target: 2026-12-31
tags: [research, gat, yeast, motifs]
---

# Preliminary Results Plan: Yeast Regulatory Motif Catalog

**Funding:** AUC DSI Startup Research Funds
**Timeline:** April 2026 - December 2026 (flexible, work fits around semester schedule)
**Target:** Preprint by Fall 2026, journal submission by December 2026
**Compute:** TACC (MSCF allocation)
**Repo:** subnetwork-motifs/

---

## Action items

- [x] Phase 0: Environment setup (venv, requirements, HPC scaffold)
- [x] Phase 1: Data acquisition (Harbison, Wapinski, Gasch, SGD proteome, paralog identity)
- [x] Phase 1: Data documentation ([docs/data.md](data.md), FAIR-aligned)
- [x] Phase 2: Build consolidated graph artifact (`yeast_graph.pt`)
- [x] Phase 3: GAT model + 5-fold CV training pipeline
- [ ] Phase 3: Hyperparameter search (Optuna)
- [ ] Phase 3: TACC deployment (Vista or Stampede 3)
- [ ] Phase 4: Baselines (uniform pi, sequence-similarity pi, RF, LR)
- [ ] Phase 5: Ablations (topology vs features, layer count, learning curves)
- [ ] Phase 6: Significance testing pipeline (Scruse 2026 closed-form formulas)
- [ ] Phase 6: Cross-species ortholog mapping (human via Ensembl Compara)
- [ ] Phase 7: Draft methods + results
- [ ] Phase 7: Arnold review
- [ ] Phase 7: Preprint submission (bioRxiv / arXiv)
- [ ] Phase 7: Journal submission
- [ ] Cross-cutting: Model card ([docs/model_card.md](model_card.md))

---

## What This Produces

A published paper presenting the first catalog of statistically significant regulatory motifs in yeast, using empirically estimated regulatory inheritance probabilities (pi values) from a graph attention network combined with the exact closed-form significance testing framework from Scruse et al. 2026.

This paper:
- Is a standalone research contribution
- Provides preliminary results for NSF IIBR proposal
- Validates the end-to-end pipeline (data, GAT, pi estimation, significance testing, catalog)
- Produces the first real numbers for the framework

---

## Target Journal

Bioinformatics (Oxford), BMC Bioinformatics, or PLOS Computational Biology.

---

## Paper Story (1 paragraph)

Network motifs have been used for two decades to identify regulatory design principles, but existing tools (FANMOD, mfinder) use randomization-based null models that ignore gene family identity and duplication history. We previously published exact closed-form expressions for the expected count and variance of subnetwork motifs under a gene duplication null model (Scruse et al. 2026), but these formulas require regulatory inheritance probabilities (pi_i) that cannot be measured directly. Here we train a graph attention network on the yeast regulatory network to estimate pi_i from network topology, expression, and sequence data, then apply the exact significance testing framework genome-wide to produce the first catalog of statistically significant regulatory motifs under an evolutionary null model. We identify [N] significant motifs out of [M] tested, map them to disease-relevant human orthologs, and show that the GAT outperforms random forest and logistic regression baselines at predicting regulatory inheritance.

---

## The Work

### Phase 0: Environment Setup

Set up the development and compute environment.

**Tasks:**
- Configure Python environment on TACC with PyTorch Geometric, networkx, pandas, scikit-learn
- Add requirements.txt to repo
- Create TACC SLURM job templates in `hpc/`
- Write data download/acquisition scripts
- Verify GPU access on TACC for GAT training

**Deliverables:**
- Working environment on TACC
- `requirements.txt`
- `hpc/` job templates
- Data download scripts

---

### Phase 1: Data Acquisition and Processing

Download, parse, and clean all datasets needed for the GAT.

**Datasets:**

| Dataset | Source | Format needed |
|---|---|---|
| Yeast TF-target regulatory network | Harbison et al. 2004 supplementary + YEASTRACT database | Edge list: TF gene -> target gene |
| Duplication fate labels | Wapinski et al. 2007 supplementary | Per duplicate pair: retained / lost / rewired |
| Gene family assignments | Wapinski et al. 2007 or Ensembl Compara | Gene -> family ID mapping |
| Gene expression | Gasch et al. 2000 (stress response) or Hughes et al. 2000 | Gene x condition matrix |
| Sequence similarity | Ensembl or BLAST between paralogs | Pairwise similarity scores within families |

**Tasks:**
- Download all datasets
- Write parsing scripts for each data source
- Resolve gene ID mappings across sources (systematic names vs common names)
- Produce clean, merged CSVs
- Compute and verify summary statistics

**Key numbers to verify:**
- ~203 TFs, ~6,000 target genes
- ~450 duplicate pairs with fate labels
- ~500 paralog groups

**Deliverables:**
- `data/harbison/` -- regulatory network edge list
- `data/wapinski/` -- duplication fate labels
- `data/families/` -- gene family assignments
- `data/expression/` -- expression matrix (optional for v1)
- `data/processed/` -- merged dataset ready for PyTorch Geometric
- Data summary statistics notebook

---

### Phase 2: Build the Graph

Convert processed data into a PyTorch Geometric Data object.

**Graph structure:**
- Nodes: all genes in the regulatory network
- Edges: directed regulatory relationships (TF -> target)
- Node features: expression profile, sequence similarity to paralogs, gene family ID (as learned embedding), in-degree, out-degree
- Labels: for genes with Wapinski duplication data, binary (kept regulation = 1, lost = 0)

**Tasks:**
- Write graph construction code
- Implement node feature extraction and normalization
- Verify edge directionality
- Check label distribution (class balance)
- Generate graph statistics (degree distribution, connected components)

**Deliverables:**
- `src/data/build_graph.py`
- `src/data/features.py`
- Notebook: graph statistics and label distribution

---

### Phase 3: Train the GAT

Train the graph attention network to predict regulatory inheritance.

**Model architecture:**
- GATConv layers (start with 2, ablate to 3-4)
- Attention heads (start with 4)
- Family-level pooling to aggregate gene embeddings per family
- MLP prediction head (sigmoid output for pi_i)
- Residual connections between layers
- Dropout for regularization

**Loss function:**
- Primary: weighted binary cross-entropy on Wapinski labels
- Secondary: self-supervised link prediction on held-out edges
- Regularization: consistency term on families with similar in-degree/out-degree and regulatory target overlap

**Training setup:**
- 5-fold cross-validation on Wapinski labels
- Bayesian hyperparameter optimization (learning rate, dropout, layer count, hidden dim)
- Early stopping on validation loss
- Track AUC-ROC, Pearson correlation between predicted and observed pi per family

**Deliverables:**
- `src/models/gat.py`
- `src/training/train.py`
- `src/training/evaluate.py`
- `hpc/train_gat.slurm`
- Results: AUC-ROC on each fold, mean and std

---

### Phase 4: Run Baselines

The GAT must beat four baselines to justify its use.

**Baselines:**
1. Uniform pi (pi_i = 0.5 for all families)
2. Sequence-similarity pi (pi_i proportional to mean paralog sequence identity)
3. Random forest on hand-crafted family features (family size, mean expression correlation, regulatory target count, sequence conservation)
4. Logistic regression on the same features

**Tasks:**
- Implement all four baselines
- Run on same CV splits as the GAT
- Compare AUC-ROC across all methods

**Note:** If the GAT advantage over random forest is small (< 5 AUC points), the paper shifts framing to "comparison of estimation methods" rather than "GAT is best." The significance testing and catalog are still valid regardless of which method estimates pi. The math doesn't care.

**Deliverables:**
- `src/baselines/` -- all four implementations
- Comparison table: AUC-ROC for GAT vs all baselines (goes directly into paper)

---

### Phase 5: Ablation Studies

Understand which components of the GAT matter.

**Ablations:**
- Topology-only (no node features) vs full model
- Features-only (no message passing) vs full model
- Layer count: 2, 3, 4, 5 layers
- Learning curves: 10%, 25%, 50%, 75%, 100% of labels

**Deliverables:**
- Ablation results table (goes directly into paper)
- Learning curve plot (goes directly into paper)

---

### Phase 6: Significance Testing and Catalog

This is where the published math meets the GAT. The core contribution of the paper.

**Steps:**
1. Take the trained GAT's predicted pi_i for every yeast gene family
2. For all subnetwork motifs at k=2, 3, 4 (k=5 with GNN pruning if feasible):
   - Compute E(M) and Var(M) using exact formulas from Scruse et al. 2026
   - Count observed motif instances in the yeast regulatory network
   - Compute Z-score: Z = (X_obs - E(M)) / sqrt(Var(M))
3. Apply multiple testing correction (Bonferroni or FDR)
4. Catalog all significant motifs (|Z| > threshold after correction)
5. Map significant motifs to human orthologs via Ensembl Compara

**Deliverables:**
- `src/significance/formulas.py` -- E(M), Var(M) from published theory
- `src/significance/enumerate.py` -- motif counting
- `src/significance/catalog.py` -- Z-score computation and catalog generation
- The catalog: table of significant motifs with gene families, Z-scores, pi values, human orthologs

---

### Phase 7: Write the Paper

**Sections:**

| Section | Content |
|---|---|
| Introduction | Problem (existing tools ignore duplication), gap (pi unmeasurable), contribution (GAT + exact formulas = first catalog) |
| Background | Subnetwork motifs, FANMOD/mfinder limitations, pi estimation problem |
| Methods: Data | Datasets, processing, graph construction |
| Methods: GAT Architecture | Model design, loss function, training procedure |
| Methods: Significance Testing | Exact formulas from Scruse et al. 2026, enumeration, multiple testing |
| Results: GAT vs Baselines | AUC-ROC comparison table |
| Results: Ablations | What matters, what doesn't |
| Results: Motif Catalog | The catalog itself, biological interpretation, MAPK example |
| Discussion | What this enables, limitations, future (cross-species, NSF proposal) |

**Writing approach:**
- Draft methods and results directly from experiment notebooks and logs
- Introduction and discussion written last (need results to frame correctly)
- Arnold reviews before submission
- Preprint on bioRxiv or arXiv, then submit to journal

---

## Realistic Timeline

Work fits around the semester. Heavy coding sprints happen during breaks and summer. Writing overlaps with later phases.

| When | What | Key output |
|---|---|---|
| April - May 2026 (semester) | Phase 0 + Phase 1 | Environment ready, data downloaded and parsed |
| June - July 2026 (summer) | Phase 2 + Phase 3 + Phase 4 | Graph built, GAT trained, baselines run, AUC-ROC numbers |
| August 2026 (late summer) | Phase 5 + Phase 6 | Ablations done, significance testing run, catalog produced |
| September - October 2026 | Phase 7 | Paper drafted, Arnold reviews |
| November - December 2026 | Revisions + submission | Preprint posted, journal submission |

If student is available, phases accelerate. If not, the work is still doable solo with AI assistance. The coding is not the bottleneck; the data parsing and biological verification are.

---

## What AI Tools Accelerate

AI (Claude, Copilot, etc.) is used as a coding and writing assistant. It does not design experiments, interpret results, or make architectural decisions.

It generates: boilerplate code, data parsing scripts, plotting code, LaTeX formatting, SLURM scripts, and speeds up debugging.

Conservative estimate: AI reduces coding time by 40-50%.

---

## What This Enables

- **NSF IIBR submission:** Proposal includes real AUC-ROC numbers and the first motif catalog
- **The science report:** The MAPK Z-score example becomes a real number
- **The full platform:** Yeast is done. Cross-species transfer is the next grant's funded work.

---

## Risks and Mitigations

| Risk | Mitigation |
|---|---|
| GAT doesn't outperform random forest | Random forest pi estimates still produce a valid catalog. Paper becomes "comparison of estimation methods." Still publishable. |
| Data harder to parse than expected | Use YEASTRACT (web-downloadable) as backup for regulatory network. Budget extra time in Phase 1. |
| Motif enumeration too slow at k=4+ | Start with k=2,3 (guaranteed fast). k=4 on TACC should be feasible. k=5+ is a stretch goal for the paper. |
| Writing takes longer than planned | Preprint can go out with k=2,3 results. k=4+ added in revision. |
