# Data Documentation: Yeast Regulatory Genomics Dataset

**Dataset:** *Saccharomyces cerevisiae* regulatory genomics dataset (regulatory edges, paralog pairs with duplication-fate labels, expression matrix, reference proteome, derived gene families and paralog sequence identity)

**Maintainer:** Ashley Scruse, Ph.D. (Morehouse College / AUC DSI)

**Last reviewed:** 2026-05-01

This document describes a curated *S. cerevisiae* regulatory genomics dataset assembled from five public, peer-reviewed sources. It covers what each input is, where and how it was obtained, what processing was applied, how the inputs combine into a single artifact, and how to reproduce every step. It is the authoritative reproducibility reference for the dataset.

---

## 1. Scope

This document describes a curated *Saccharomyces cerevisiae* regulatory-genomics dataset assembled from five public, peer-reviewed sources. The combined dataset comprises a transcription-factor-to-target regulatory network, paralog pair classifications with duplication-fate labels, an environmental stress expression matrix, a reference proteome, and two derived artifacts (gene family assignments and paralog sequence-identity scores).

Every input is keyed by systematic ORF names from the Saccharomyces Genome Database, so the sources combine via simple equijoins with no ID mapping. The processing pipeline lives in [`src/data/`](../src/data/) and writes deterministic outputs from deterministic inputs.

How any specific downstream model uses these data is out of scope for this document and belongs in a model card.

---

## 2. Data sources

Each subsection below covers one dataset. The structure is identical across sources: **what it is, where we got it, what the file looks like, what we extracted, what to watch out for**.

### 2.1 Harbison 2004: regulatory edges

| | |
|---|---|
| **What it is** | Genome-wide ChIP-chip binding profiles for ~203 *S. cerevisiae* transcription factors against ~6,300 ORFs, profiled under 1+ growth conditions per TF. The canonical yeast regulatory network. |
| **Citation** | Harbison CT, Gordon DB, Lee TI, *et al.* (2004). Transcriptional regulatory code of a eukaryotic genome. *Nature* 431:99-104. DOI [10.1038/nature02800](https://doi.org/10.1038/nature02800) |
| **Original hosting** | Young Lab, MIT (`younglab.wi.mit.edu/regulatory_code/`, mirrored at `jura.wi.mit.edu/young_public/regulatory_code/`). **Both URLs are dead.** |
| **Current source** | Internet Archive Wayback snapshot from 2015-09-16: `http://web.archive.org/web/20150916025747if_/http://jura.wi.mit.edu/young_public/regulatory_code/files_for_paper_txt.zip` |
| **Raw files** | `pvalbygene_forpaper_ypd_txt_withheader.txt` (TFs profiled in YPD) and `pvalbygene_forpaper_otherconditions_txt_withheader.txt` (TFs profiled under stress). Together they cover all assayed TF -> ORF p-values. |
| **Raw format** | Tab-separated. Two header rows (experiment IDs, then TF_condition labels like `ABF1_YPD` or `ARG81_SM`). Three metadata columns (ORF, common name, description) then one column of ChIP-chip binding p-values per (TF, condition). Missing values: `NaN`. |
| **Output** | [`data/harbison/edges.csv`](../data/harbison/edges.csv): one row per (TF, target) pair where any condition crossed `p <= 0.001`. Columns: `tf`, `target`, `min_pvalue`, `n_conditions_bound`. |
| **Threshold** | `p <= 0.001`, the cutoff used in Harbison 2004 itself and standard in the downstream literature (Milo, MacIsaac). Override via `--threshold` if a looser cutoff is needed for sensitivity analysis. |
| **Aggregation** | TFs profiled under multiple conditions are unioned: a (TF, target) pair appears if **any** condition shows `p <= 0.001`. `min_pvalue` is the best p-value across conditions; `n_conditions_bound` counts conditions where the threshold held. |
| **Caveats** | (i) Multi-condition aggregation slightly inflates false positives versus per-condition analysis. (ii) TF labels retain spaces and parentheses where Harbison used them (e.g. `A1 (MATA1)`); preserve exactly when joining. (iii) MacIsaac 2006 (DOI [10.1186/1471-2105-7-113](https://doi.org/10.1186/1471-2105-7-113)) provides a conservation-filtered refinement; we use original Harbison for direct comparability with the network-motif literature. |
| **Output stats** | 10,853 edges, 177 TFs with at least one significant target, 3,839 unique targets. |

### 2.2 Wapinski 2007: paralog pairs with duplication-fate labels

| | |
|---|---|
| **What it is** | Per-pair classifications of *S. cerevisiae* paralogs derived from the SYNERGY phylogenomic reconstruction across 17 fungal genomes. For each pair: ancestral duplication point, GO category divergence, transcription-module divergence, and protein-interaction-network conservation/divergence statistics. |
| **Citation** | Wapinski I, Pfeffer A, Friedman N, Regev A (2007). Natural history and evolutionary principles of gene duplication in fungi. *Nature* 449:54-61. DOI [10.1038/nature06107](https://doi.org/10.1038/nature06107) |
| **Original hosting** | Broad Institute (`broadinstitute.org/regev/orthogroups/`). **Hosting is offline; the bulk text dump is not preserved on the Wayback Machine** (only fragments of per-orthogroup HTML pages remain). |
| **Current source** | Nature publisher static-content host: `https://static-content.springer.com/esm/art%3A10.1038%2Fnature06107/MediaObjects/`. We pull MOESM276 (supplementary methods PDF), MOESM278 (transcription modules), and MOESM281 (the paralog pair table that drives this work). `--all` retrieves all seven supplements. |
| **Raw format** | MOESM281.xls / sheet "PARALOGS": two header rows (group titles in row 0, sub-column names in row 1), then 801 paralog-pair rows. The sub-columns relevant for us are `GENE1`, `GENE2`, `DUPLICATION POINT`, four GO/transcription-module migration flags, and six protein-interaction-network statistics. |
| **Output** | [`data/wapinski/paralogs.csv`](../data/wapinski/paralogs.csv): 801 rows, 14 columns, including a derived `regulation_retained` flag (defined below). |
| **Label derivation** | `regulation_retained = 1 - transcription_module`. The raw `transcription_module` column is `1` when paralogs migrated apart between transcription modules (regulation diverged) and `0` when they share a module (regulation kept). The flip makes the label intuitive: `1 = regulation kept`, `0 = regulation lost`. |
| **Caveats** | (i) Class balance is heavily skewed: 702 retained / 99 lost (88% / 12%). (ii) Wapinski's transcription-module assignment is itself a derived clustering (MOESM278); col 6 is ground truth only insofar as that clustering is reliable. (iii) The sheet provides paralog **pairs**, not multi-gene families; family IDs are derived in section 2.5 below. (iv) 801 total pairs include WGD and older small-scale duplications; filter on `duplication_point == "WGD"` for a WGD-only subset (438 pairs). |
| **Output stats** | 801 pairs (438 WGD), 1,297 distinct genes, label balance 702/99. |

### 2.3 Gasch 2000: environmental stress expression matrix

| | |
|---|---|
| **What it is** | Genome-wide microarray expression for *S. cerevisiae* under 173 stress and growth-transition conditions: heat shock, hydrogen peroxide, menadione, diamide, dithiothreitol, hyper- and hypo-osmotic shock, amino-acid starvation, nitrogen depletion, stationary phase. The canonical environmental stress response (ESR) reference dataset. |
| **Citation** | Gasch AP, Spellman PT, Kao CM, *et al.* (2000). Genomic expression programs in the response of yeast cells to environmental changes. *Mol. Biol. Cell* 11:4241-4257. DOI [10.1091/mbc.11.12.4241](https://doi.org/10.1091/mbc.11.12.4241) |
| **Original hosting** | `genome-www.stanford.edu/yeast_stress/`. **Offline.** |
| **Current source** | Internet Archive Wayback snapshot from 2021-10-19: `http://web.archive.org/web/20211019175531if_/http://genome-www.stanford.edu/yeast_stress/data/rawdata/complete_dataset.txt`. The same measurements are also in NCBI GEO; we use the original Stanford companion file because it was the canonical published export. |
| **Raw format** | Tab-separated, single header row. Columns: `UID` (ORF systematic name), `NAME` (free-text annotation embedding ORF + gene name + ontology + SGD ID), `GWEIGHT` (always 1, legacy from clustering software), then 173 condition-label columns of log2 expression ratios (treated / control). Missing cells are blank. |
| **Output** | [`data/expression/expression.csv`](../data/expression/expression.csv) (wide: 6,152 genes x 173 conditions) and [`data/expression/expression_long.csv`](../data/expression/expression_long.csv) (tidy: gene, condition, log2_ratio, NaN dropped). |
| **Caveats** | (i) ~3% of cells are missing (1,032,274 of 1,064,296); imputation strategy is left to the data consumer. (ii) Condition labels are verbatim Gasch strings (e.g. `Heat Shock 05 minutes hs-1`); they are stable identifiers within this file but are not standardized to any community vocabulary. (iii) Values are log2 ratios, not absolute expression; sign is meaningful. |
| **Output stats** | 6,152 ORFs, 173 conditions, 97% cell coverage. |

### 2.4 SGD reference proteome

| | |
|---|---|
| **What it is** | *S. cerevisiae* S288C reference protein FASTA, used to compute paralog sequence identity. |
| **Source** | Saccharomyces Genome Database, file `sequence/S288C_reference/orf_protein/orf_trans.fasta.gz` from `sgd-archive.yeastgenome.org`. |
| **License** | CC-BY 4.0 (see [SGD use policy](https://www.yeastgenome.org/help)). |
| **Citation** | Cherry JM, Hong EL, Amundsen C, *et al.* (2012). Saccharomyces Genome Database: the genomics resource of budding yeast. *Nucleic Acids Research* 40(D1):D700-D705. DOI [10.1093/nar/gkr1029](https://doi.org/10.1093/nar/gkr1029) |
| **Format** | FASTA. Headers begin with the systematic ORF name as the first whitespace token (e.g. `>YAL001C TFC3 SGDID:S0000001 ...`), which allows direct joining against Wapinski paralogs and Harbison edges without ID mapping. |
| **Output** | [`data/families/raw/orf_trans.fasta`](../data/families/raw/) (~5.2 MB uncompressed, ~6,039 sequences). |

### 2.5 Derived: gene families and paralog sequence identity

Two artifacts produced from the data above without an external download:

#### Gene families ([`data/families/families.csv`](../data/families/families.csv))

Family IDs are connected components of the *Saccharomyces cerevisiae* paralog graph in Wapinski 2007 supplementary table S5 (MOESM281). Two genes are in the same family if they appear together in any paralog pair, transitively. Singleton ORFs (no paralog) do not appear in this output.

Columns: `gene` (systematic name), `family_id` (synthetic, format `FAM####`), `family_size`.

Output: 599 multi-gene families covering 1,297 genes.

#### Paralog sequence identity ([`data/families/sequence_similarity.csv`](../data/families/sequence_similarity.csv))

For each Wapinski paralog pair, both protein sequences are looked up in `orf_trans.fasta` and aligned globally with Biopython's `PairwiseAligner` (BLOSUM62, gap-open `-10`, gap-extend `-0.5`). Percent identity is `100 * n_matches / aln_length`, where matches exclude gap positions.

Columns: `gene1`, `gene2`, `length1`, `length2`, `aln_length`, `n_matches`, `percent_identity`.

Output stats: 791 of 801 paralog pairs aligned. Mean percent identity 45.2%, median 40.4% -- consistent with the ~100 Mya age of the yeast WGD. Ten pairs were skipped because at least one ORF (e.g. `YKR004C-A`, `YIL168W`, `YER186W-A`, `YAR062W`) is not in SGD's standard reference proteome; these tend to be dubious or short ORFs and are absent from this output.

---

## 3. Combined dataset

The five sources above are joined into a single graph artifact at [`data/processed/yeast_graph.pt`](../data/processed/), built deterministically by [`src/data/build_graph.py`](../src/data/build_graph.py) and serialized as a PyTorch Geometric `Data` object.

### 3.1 Composition

**Nodes** are yeast genes. The node set is the union of:

- TFs and targets appearing in `harbison/edges.csv`
- `gene1` and `gene2` from `wapinski/paralogs.csv`
- ORFs in `expression/expression.csv` with at least one non-missing value

Each node is assigned an integer index in `0..N-1`. Genes appearing in only a subset of sources still become nodes; per-source attributes that are not available are filled per the rules in section 3.2. The `gene_names` field preserves the systematic ORF name at each index.

**Edges** are directed regulatory relationships from `harbison/edges.csv`: TF -> target. The `edge_index` field has shape `[2, num_edges]` with source node indices in row 0 and destination node indices in row 1. Per-edge attributes `min_pvalue` and `n_conditions_bound` from the Harbison source are preserved on each edge.

### 3.2 Per-node attributes

| Attribute | Source | Type | Missing-value rule |
|---|---|---|---|
| `family_id` | `families/families.csv` | int | Genes not in any Wapinski paralog pair receive a unique singleton family ID |
| `in_degree` | `harbison/edges.csv` | int | 0 if no incoming edges |
| `out_degree` | `harbison/edges.csv` | int | 0 if no outgoing edges |
| `expression` | `expression/expression.csv` | float[173] | Missing condition values are imputed to 0.0 (no-change baseline); genes absent from the matrix get an all-zero vector |
| `mean_paralog_identity` | `families/sequence_similarity.csv` | float | 0.0 for singletons or pairs missing from the FASTA-derived alignment |

### 3.3 Labels

Labels are stored at the paralog pair level rather than projected onto individual genes. A single gene can belong to multiple paralog pairs (multi-paralog families), and projection would force aggregation that discards evolutionary detail the Scruse 2026 framework operates on. The artifact carries pair-level labels directly:

- `paralog_edges`: integer tensor of shape `[2, num_pairs]`, in the same PyG edge format as `edge_index`. Each column is a paralog pair from `wapinski/paralogs.csv`. Both rows index into the same node set as the regulatory edges, so paralog pairs and regulatory edges live in the same node space.
- `paralog_labels`: integer tensor of shape `[num_pairs]` carrying `regulation_retained` per pair:
  - `1` = paralogs share their transcription module (regulation retained)
  - `0` = paralogs migrated to different transcription modules (regulation lost or rewired)

Class balance across pairs follows the underlying Wapinski distribution (702 retained / 99 lost, ~88% / 12%).

Per-pair storage preserves every paralog relationship individually. Consumers that need a per-node aggregate can derive one from `paralog_edges` and `paralog_labels` (e.g. majority vote, any-divergence-wins, or a soft mean), with the choice of aggregation rule made explicit at the consumer side rather than baked into the artifact.

### 3.4 Splits

The artifact ships with precomputed 5-fold cross-validation assignments at the **family level**, stored in `family_fold`: every gene in the same family is assigned to the same fold. This prevents leakage between folds through paralog membership. Fold IDs are stable across runs of the build pipeline.

### 3.5 Format

```
data/processed/yeast_graph.pt
```

A serialized PyTorch Geometric `Data` object, loadable with:

```python
import torch
data = torch.load("data/processed/yeast_graph.pt", weights_only=False)
```

The object exposes the following fields. `N` is the node count, `E` is the regulatory edge count, `P` is the paralog-pair count.

| Field | Shape | dtype | Description |
|---|---|---|---|
| `x` | `(N, 176)` | `float32` | Continuous per-node attributes packed in a single tensor. Column order is `[in_degree, out_degree, expression_0, ..., expression_172, mean_paralog_identity]`. The exact column meaning is recorded in `feature_names` so downstream code does not have to hard-code positions. |
| `feature_names` | list of length 176 | `str` | Per-column label for `x`. Use to look up which column holds which attribute. |
| `edge_index` | `(2, E)` | `int64` | Regulatory edges in PyG convention (row 0 = source, row 1 = destination). |
| `edge_attr` | `(E, 2)` | `float32` | Per-regulatory-edge attributes: `[min_pvalue, n_conditions_bound]` from the Harbison source. |
| `paralog_edges` | `(2, P)` | `int64` | Paralog pairs as edges (same PyG format as `edge_index`), indexing into the same node set. |
| `paralog_labels` | `(P,)` | `int64` | Per-pair `regulation_retained`. See section 3.3. |
| `family_id` | `(N,)` | `int64` | Integer family ID per gene. Multi-gene Wapinski families occupy `0..M-1`; singletons get unique IDs `M..N-1`. |
| `family_fold` | `(N,)` | `int64` | 5-fold CV assignment via family. Genes whose family contains no paralog pair receive `-1`. |
| `gene_names` | list of length `N` | `str` | Systematic ORF name at each node index. |

`x` and `family_id` are stored separately on purpose: `x` is the standard PyG entry point for continuous features, while `family_id` is categorical and intended to be looked up rather than concatenated. Combining or transforming them further is a downstream consumer concern.

An optional human-readable sidecar `data/processed/node_table.csv` (one row per gene, one column per scalar attribute, omitting the 173 expression columns for readability) can be produced with `python -m src.data.build_graph --node-table` for inspection.

---

## 4. Reproducibility

### 4.1 Software environment

Recorded in [`requirements.txt`](../requirements.txt). Pipeline developed and tested on:

- macOS 15 (darwin 25.4) on Apple Silicon
- Python 3.13.9 (3.14 is too new for `torch` wheels at time of writing)
- `pandas`, `numpy`, `biopython>=1.83`, `xlrd>=2.0`, `openpyxl>=3.1`, `requests`
- `torch>=2.0`, `torch-geometric>=2.4` (required to load `data/processed/yeast_graph.pt`)

To reproduce the environment exactly:

```bash
python3.13 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

### 4.2 End-to-end pipeline

Every script in [`src/data/`](../src/data/) is invoked as a Python module from the repo root:

```bash
# Step 1: download and parse each source
python -m src.data.download_harbison
python -m src.data.parse_harbison

python -m src.data.download_wapinski
python -m src.data.parse_wapinski

python -m src.data.download_gasch
python -m src.data.parse_gasch --long

python -m src.data.download_yeast_proteome
python -m src.data.compute_paralog_similarity

# Step 2: assemble the consolidated graph artifact
python -m src.data.build_graph
# add --node-table to also write the human-readable sidecar CSV
```

Every download script has stable defaults (Wayback snapshots where original hosting is dead, current canonical URLs otherwise). All defaults are overridable with `--url` if a more recent mirror appears.

### 4.3 Sanity-check numbers

After running the pipeline, the following quantities should match. Deviations indicate something has changed in a source file, in parsing logic, or in the build script.

**Per-source outputs (after step 1):**

| Quantity | Expected | Observed |
|---|---|---|
| Distinct TFs in Harbison | ~203 | 204 |
| TF -> target edges (`p <= 0.001`) | ~10,000 | 10,853 |
| Paralog pairs in Wapinski | ~800 (450 WGD) | 801 (438 WGD) |
| Multi-gene families | ~500 | 599 |
| ORFs in Gasch matrix | ~6,200 | 6,152 |
| Conditions in Gasch matrix | ~173 | 173 |
| Paralog pairs aligned | ~800 | 791 (10 missing in SGD proteome) |

**Combined artifact (after step 2):**

| Quantity | Observed |
|---|---|
| Nodes (`N`) | 6,422 |
| Regulatory edges (`E`) | 10,853 |
| Paralog pairs (`P`) | 801 (702 retained, 99 lost) |
| Continuous feature dim (`x.shape[1]`) | 176 (= 2 degree + 173 expression + 1 identity) |
| Unique family IDs | 5,724 (599 multi-gene + 5,125 singletons) |
| Genes assigned to a CV fold | 1,297 (20.2%) |
| Artifact size on disk | ~4.8 MB |

---

## 5. FAIR alignment

This data collection is structured to support the FAIR principles for research data (Wilkinson *et al.* 2016, DOI [10.1038/sdata.2016.18](https://doi.org/10.1038/sdata.2016.18)):

**Findable.** Each dataset is identified by the DOI of its source publication. Persistent identifiers are used throughout (DOIs, SGD systematic ORF names). This document indexes every artifact.

**Accessible.** Every retrieval step is a self-contained, deterministic Python script. Where original hosting is offline (Harbison's Young Lab page, Wapinski's Broad Institute portal, Gasch's Stanford companion site), scripts pull from Internet Archive snapshots whose timestamps are pinned in the script source. Authentication is not required for any source.

**Interoperable.** Every output is CSV with a documented schema. All five datasets share a single gene-identifier convention (systematic ORF names from SGD), which means joins between them are simple equijoins with no ID mapping. The protein FASTA uses the same key. PyTorch Geometric is an open-standard format readable by any PyG-compatible library.

**Reusable.** Provenance (source, retrieval method, processing steps) is recorded per dataset in section 2 and per-script as docstrings. Source data is from peer-reviewed publications; licensing is documented in section 6. The pipeline is reproducible from the commands in section 4.

---

## 6. Licensing and citation

### 6.1 Source licenses

To the best of our knowledge:

- **Harbison 2004** -- published in *Nature*; supplementary text data is widely redistributed for academic research. Verify license at the [article page](https://www.nature.com/articles/nature02800).
- **Wapinski 2007** -- published in *Nature*; supplementary tables are available via the publisher's static-content host. Verify at the [article page](https://www.nature.com/articles/nature06107).
- **Gasch 2000** -- published open-access in *Molecular Biology of the Cell* under an ASCB license; data redistribution is standard.
- **SGD `orf_trans.fasta`** -- distributed by the Saccharomyces Genome Database under [CC-BY 4.0](https://www.yeastgenome.org/help).

Derived files in this directory are released under the project's top-level license.

### 6.2 Citations (BibTeX)

```bibtex
@article{harbison2004,
  author    = {Harbison, Christopher T. and Gordon, D. Benjamin and Lee, Tong Ihn and others},
  title     = {Transcriptional regulatory code of a eukaryotic genome},
  journal   = {Nature},
  volume    = {431},
  pages     = {99--104},
  year      = {2004},
  doi       = {10.1038/nature02800}
}

@article{wapinski2007,
  author    = {Wapinski, Ilan and Pfeffer, Avi and Friedman, Nir and Regev, Aviv},
  title     = {Natural history and evolutionary principles of gene duplication in fungi},
  journal   = {Nature},
  volume    = {449},
  pages     = {54--61},
  year      = {2007},
  doi       = {10.1038/nature06107}
}

@article{gasch2000,
  author    = {Gasch, Audrey P. and Spellman, Paul T. and Kao, Camilla M. and others},
  title     = {Genomic expression programs in the response of yeast cells to environmental changes},
  journal   = {Molecular Biology of the Cell},
  volume    = {11},
  number    = {12},
  pages     = {4241--4257},
  year      = {2000},
  doi       = {10.1091/mbc.11.12.4241}
}

@article{cherry2012sgd,
  author    = {Cherry, J. Michael and Hong, Eurie L. and Amundsen, Craig and others},
  title     = {Saccharomyces Genome Database: the genomics resource of budding yeast},
  journal   = {Nucleic Acids Research},
  volume    = {40},
  number    = {D1},
  pages     = {D700--D705},
  year      = {2012},
  doi       = {10.1093/nar/gkr1029}
}

@article{scruse2026,
  author    = {Scruse, Ashley and Arnold, Jonathan and Robinson, Robert},
  title     = {Counting subnetworks under gene duplication in genetic regulatory networks},
  journal   = {Bulletin of Mathematical Biology},
  volume    = {88},
  pages     = {31},
  year      = {2026},
  doi       = {10.1007/s11538-025-01592-1}
}

@article{wilkinson2016fair,
  author    = {Wilkinson, Mark D. and Dumontier, Michel and Aalbersberg, IJsbrand Jan and others},
  title     = {The {FAIR} {G}uiding {P}rinciples for scientific data management and stewardship},
  journal   = {Scientific Data},
  volume    = {3},
  pages     = {160018},
  year      = {2016},
  doi       = {10.1038/sdata.2016.18}
}
```
