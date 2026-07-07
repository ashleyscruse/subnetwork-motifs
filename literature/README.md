# Literature

Download PDFs to this folder for local reference. PDFs are gitignored (copyrighted, large files). This README is the reading list.

## Must-Read

These papers are directly used in the project. Read them before writing any code.

1. **Scruse, Arnold & Robinson 2026** -- "Counting Subnetworks Under Gene Duplication in Genetic Regulatory Networks." *Bulletin of Mathematical Biology*, 88:31.
   - [Open access](https://link.springer.com/article/10.1007/s11538-025-01592-1)
   - The published framework. The math this entire project implements. Start here.

2. **Harbison et al. 2004** -- "Transcriptional regulatory code of a eukaryotic genome." *Nature*, 431:99-104.
   - [DOI: 10.1038/nature02800](https://doi.org/10.1038/nature02800)
   - The yeast TF-target binding data. 203 TFs mapped to ~6,000 targets. This is the regulatory network.

3. **Wapinski et al. 2007** -- "Natural history and evolutionary principles of gene duplication in fungi." *Nature*, 449:54-61.
   - [DOI: 10.1038/nature06107](https://doi.org/10.1038/nature06107)
   - ~450 duplicate gene pairs with regulatory fate classifications across 17 fungal genomes. These are the training labels.

4. **Milo et al. 2002** -- "Network motifs: simple building blocks of complex networks." *Science*, 298:824-827.
   - [DOI: 10.1126/science.298.5594.824](https://doi.org/10.1126/science.298.5594.824)
   - The original network motifs paper. Defines the concept that subnetwork motifs improve upon.

5. **Cordero & Hogeweg 2006** -- "Feed-forward loop circuits as a side effect of genome evolution." *Molecular Biology and Evolution*, 23:1931-1936.
   - [DOI: 10.1093/molbev/msl060](https://doi.org/10.1093/molbev/msl060)
   - Showed feedforward loops can arise from duplication, not selection. The biological motivation for this project.

## Should-Read

Context and methods. Read these to understand the tools and approaches.

6. **Teichmann & Babu 2004** -- "Gene regulatory network growth by duplication." *Nature Genetics*, 36:492-496.
   - [DOI: 10.1038/ng1340](https://doi.org/10.1038/ng1340)
   - Gene duplication as the dominant mechanism of regulatory network growth.

7. **Velickovic et al. 2018** -- "Graph Attention Networks." *ICLR 2018*.
   - [arXiv: 1710.10903](https://arxiv.org/abs/1710.10903)
   - The GAT architecture we use. Read sections 1-3 for the core idea.

8. **Kipf & Welling 2017** -- "Semi-Supervised Classification with Graph Convolutional Networks." *ICLR 2017*.
   - [arXiv: 1609.02907](https://arxiv.org/abs/1609.02907)
   - GCN, the fallback architecture. Simpler than GAT, no attention mechanism.

9. **MacIsaac et al. 2006** -- "An improved map of conserved regulatory sites for Saccharomyces cerevisiae." *BMC Bioinformatics*, 7:113.
   - [DOI: 10.1186/1471-2105-7-113](https://doi.org/10.1186/1471-2105-7-113)
   - Reanalysis and refinement of the Harbison regulatory network.

10. **Wernicke & Rasche 2006** -- "FANMOD: a tool for fast network motif detection." *Bioinformatics*, 22:1152-1153.
    - [DOI: 10.1093/bioinformatics/btl038](https://doi.org/10.1093/bioinformatics/btl038)
    - The tool we benchmark against. Understand what it does and why its null model is limited.

## Nice-to-Read

Deeper context. Not required to start working but useful for the paper and discussion.

11. **Kashtan et al. 2004** -- "Efficient sampling algorithm for estimating subgraph concentrations and detecting network motifs." *Bioinformatics*, 20:1746-1758.
    - [DOI: 10.1093/bioinformatics/bth163](https://doi.org/10.1093/bioinformatics/bth163)
    - mfinder, the other motif detection tool we benchmark against.

12. **Gasch et al. 2000** -- "Genomic expression programs in the response of yeast cells to environmental changes." *Molecular Biology of the Cell*, 11:4241-4257.
    - [DOI: 10.1091/mbc.11.12.4241](https://doi.org/10.1091/mbc.11.12.4241)
    - Yeast stress response expression data. Used as optional node features.

13. **Force et al. 1999** -- "Preservation of duplicate genes by complementary, degenerative mutations." *Genetics*, 151:1531-1545.
    - [DOI: 10.1093/genetics/151.4.1531](https://doi.org/10.1093/genetics/151.4.1531)
    - The DDC (duplication-degeneration-complementation) model. Theoretical background for why duplicates are retained.
