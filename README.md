# xeno-allo-pbmc

Analysis code for a single-cell transcriptomic comparison of peripheral blood
mononuclear cells (PBMCs) between pig-to-human kidney xenotransplantation and
human kidney allograft rejection.

This repository contains the notebooks that produce the main and supplementary
figures of the accompanying manuscript. It is intended to document exactly how
each panel was generated, not to be a one-command reproducible pipeline: the
input objects are derived from proprietary sequencing data and are not
distributed here.

---

## Repository layout

| File | Panels | Content |
|---|---|---|
| `01_Fig1_S1.ipynb` | Fig. 1A–D, Fig. S1A–D | Quality control, filtering, normalisation, HVG selection, PCA, Harmony batch integration, Leiden clustering, major lineage annotation (monocytes/DC/megakaryocytes, T/NK, B cells), marker dot plots and cell-type composition |
| `02_Fig2.ipynb` | Fig. 2A–E | Monocyte subclustering, marker definition of monocyte states, differential composition testing (scCODA and DCT), subset-level marker summaries |
| `03_Fig3.ipynb` | Fig. 3A–E | Time-resolved analysis of the xenograft monocyte compartment: pseudobulk DESeq2 (`ov.bulk.deseq`), differential expression (`ov.bulk.pyDEG`) and temporal (fuzzy c-means) gene clustering (`pymfuzz`, `ov.bulk.temporal_clusters`) |
| `04_Fig4_S3.ipynb` | Fig. 4A–F, Fig. S3A–D | Trajectory inference (PAGA, `TrajInfer`), dynamic gene trends along pseudotime, SCENIC regulon inference and CellOracle in silico transcription-factor perturbation |
| `05_Fig5_S4.ipynb` | Fig. 5A–F, Fig. S4A–H | Cell–cell communication: CellPhoneDB run through `omicverse`, communication heatmaps and networks, `ktplotspy` visualisation |
| `06_Fig6.ipynb` | Fig. 6A–B | NicheNet ligand–activity analysis of the PF4⁺PPBP⁺ monocyte state, ligand–receptor signalling and target-gene programmes |
| `02_Fig2F.R` | Fig. 2F | Volcano plot of differentially expressed genes in PF4⁺ monocytes versus the remaining monocytes (`ggplot2`/`ggrepel`) |

Notebooks are numbered in figure order. Cells are annotated with the panel they
generate (for example `# Figure 2D`), so a given panel can be located by
searching the notebook for its label.

All notebooks were run with Python 3.10 and the packages listed in
`requirements.txt`; `02_Fig2F.R` was run with R 4.x and the CRAN packages
`tidyverse`, `RColorBrewer` and `ggrepel`.

---

## Input data

The analysis starts from AnnData objects. These files are **not** included in
this repository; they are large and derived from controlled-access data.
The objects expected by each notebook are:

| File | Used by |
|---|---|
| `20260909GSE2Xeno_GK.h5ad` | `01` |
| `20260909AllO_PBMC.h5ad` | `01` |
| `20260911Allo_Xeno_PBMC.h5ad` | `02`, `03`, `04`, `05` |
| `20260917Xeno_Allo.h5ad` | `05` |
| `adata_mono.h5ad` | `02` onwards (monocyte subset object) |
| `adata_xeno_pf4_major.h5ad`, `adata_allo_pf4_major.h5ad` | `06` |
| `PF4_Mono_vs_rest_DEGs_clean.csv` | `02_Fig2F.R` |
| `CellOracle_TF_KO_Ranking.csv`, `CellOracle_9TF_AllSubtypes_*.csv` | `04` (intermediate tables) |
| `GPL30305_family.soft.gz` | `06` (public GEO platform annotation, used to map probe IDs) |

Two external resources are also required by `06`: the NicheNet model files
`ligand_target_matrix_human.bin`, `lr_network_human.bin` and `lr_sig_human.bin`,
available from the NicheNet distribution.

> **Data availability.** The raw sequencing data supporting this study are
> deposited in _[accession to be added]_. Accession numbers and the conditions of
> use are described in the manuscript.

---

## Execution order and state dependencies

The notebooks are not independent. `02` reads its input object directly, but
`03` continues from objects (`adata_mono`, `adata_pf4`) left in memory by `02`,
and `04` expects a monocyte object loaded into `adata_mono2`. In practice they
were run as one sequence: `01` → `02` → `03` → `04`, with `05` and `06` run as
separate sessions after their input objects had been written to disk.

| Notebook | Starts from |
|---|---|
| `01_Fig1_S1` | reads the two raw `h5ad` files (see below) |
| `02_Fig2` | reads `20260911Allo_Xeno_PBMC.h5ad`; leaves `adata_mono` in memory |
| `03_Fig3` | continues from `02` (`adata_mono`) |
| `04_Fig4_S3` | reads the monocyte object into `adata_mono2` |
| `05_Fig5_S4` | reads `20260917Xeno_Allo.h5ad` and `adata_mono.h5ad` |
| `06_Fig6` | reads `adata_xeno_pf4_major.h5ad` and `adata_allo_pf4_major.h5ad` |

## Required edits before running

The notebooks were developed on an internal server and correspond exactly to
what was executed. Absolute and working-directory paths were left as they were,
and several working-directory strings were blanked to empty strings (`""`)
before publication. All of them must be set before the notebooks will run
elsewhere.

| Literal | Where | Meaning |
|---|---|---|
| `os.chdir("")` | `01`, `02`, `04`, `05` | working directory |
| `save_dir=""` (default of `save_fig`) | `01`, `02`, `04`, `05` | figure output directory |
| `SAVE_DIR = ""` | `04` (CellOracle cells), `03` (`/20260922`) | figure/table output directory |
| `/20260909GSE2Xeno_GK.h5ad`, `/20260909AllO_PBMC.h5ad` | `01` | raw objects, at filesystem root |
| `/20260911Allo_Xeno_PBMC.h5ad` | `02` | merged object, at filesystem root |
| `ov.read("")` | `04` | monocyte object |
| `20260917Xeno_Allo.h5ad`, `adata_mono.h5ad` (relative) | `05` | input objects |
| `scenic/cellphonedb.zip`, `/scenic/cellphonedb.zip` | `05` | CellPhoneDB database archive |
| `/data/zdh/ZDH_Spatial/Python包` | `04`, `06` | SCENIC resources (hg38 motif databases) and NicheNet model binaries |
| `/data/zdh/ZDH_Spatial/Kidney_Spatial_Project/ceshi/Figure/20260922` | `06`, `02_Fig2F.R` | working/figure output directory |

The simplest approach is to define one variable at the top of each notebook
(for example `DATA_DIR = "/path/to/data"`) and replace the literals below it.
`save_fig` and `SAVE_DIR` control where figures and intermediate tables are
written; setting those two strings redirects almost all output.

---

## Notes

- Notebooks are stored with outputs cleared, so diffs stay readable. Rendering
  a notebook on GitHub will therefore show code only.
- The analysis was performed with several packages under active development
  (`omicverse`, `pySTIM`, `pymfuzz`, `nichenetpy`, `CellOracle`). Behaviour and
  defaults may differ between versions; pin the versions if exact numerical
  reproduction matters.
- Random seeds are set explicitly where stochastic steps are used
  (for example `rng_key=42` for scCODA, `np.random.seed` in `04`).

---

## Citation

If you use this code, please cite the accompanying manuscript:

> _[authors]. [title]. [journal], [year]. doi:[DOI]_

## Contact

_[name], [affiliation]_ — _[email]_

Correspondence and questions about the analysis are welcome via GitHub issues.

## License

No license has been granted for this repository, which means the default
copyright applies: the code is publicly readable but may not be reused,
modified or redistributed without permission. If you would like to reuse it,
please get in touch.
