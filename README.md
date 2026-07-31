# ML4GutHealth

### A Data-Driven Universal Gut Microbiome Health Assessment: A Machine Learning Framework Trained on Large Metagenomic Data

**Bablu Kumar**¹ ², **Erika Lorusso**² ³, **Bruno Fosso**² \*, **Graziano Pesole**² ³ \*

¹ Department of Oncology and Hematology-Oncology, Università degli Studi di Milano, Milan, Italy
² Department of Biosciences, Biotechnology and Environment, University of Bari A. Moro, Bari, Italy
³ National Research Council, Institute of Biomembranes, Bioenergetics and Molecular Biotechnologies, Bari, Italy

\* Corresponding authors — Graziano Pesole (graziano.pesole@uniba.it, graziano.pesole@cnr.it), Bruno Fosso (bruno.fosso@uniba.it)

---

## About

This repository contains the code and analysis scripts for a machine learning classifier that distinguishes healthy from non-healthy individuals using gut microbiome species-level profiles. The classifier was trained on 7,452 publicly available stool shotgun metagenomes from 32 independent studies, covering 12 disease conditions and healthy controls, and externally validated on 642 samples from 6 independent studies, including disease types not present in the training data.

---

## Repository Structure

```
ML4GutHealth/
├── preprocessing/
│   └── computational_pipeline.sh
│
├── diversity_analysis/
│   ├── alpha_diversity/
│   │   └── alpha_diversity_calculation.R
│   └── beta_diversity/
│
├── batch_correction/
│
├── model_training/
│   ├── rf_model_training.py
│   ├── svm_model_training.py
│   └── LR_ElasticNet_training.py
│
├── sensitivity_analysis/
│   ├── normalization_strategy/     # Raw counts, TSS, and CLR comparison
│   └── clr_pseudocount/            # CLR pseudocount sensitivity (1.0, 0.5, 0.1)
│
└── figures/                        # Scripts used to generate manuscript figures
```

---

## Pipeline Overview

| Step | Description | Tool(s) |
|---|---|---|
| 1 | Read quality trimming and adapter removal | Fastp (v0.24.0) |
| 2 | Host-read removal (GRCh38, no-alt) | Bowtie2 (v2.5.4) |
| 3 | Taxonomic classification and abundance estimation | Kraken2 (v2.1.3), Bracken (v2.9) |
| 4 | Batch correction (BioProject as batch variable) | MMUPHin |
| 5 | Normalization | Centered log-ratio (CLR) transformation |
| 6 | Model training (4 algorithms × 3 feature-selection strategies) | scikit-learn |
| 7 | Feature selection | Permutation feature importance (PFI) |
| 8 | Interpretability | MaAsLin2 |

The Kraken2 database used for taxonomic classification is built from the UHGG human-gut reference (v2.0.2).

---

## Data

Species-level abundance matrices, sample metadata, and PFI-selected feature lists are described in the manuscript's Data Availability section. Raw sequencing data were obtained from public repositories (SRA/ENA) under the accessions listed in the manuscript's supplementary materials.

---

## Requirements

- **Python** 3.10 — conda environment specification in `environment.yml`
- **R** — version specified in `renv.lock`; required packages: MMUPHin, MaAsLin2, iNEXT, ggplot2, patchwork
- **Key Python packages**: scikit-learn, scikit-bio, pandas, numpy

---

## Usage

Each subdirectory contains scripts intended to be run in the order listed in its own README, or as indicated by numbered filenames. Input paths are set at the top of each script and should be adjusted to match the local data location.

---

## Citation

If you use this code, please cite the associated manuscript (citation to be added upon publication).

---

## Contact

For questions regarding this repository, please open an issue or contact the corresponding authors:

- Graziano Pesole — graziano.pesole@uniba.it, graziano.pesole@cnr.it
- Bruno Fosso — bruno.fosso@uniba.it
