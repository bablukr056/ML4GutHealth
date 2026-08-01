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
├── Computational_pipeline/     # Raw read download and preprocessing (QC, host removal, taxonomic classification)
├── Data/                       # Species-level abundance matrices and sample metadata
├── R_Code_Analysis/            # Batch correction, diversity analyses, and figure-generation scripts (R)
├── ml_model_training_codes/    # Model training scripts (Random Forest, SVM, Logistic Regression ElasticNet)
├── Manuscript_Figures/         # Main text and supplementary figures
```

---

## Pipeline Overview

![Pipeline Overview](https://raw.githubusercontent.com/bablukr056/ML4GutHealth/main/Figure_1.jpeg)
---

## Citation

If you use this code, please cite the associated manuscript (citation to be added upon publication).

---

## Contact

For questions regarding this repository, please open an issue or contact the corresponding authors:

- Graziano Pesole — graziano.pesole@uniba.it, graziano.pesole@cnr.it
- Bruno Fosso — bruno.fosso@uniba.it
