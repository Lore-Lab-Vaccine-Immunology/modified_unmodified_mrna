# Immunogenicity and Innate Immunity to High-Dose and Repeated Vaccination of Modified mRNA versus Unmodified mRNA
> Engstrand, Olivia et al. (2025) *Molecular Therapy – Nucleic Acids*
> 
This repository contains the analysis scripts used in the study.


## Contents

- `scripts/` – All R and Quarto scripts used for data analysis
- `data/` – Expected data input folder
- `output/` – Results generated from the scripts
- `renv/` and `renv.lock` – Environment for package reproducibility
- `.Rproj` – R project file for ease of use in RStudio


This project uses [`renv`](https://rstudio.github.io/renv/) to ensure reproducible environments.
To restore the exact R package versions used in the analysis:

```r
renv::restore()
