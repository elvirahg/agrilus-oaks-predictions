# Predicting future <em>Agrilus</em> – <em>Quercus</em> interactions

This repository contains code and analyses associated with Hernandez-Gutierrez et al. (2024), “Combined phylogenetic and geographic data can predict plant–pest interactions with high accuracy” (DOI: 10.1101/2024.10.30.619026). These analyses explore how readily available phylogenetic and geographical information can be used to predict previously unobserved host–pest interactions, through relatively simple multilevel Bayesian models. Our aim is to implement a broadly applicable approach for identifying potentially harmful host – pest interactions, and help to prioritise counter measures against threats worldwide.

## Repository structure and contents

Please note that the original code used to generate the data in the paper is under the `original-analyses` branch. The current branch mirrors these same analyses, but with refactored code for better readability, reusability, and efficiency.

* `analysis/`: R scripts implementing the analysis workflow. The scripts are numbered to indicate execution order (e.g. `02_*` follows `01_*`). These can be used to reproduce the published analyses or adapted for related work.
* `R/`: helper functions used across the analyses.
* `docs/`: HTML documentation for the functions in R/. If you clone or download the repository locally, you can open index.html to browse the documentation.
* `data/`: input data required for the analyses.
