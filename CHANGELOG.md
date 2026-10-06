# Changelog

## 0.1.1
DIA-NN: descriptive columns (Proteotypic, Precursor.Charge, sequences) are no longer mistaken for samples; a precursor matrix (`pr_matrix`) is refused with advice to load `pg_matrix`; a table with a single run explains that replicates are needed. The macOS app carries its own copy of the code (fixes it not opening).

## 0.1.0
First version: MaxQuant / DIA-NN / generic table import, group assignment, filtering, median and quantile normalisation, imputation, QC plots, moderated and Welch t-tests with BH FDR, volcano, heatmap, gene-set enrichment from `.gmt`, CSV and settings export. Quits completely when its window is closed.
