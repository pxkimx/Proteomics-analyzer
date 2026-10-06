# Changelog

## 0.2.1
The message shown when a DIA-NN precursor matrix is loaded on the Quantification page now points to Database check mode.

## 0.2.0
**Database check mode** (switch at the top of the page): which version of a protein did a search see? Load a peptide report and a reference FASTA, plus optionally an alternative-coding FASTA (or one combined FASTA whose variant entries share an ID marker) and the report from a standard-database search. The program pairs each variant with its reference protein, finds every differing residue, and lists the identified peptides that span it, with intensities, missed cleavages and an ambiguity check, and compares with the standard-database search. CSV export of sites and peptides.

## 0.1.1
DIA-NN: descriptive columns (Proteotypic, Precursor.Charge, sequences) are no longer mistaken for samples; a precursor matrix (`pr_matrix`) is refused with advice to load `pg_matrix`; a table with a single run explains that replicates are needed. The macOS app carries its own copy of the code (fixes it not opening).

## 0.1.0
First version: MaxQuant / DIA-NN / generic table import, group assignment, filtering, median and quantile normalisation, imputation, QC plots, moderated and Welch t-tests with BH FDR, volcano, heatmap, gene-set enrichment from `.gmt`, CSV and settings export. Quits completely when its window is closed.
