# Changelog

## 0.5.0
- **Saved runs.** Every analysis is saved automatically (inputs, settings, a note you can edit) and can be reopened from the new Runs page. Results are recomputed on opening, so nothing is stale. Rename, annotate and delete runs.
- **Raw spectrum check (Database check mode).** Reads Thermo `.raw` files in place (no conversion, via ThermoRawFileParser's query and xic modes) and mzML. For each peptide it finds the precursor in MS1, takes the fragmentation scans at its elution peak, matches b and y ions, and scores them with a chance probability corrected for choosing the best scan and checked against 20 shuffled decoy peptides of the same mass. Each peptide is compared with its counterpart carrying the other residue. Annotated spectra, MS1 chromatograms and fragment co-elution plots. A "check any peptide" box works on any sequence.
- **Detectability.** For every changed position the program works out which tryptic peptides each version would make and whether a search could see them, so "no coverage" can be interpreted.
- **QC.** A QC page for the database check (digestion, charge states, lengths, single-peptide proteins, run-to-run spread, proteotypic share) and a checks panel with per-sample table and outlier flags for quantification. Thresholds are shown next to each result.
- **Export every table on its own** as CSV, TSV or Excel, from a button on each table or from the export page, plus a ZIP of everything. New tables: QC checks, per-run and per-sample QC, spectrum checks, matched fragment ions, saved runs.
- DIA-NN long-format reports (`report.tsv`) are read too, giving retention times and q-values.
- The PDF report gains a quality-checks section, a detectability explanation and a spectrum-evidence section with annotated spectra.
- Simulated examples are now identical on every run (they depended on Python's set ordering before).

## 0.4.0
Renamed from Proteomics Analyzer to **Recode Detector**. New look: "Paper Lab" (cream paper, ink blue and vermilion, serif headings, an animated ridge-line plot) is the default; nine other looks (Ultraviolet, Solar Flare, Neon Grid, Mycelium, Black Box, Rose Quartz, Copper Helix, Sky Lab, Spectrum Sunset) are in the "Look" menu in the header. The Python package is now `recode_detector`, the port variable `RD_PORT`, the Mac app `Recode Detector.app`.

## 0.3.1
The PDF report always has a section on the standard-database search: a figure, a table of positions and a plain-English reading when a standard-search report is loaded, and a "not included" note with instructions when it is not.

## 0.3.0
Database check mode: a **Coverage** page (proteins detected, sequence coverage, detection by protein size, peptide length and charge, per-protein table) and a **plain-English PDF report** (short version, what was measured, coverage, the question, results table, caveats, next steps, glossary, methods). New dependency: reportlab; the Mac app installs it itself on first start.

## 0.2.1
The message shown when a DIA-NN precursor matrix is loaded on the Quantification page now points to Database check mode.

## 0.2.0
**Database check mode** (switch at the top of the page): which version of a protein did a search see? Load a peptide report and a reference FASTA, plus optionally an alternative-coding FASTA (or one combined FASTA whose variant entries share an ID marker) and the report from a standard-database search. The program pairs each variant with its reference protein, finds every differing residue, and lists the identified peptides that span it, with intensities, missed cleavages and an ambiguity check, and compares with the standard-database search. CSV export of sites and peptides.

## 0.1.1
DIA-NN: descriptive columns (Proteotypic, Precursor.Charge, sequences) are no longer mistaken for samples; a precursor matrix (`pr_matrix`) is refused with advice to load `pg_matrix`; a table with a single run explains that replicates are needed. The macOS app carries its own copy of the code (fixes it not opening).

## 0.1.0
First version: MaxQuant / DIA-NN / generic table import, group assignment, filtering, median and quantile normalisation, imputation, QC plots, moderated and Welch t-tests with BH FDR, volcano, heatmap, gene-set enrichment from `.gmt`, CSV and settings export. Quits completely when its window is closed.
