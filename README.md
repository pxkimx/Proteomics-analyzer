# Recode Detector

*(formerly Proteomics Analyzer)*

A local, offline program for proteomics. Its headline job is to **detect genetic recoding**: when the same data are searched against a standard protein database and an alternative one (for example a codon read as a different amino acid), it shows which version of each protein the data support, with coverage statistics and a plain-English PDF report. It also does the protein-level part of an ordinary quantitative experiment: load a quantification table, check its quality, normalise, impute, find differentially abundant proteins, and test gene-set enrichment. It runs on your own computer; nothing is uploaded anywhere.

It starts after the search engine. Use MaxQuant, DIA-NN, Spectronaut or similar to identify and quantify proteins, then bring the protein table here.

## What it does

| Step | |
|---|---|
| **Data** | MaxQuant `proteinGroups.txt` (LFQ, raw intensity or iBAQ; reverse hits, contaminants and site-only hits are removed), DIA-NN `pg_matrix`, or any table with one row per protein and one numeric column per sample (TSV, CSV, Excel). Drag and drop or browse. |
| **Groups** | Replicate groups are guessed from sample names (`Ctrl_1`, `Ctrl_2`, `KO rep3` …) and are editable. Leave a group blank to exclude a sample. |
| **Process** | log2 transform, valid-value filter (per group), median or quantile normalisation, imputation (Perseus-style down-shifted normal, half-minimum, or none). Every choice is recorded. |
| **Quality** | Proteins per sample, intensity distributions before/after normalisation, PCA, sample correlation (clustered), CV within groups, dynamic range. |
| **Differential** | Moderated t-test (limma-style empirical Bayes) or Welch t-test, Benjamini–Hochberg FDR, adjustable FDR and fold-change cut-offs. Interactive volcano plot, per-protein dot plot with imputed values marked, clustered heatmap, sortable table. |
| **Enrichment** | Over-representation (hypergeometric) of significant proteins in gene sets from a `.gmt` file you provide, against the proteins quantified in the comparison as background. |
| **Export** | Results CSV, processed matrix CSV, and a JSON file with the groups, parameters and processing steps for your methods section. |

## Database check mode

A second mode (switch at the top of the page) for a different question: the same data were searched against a reference protein database and an alternative one, for example a genome translated with a codon read differently (CGG as Trp instead of Arg, an alternative start site, a point-mutation set). **Which version of each protein did the search actually see?**

You load the peptide report (DIA-NN `report.pr_matrix.tsv`), the reference FASTA, and, optionally, the alternative FASTA. If one FASTA holds both standard and variant proteins, give the text that marks the variant IDs. The program pairs each variant with its reference protein, finds every position where they differ, and lists the identified peptides that span that position. Only such a peptide can tell the two versions apart. Peptides that also occur elsewhere in the database are shown as ambiguous and not counted. A report from a search against the standard database alone can be added for comparison.

A **Coverage** page reports how much of the reference proteome the run saw: proteins detected, how many with two or more peptides, the share of all amino acids covered, detection by protein size, and a per-protein table. The **Report** page downloads a **PDF written in plain English** (summary, what was measured, coverage, the question, results, caveats, next steps, glossary, methods) and CSVs of sites, peptides and coverage.

Each site gets a verdict: *alternative seen*, *reference seen*, *both seen*, or *no coverage*. Treat a verdict as a lead to check against the spectra, not as proof: it is only as reliable as the search's own false-discovery control, and a single peptide is thin evidence. *No coverage* is not evidence against either version.

## Run it

Needs Python 3.9 or newer.

```bash
git clone https://github.com/<you>/Recode-Detector.git
cd Recode-Detector
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
python -m recode_detector          # opens http://localhost:8775
```

Closing the browser tab stops the program. (A page reload does not.)

**macOS app:** `macos/build_app.sh` builds a double-clickable `Recode Detector.app` with an icon. On first start it creates its own Python environment under `~/Library/Application Support/RecodeDetector`. The app carries its own copy of the code, so rebuild it after you change the code.

No data of your own yet? Use **Load simulated example** on the first page: a MaxQuant-style table with 3,000 invented proteins, two groups of four runs, left-censored missing values and flagged contaminants. It is simulated, not real biology, and the true changes are known, which the test suite uses.

## Methods notes

- **Missing values.** Zeros and blanks are treated as "not quantified". The default imputation draws from a normal distribution shifted 1.8 SD below each sample's observed distribution (width 0.3 SD), the convention popularised by Perseus. It assumes missing values are mostly low-abundance proteins below detection. Imputed values are flagged in the dot plot, and the table reports how many values per group were actually observed. Choose "observed values only" in the comparison to test without imputed values.
- **Moderated t-test.** Per-protein variances are shrunk toward a common prior following Smyth (2004), as in limma, with the prior degrees of freedom estimated from the data and shown on the page. For two-group comparisons it matches limma's `eBayes` approach; it is not a full linear-model framework (no covariates, paired designs or multi-group contrasts yet).
- **Fold change** is group A minus group B on the log2 scale, so positive means higher in A.
- **Multiple testing** is Benjamini–Hochberg across all tested proteins in the comparison.
- The valid-value filter, fold-change cut-off and FDR are analysis choices, not defaults of nature: set them for your experiment and report them. The exported settings file does this for you.

## Tests

```bash
pip install pytest
pytest
```

The tests check the file readers, normalisation and imputation, the statistics against simulated truth (sensitivity, false-discovery proportion, calibration under the null, equality with SciPy's Welch test), the HTTP API, and that the program quits when its window closes.

## Not yet included

Paired and multi-factor designs, database check from long-format DIA-NN reports (only the precursor matrix is read), ANOVA across more than two groups, phospho-site level analysis, built-in annotation databases (gene sets are supplied by you), and PTM or peptide-level views.

## License

MIT. See `LICENSE`.
