# Proteomics Analyzer

A local, offline program for the protein-level part of a proteomics experiment: load a quantification table, check its quality, normalise, impute, find differentially abundant proteins, and test gene-set enrichment. It runs on your own computer; nothing is uploaded anywhere.

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

## Run it

Needs Python 3.9 or newer.

```bash
git clone https://github.com/<you>/Proteomics-analyzer.git
cd Proteomics-analyzer
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
python -m proteomics_analyzer          # opens http://localhost:8775
```

Closing the browser tab stops the program. (A page reload does not.)

**macOS app:** `macos/build_app.sh` builds a double-clickable `Proteomics Analyzer.app` with an icon. On first start it creates its own Python environment under `~/Library/Application Support/ProteomicsAnalyzer`. The app points at this folder, so rebuild it if you move the folder.

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

Paired and multi-factor designs, ANOVA across more than two groups, phospho-site level analysis, built-in annotation databases (gene sets are supplied by you), and PTM or peptide-level views.

## License

MIT. See `LICENSE`.
