# DNA Methylation 450K Analysis Pipeline

An R pipeline for processing and analyzing Illumina Infinium HumanMethylation450 (450K) array data using the `minfi` Bioconductor package. The pipeline takes raw IDAT files through quality control, normalization, and differential methylation analysis between two sample groups (control vs. disease).

## What it does

1. **Data import** — reads a sample sheet and loads raw IDAT files into an `RGChannelSet`
2. **Raw signal inspection** — extracts Red/Green fluorescence intensities and checks a specific probe address against the Illumina manifest (Type I/II design, color channel)
3. **Quality control** — computes median intensities (`getQC`), plots a QC summary, and checks negative control probes
4. **Detection p-values** — flags probes that fail a detection p-value threshold (default `0.01`)
5. **Raw signal distributions** — compares mean Beta/M-values between groups before normalization
6. **Normalization** — SWAN normalization (`preprocessSWAN`), with before/after comparison split by probe chemistry (Type I vs Type II)
7. **Exploratory analysis** — PCA on normalized Beta values, colored by group, sex, and batch (Sentrix ID) to check for confounders
8. **Differential methylation testing** — `dmpFinder` for nominal p-values, with Benjamini-Hochberg and Bonferroni correction
9. **Visualization** — volcano plot (effect size vs. significance), Manhattan plot (genomic position vs. significance) via `qqman`, and a hierarchically clustered heatmap of the top differentially methylated probes

## Requirements

- R (≥ 4.0 recommended)
- Bioconductor packages: `minfi`, `qqman` (via CRAN), `gplots`

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("minfi")

install.packages(c("qqman", "gplots"))
```

## Data

A small demo dataset (8 samples, WT vs. MUT) is included in [`data/Input_Data/`](data/Input_Data/) so the pipeline can be run out of the box — see [`data/README.md`](data/README.md) for details.

You'll additionally need a cleaned Illumina 450K manifest file, `Illumina450Manifest_clean.RData`, with probe annotation columns `IlmnID`, `AddressA_ID`, `AddressB_ID`, `Infinium_Design_Type`, `Color_Channel`, `CHR`, `MAPINFO`. This isn't included in the repo (it's a large derived file) — it can be built from Illumina's publicly available [450K manifest files](https://support.illumina.com/downloads/infinium_humanmethylation450_product_files.html). Place it in the project root, or update the `load(...)` path in the script.

## Usage

```r
source("scripts/methylation_pipeline.R")
```

Or open it as an R Markdown document if you're working from the `.Rmd` source, and knit it to produce a full report with all figures inline.

## Example output

See the [`figures/`](figures/) folder for example plots produced by the pipeline (PCA, volcano plot, Manhattan plot, heatmap of top differentially methylated probes).

## License

This project is licensed under the MIT License — see [LICENSE](LICENSE) for details.

