# Example Test Data

Synthetic CSV and metadata files demonstrating expected input formats for each pipeline and script. These files contain fabricated data for testing column detection, parsing, and pipeline integration -- they are **not** real cytometry measurements.

## Files

| File | Used By | Description |
|------|---------|-------------|
| `expression_matrix.csv` | Pipeline B, `load_csv_expression.R` | Pure expression data (markers only) |
| `cytolution_gating_export.csv` | Pipeline C, `load_csv_gating.R` | Cytolution export with all column types: index, scatter, markers, gating, dim-reduction, autofluorescence |
| `conditions.csv` | `differential_abundance.R`, `cluster_enrichment.R` | Sample-to-condition mapping |
| `metadata.csv` | `check_metadata.R` | Clinical metadata with mixed variable types |
| `cluster_tracking_matrix.csv` | `abundance_plotting.R`, `cell_ratio.R` | Cluster counts per sample |
| `manual_gating_labels.csv` | `load_manual_gating.R` | Binary manual gating columns |
| `config_example.R` | All pipelines | Minimal parameter config example |

## Column Types in Cytolution Exports

The `cytolution_gating_export.csv` file demonstrates all 6 column types that `classify_csv_columns()` auto-detects:

1. **Index**: `index` -- row identifier
2. **Scatter/Time**: `SSC.H | SSC.H`, `FSC.A | FSC.A`, `TIME | Time` -- physical parameters
3. **Expression markers**: `CD4 | BV421.A`, `CD8 | APC.A`, etc. -- marker | channel format
4. **Gating/Clustering**: columns ending in `(Clustering annotated)` or `Clustering raw`
5. **Dimensionality reduction**: `cyto_2d_Samples_All_X`, `cyto_2d_Samples_All_Y`
6. **Autofluorescence**: `AF.P1 | AF.P1.A` -- channels without a meaningful marker name

## Usage

```r
# Test column classification
source("R_shared/io_utils.R")
header <- names(data.table::fread("test_data/examples/cytolution_gating_export.csv", nrows = 0))
classes <- classify_csv_columns(header)
str(classes)

# Test expression loading
source("loading/load_csv_expression.R")
data <- load_csv_expression("test_data/examples/")
```
