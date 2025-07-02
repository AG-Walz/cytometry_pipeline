# Test Data

This folder contains sample FCS files and associated data for testing the cytometry pipeline.

## Contents

| File | Description |
|------|-------------|
| `ICK-049-01_Unmixed.fcs` | Unmixed spectral flow cytometry data (sample 1) |
| `ICK-049-03_Unmixed.fcs` | Unmixed spectral flow cytometry data (sample 3) |
| `ICK-049-051.wsp` | FlowJo workspace file with gating strategy |
| `ICK-049-01_Unmixed.fcs_2128968_populations.csv` | Population labels for sample 1 |
| `ICK-049-03_Unmixed.fcs_1445584_populations.csv` | Population labels for sample 3 |

## Usage

These files can be used to test:
- FCS file loading (`Load_fcs.R`, `Load_raw.R`)
- CSV gating import (`Load_csv_gating.R`)
- FlowJo workspace parsing (`Load_wsp.R`, `FlowJo_wsp.R`)
- UMAP visualization
- Gating comparison workflows

## Data Source

ICK study samples (anonymized spectral flow cytometry data).
