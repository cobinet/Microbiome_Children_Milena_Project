# Data Directory

**Raw data and result files are NOT included in this repository.**

This directory contains only this README describing the expected data files,
their formats, and how to obtain them for replication.

---

## Required Input Files

All files should be placed in the directory structure shown below before running
the analysis scripts.

```
data/
├── ASV_table.tsv
├── ASV_tax_species.user.tsv
├── EC_pred_metagenome_unstrat_descrip.tsv
├── KO_pred_metagenome_unstrat_descrip.tsv
├── METACYC_path_abun_unstrat_descrip.tsv
└── data_formicrobiome.xlsx
```

---

## File Specifications

### `ASV_table.tsv`
| Property | Value |
|---|---|
| Format | Tab-separated values |
| Encoding | UTF-8 |
| Rows | ASV identifiers (e.g. `ASV0001`) |
| Columns | Sample identifiers (`PC_code`, e.g. `PC001`) |
| Values | Relative abundance (non-negative floats) |
| First column header | `ASV_ID` |

### `ASV_tax_species.user.tsv`
| Property | Value |
|---|---|
| Format | Tab-separated values |
| Encoding | UTF-8 |
| Columns | `ASV_ID`, `Kingdom`, `Phylum`, `Class`, `Order`, `Family`, `Genus`, `Species` |
| Notes | Output from QIIME2 taxonomy classifier |

### `EC_pred_metagenome_unstrat_descrip.tsv`
| Property | Value |
|---|---|
| Format | Tab-separated values |
| Encoding | UTF-8 |
| First column | `function` — EC numbers (format: `ec:1.1.1.1`) |
| Second column | `description` — enzyme name |
| Remaining columns | Sample identifiers (`PC_code`), values = predicted abundance |
| Source | PICRUSt2 v2.5+ `EC_metagenome_out/pred_metagenome_unstrat_described.tsv` |

### `KO_pred_metagenome_unstrat_descrip.tsv`
| Property | Value |
|---|---|
| Format | Tab-separated values |
| Encoding | UTF-8 |
| First column | `function` — KEGG Orthology IDs (format: `K00001`) |
| Second column | `description` — KO name |
| Remaining columns | Sample identifiers, values = predicted abundance |
| Source | PICRUSt2 v2.5+ `KO_metagenome_out/pred_metagenome_unstrat_described.tsv` |

### `METACYC_path_abun_unstrat_descrip.tsv`
| Property | Value |
|---|---|
| Format | Tab-separated values |
| Encoding | UTF-8 |
| First column | `pathway` — MetaCyc pathway IDs |
| Second column | `description` — pathway description |
| Remaining columns | Sample identifiers, values = predicted abundance |
| Source | PICRUSt2 v2.5+ `pathways_out/path_abun_unstrat_described.tsv` |

### `data_formicrobiome.xlsx`
| Property | Value |
|---|---|
| Format | Microsoft Excel (.xlsx) |
| Sheet | Sheet 1 (default) |
| Key columns | `PC_code`, `AgeM`, `AgeGroup`, `Sexchildren`, `KGcode` |
| Clinical columns | `Dietary_pattern`, `Nutritionalstatus`, `HEI`, `WAZ`, `HAZ`, `BAZ` |
| Infection columns | `infectioncoinfection`, `BH`, `GL`, `DF` |
| Notes | `BH` = *Blastocystis hominis*, `GL` = *Giardia lamblia*, `DF` = *Dientamoeba fragilis* (binary 0/1) |

---

## Preprocessed Intermediate Files

The Jupyter notebook (`scripts/funAnalysis_metagPastoKids.ipynb`) generates the
following intermediate files consumed by the R script. These are written to
`data/microbiome_analysis_output/` and are also NOT committed to git.

```
data/microbiome_analysis_output/
├── preprocessing_output/
│   ├── asv_table_clr_transformed.tsv
│   ├── KO_Functions_clr_transformed.tsv
│   ├── EC_Functions_clr_transformed.tsv
│   └── MetaCyc_Pathways_clr_transformed.tsv
└── analysis_ready_data/
    ├── asv_taxonomy_clean.csv
    └── metadata_clean.csv
```

All CLR-transformed matrices have:
- Rows = feature IDs (`ASV_ID` / `function` / `pathway`)
- Columns = sample IDs (`PC_code`)
- Values = Centered Log-Ratio transformed abundances (float64)

---

## Data Availability

The primary cohort data supporting this study are available upon reasonable
request from the corresponding author (Dr. rer. nat. Guillermo G. Torres,
guigotoe@gmail.com) subject to ethical approval requirements.

If the dataset has been deposited in a public repository, accession numbers
will be listed here upon publication.

| Archive | Accession | Notes |
|---|---|---|
| ENA / SRA | TBD | 16S rRNA amplicon sequences |
| Zenodo | TBD | Processed metadata (anonymized) |

---

## PICRUSt2 Functional Prediction

Functional annotations (KO, EC, MetaCyc) were predicted from 16S rRNA
amplicon data using **PICRUSt2 v2.5.2**:

```bash
picrust2_pipeline.py \
  -s asv_sequences.fna \
  -i ASV_table.biom \
  -o picrust2_output/ \
  --stratified \
  -p 4
```

Reference: Douglas et al. (2020) *Nature Biotechnology* 38, 685–688.
DOI: [10.1038/s41587-020-0548-6](https://doi.org/10.1038/s41587-020-0548-6)
