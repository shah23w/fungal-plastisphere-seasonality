# Season shapes fungal enrichment on riverine microplastics relative to co-collected natural particles

Code and derived data for the paper. The analysis is a reanalysis of publicly
archived shotgun metagenomes. No new sequencing was performed.

Each tagged release of this repository is archived on Zenodo. The DOI is
given in the Data and Code Availability statement of the paper.

Sequence data are archived under NCBI BioProjects PRJNA1159715 and PRJNA1330252.
Accessions for the wider corpus of 859 archived runs are listed in
`tables/TableS2_detection_envelope.tsv`.

## Directory structure

```
tables/            Tables S1-S25 and the external inputs they are built from
figures/           Figures 1-4 and S1-S12, as PNG and PDF, main figures also TIFF
scripts/analysis/  runs from tables/ alone, on a laptop
scripts/pipeline/  read processing and classification, needs the raw corpus
MANIFEST_md5.txt   checksum of every file
```

## Requirements

Python 3.9 or later with NumPy, SciPy and Matplotlib. The deposit was last run
with Python 3.12.7, NumPy 1.26.4, SciPy 1.13.1 and Matplotlib 3.10.0. The
published analysis used Python 3.13 and SciPy 1.16.0. Results are identical on
both. No other packages are needed and nothing is installed by the scripts.

The pipeline scripts also require Kraken2 2.1.3, fastp 0.23.4, seqkit,
and wgsim as distributed with SAMtools 1.19.

## scripts/analysis

Each script takes an optional first argument giving the table directory and, for
the figure scripts, an optional second argument giving the output directory.
Run from the repository root:

```
python3 scripts/analysis/01_check_main_values.py
```

| script | what it does |
|---|---|
| `01_check_main_values.py` | recomputes the headline values of the paper from the tables and compares each with the value in the text |
| `02_check_remaining_values.py` | recomputes the remaining values in the Results, so that between 01 and 02 every number is derived rather than asserted |
| `03_season_summary.py` | builds Table S13, with Hodges-Lehmann estimates and exact Walsh intervals, and the leave-one-site-out F range in Table S10 |
| `04_temperature_slopes.py` | every value in the Results subsection on water temperature, including the covariate screen and the multivariable fit |
| `05_antecedent_rainfall.py` | the antecedent-rainfall variables and every rainfall statistic in the Results and in Tables S9 and S10 |
| `06_genus_slopes.py` | per-genus temperature slopes with Holm adjustment across the 12 tests, writes Table S24 |
| `07_crossstudy_dates.py` | joins sampling months and season into Table S3 and checks the single-timepoint claim |
| `08_bacterial_contrast.py` | bounds the technical contribution to the winter result using the bacterial contrast on the same blocks |
| `09_residual_diagnostics.py` | residual normality and equality of variance, and refits on ranks and site-demeaned differences, appended to Table S10 |
| `10_figure1_figure4.py` | Figure 1 and Figure 4 |
| `11_figure2.py` | Figure 2 |
| `12_figure3.py` | Figure 3 |
| `13_supporting_figures.py` | Figures S1 to S12, and Tables S11 and S23 |

Scripts 01 to 04, 07 to 09, 12 and 13 recompute each value they report and
compare it with the published one. Every one of them prints a count of checks and
mismatches and exits non-zero if any check fails. Scripts 05, 06, 10 and 11 derive
and plot values that the other scripts check, and carry no checks of their own.

Scripts 03, 06, 07, 09 and 13 rewrite a table they build. The content does not
change, so re-running them is safe.

## scripts/pipeline

These produced the classification outputs the tables are built from. They need
the raw read corpus and a compute cluster and will not run from this deposit.
They are included so the processing can be inspected.

| script | what it does |
|---|---|
| `121_get_full_db.sh` | downloads the uncapped Kraken2 PlusPF database, build 2026-06-26 |
| `161_fungal_db_fetch.sh` | builds the fungal-only index from that build's own fungal library |
| `129_corpus_process.sh` | per-run trimming and classification at the settings used throughout |
| `123_build_mixtures.sh` | builds the fungus-free mock community |
| `175_mock2x2.sh` | classifies the mock libraries under the five configurations |
| `166_spikein_build.sh` | builds the spiked libraries for the recovery experiment |
| `168_spikein_classify.sh` | classifies them at the same settings |
| `134_reconcile_1159715.py` | reconstructs arm assignment from archive metadata |
| `204_genus_matrix_180.py` | builds Table S25 from the Kraken2 reports |

Classification throughout used the pre-built 16 GB capped PlusPF index of
26 June 2026 (`k2_pluspf_16_GB_20260626`), taken as distributed, at a confidence
threshold of 0.1 with a minimum of three hit groups. Both indexes come from the
Kraken2 index collection catalogued at https://benlangmead.github.io/aws-indexes/k2,
which serves them as `k2_pluspf_16_GB_20260626.tar.gz` and
`k2_pluspf_20260626.tar.gz`. The checksum manifest published with this build is at
https://genome-idx.s3.amazonaws.com/kraken/pluspf_20260626/pluspf.md5.

## tables

Tables S1 to S25 are the Supporting Information tables of the paper. The
remaining files are inputs they are built from.

| file | contents |
|---|---|
| `bao2026_covariates.tsv` | environmental variables published with the source study |
| `bao2026_water_temperature.tsv` | water temperature per site and season |
| `bao2026_phylum_abundance_FROM_SI.tsv` | genome-resolved phylum abundances and genome counts from the source study |
| `block_sampling_dates.tsv` | sampling date of each of the 60 blocks |
| `crossstudy_season.tsv` | sampling months and inferred season for the cross-study set |
| `HKO_daily_rainfall_2023.tsv` | daily rainfall averaged across the Hong Kong Observatory stations reporting each day |
| `nz_collection_dates.tsv` | collection dates for the two southern-hemisphere studies |
| `SI_TABLE_INDEX.tsv` | file name and row count of each table |

## Checksums

`MANIFEST_md5.txt` lists the MD5 of every file. To verify:

```
md5sum -c MANIFEST_md5.txt
```

On macOS, `md5 -q` replaces `md5sum`.

## Licence

Code is released under the MIT Licence. Derived data tables and figures are
released under CC BY 4.0.
