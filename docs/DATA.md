# Data provenance

Raw data is not stored in this repository. `scripts/00_download_data.py` retrieves it
and writes `data/raw/PROVENANCE.json`, recording the URL, retrieval date, file size and
SHA-256 checksum of every file. If an upstream source is revised, the checksum changes
and the difference is visible rather than silently altering results.

## Case counts — OpenDengue v1.3

Clarke, J. et al. *A global dataset of publicly available dengue case count data.*
Scientific Data 11, 296 (2024). https://doi.org/10.1038/s41597-024-03120-7
Licence: CC-BY 4.0. Downloaded from the project's GitHub release.

The "temporal extract" is used, which retains the finest temporal resolution available
for each location. The Pakistan subset contains 922 records spanning 1994–2025, of
which 795 are weekly.

### What is usable, and what is not

Only unbroken weekly runs can be modelled: a gap shifts the alignment between cases and
climate covariates for every week that follows. Ranking all Pakistani series by their
longest unbroken run gives:

| Series | Level | Longest run | Peak week | Usable |
|---|---|---|---|---|
| Pakistan | national | **77 weeks** (2012-12-23 → 2014-06-08) | 2,184 | yes |
| Sindh | province | 19 weeks (2021-09-19 → 2022-01-23) | 284 | yes |
| Khyber Pakhtunkhwa | province | 15 weeks (2021-08-29 → 2021-12-05) | 304 | yes |
| Tharparkar | district | 19 weeks | 148 | marginal |
| Haripur | district | 15 weeks | 204 | marginal |
| Hyderabad | district | 11 weeks | 154 | no — too short |
| **Punjab** | province | **6 weeks, 6 cases in total** | 6 | **no** |

**Punjab has no usable series.** This is not an artefact of the extraction: the province
contributes 28 weekly records totalling six cases, despite Punjab reporting the largest
share of Pakistan's confirmed dengue cases in several seasons. Incomplete subnational
coverage for South Asia is acknowledged as a known limitation by the OpenDengue authors,
who list greater disaggregation for India, Bangladesh, Nepal and Pakistan as a priority
for future data collection. Any Punjab-specific analysis would require obtaining data
directly from the Punjab health authorities or digitising published figures.

### Known limitations of the case data

- Counts are **reported** cases, not infections. Under-ascertainment for dengue is
  substantial and varies with healthcare access, awareness and testing capacity, so a
  reporting fraction is estimated rather than assumed to be one.
- Case definitions are not uniform across the record; the `case_definition_standardised`
  field distinguishes suspected from confirmed cases.
- The national series aggregates regional epidemics with different timing. This is the
  motivation for the aggregation-bias comparison against the 2021 provincial fits, not
  something the analysis can ignore.

## Climate — NASA POWER

NASA Prediction Of Worldwide Energy Resources, daily point data from the MERRA-2
reanalysis. Public domain; no key or registration required.

Parameters retrieved: `T2M` (mean temperature at 2 m), `T2M_MAX`, `T2M_MIN`,
`PRECTOTCORR` (corrected precipitation), `RH2M` (relative humidity at 2 m).

One series is retrieved per study window, at a representative location:

| Window | Location | Coordinates |
|---|---|---|
| national_2013 | Lahore | 31.55 N, 74.35 E |
| sindh_2021 | Karachi | 24.86 N, 67.01 E |
| kp_2021 | Peshawar | 34.01 N, 71.58 E |

Each request begins one year before the modelled window so that the lagged, smoothed
rainfall covariate is available from the first modelled week without back-filling.

Representing a national climate by a single city is a real simplification. Sensitivity
to this choice is tested by refitting with alternative locations; the result is reported
whether or not it is favourable.

## Processed outputs

`scripts/01_build_dataset.py` writes one tidy CSV per window to `data/processed/`, with
columns: `week_start`, `week_index`, `days_from_start`, `cases`, `T2M`, `T2M_MAX`,
`T2M_MIN`, `RH2M`, `PRECTOTCORR`, `rain_lagged`, `population`.

The script refuses to write a window whose weekly grid is broken or whose climate
coverage is too short for the configured lag, so a silently misaligned series cannot
reach the analysis.

| Window | Weeks | Total cases | Peak week | Mean temperature |
|---|---|---|---|---|
| national_2013 | 77 | 17,894 | 2,184 | 24.4 °C |
| sindh_2021 | 19 | 2,123 | 284 | 23.4 °C |
| kp_2021 | 15 | 1,899 | 304 | 21.4 °C |
