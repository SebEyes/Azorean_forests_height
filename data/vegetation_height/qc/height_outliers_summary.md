# Vegetation-height distribution and outliers, by island

Source: the **merged per-island mosaics** (`data/vegetation_height/*_vegetation_height.tif`),
not the individual tiles. Analysis script: `code/raster_quality_check.jl`.

## Method

- Every valid (non-NoData) 1 m pixel was accumulated into a 0.25 m histogram,
  scanned in row windows so the largest mosaic is never loaded whole.
- Reported statistics: percentiles, mean/std, Tukey fences (Q3 ± 1.5×IQR and
  Q3 ± 3×IQR), and exact pixel counts above 20/25/30/35/40/50 m.
- "Forest" statistics are restricted to pixels > 2 m (the high-vegetation class
  in the source LiDAR).
- No negative heights were found in any island (the CHM clips at 0).

Because 15–49 % of the pixels are exactly zero (no vegetation), Tukey fences are
strongly affected by the zero inflation: they range from 11 m (Corvo) to 54 m
(São Miguel) and are best treated as descriptive. The fixed 30/40/50 m
thresholds are more useful for flagging physically implausible values.

## Results

| Island | Valid pixels | % zero | Median (m) | P99 (m) | P99.9 (m) | Max (m) | Tukey 1.5 IQR (m) | > 30 m | > 40 m | > 50 m |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Corvo (COR) | 17,269,048 | 24.8 | 0.40 | 11.47 | 20.15 | 107.46 | 11.48 | 2,655 | 843 | 374 |
| Faial (FAI) | 174,049,600 | 32.7 | 0.25 | 23.33 | 30.41 | 93.93 | 38.67 | 203,877 | 1,588 | 113 |
| Flores (FLO) | 141,709,372 | 15.5 | 0.66 | 21.10 | 29.53 | 130.38 | 35.70 | 125,635 | 22,372 | 10,311 |
| Graciosa (GRA) | 61,676,711 | 42.1 | 0.22 | 20.51 | 27.19 | 75.14 | 30.75 | 14,042 | 956 | 369 |
| Pico (PIC) | 446,482,356 | 23.9 | 0.58 | 19.99 | 26.20 | 85.55 | 35.76 | 114,722 | 8,574 | 1,496 |
| São Jorge (SJO) | 245,102,934 | 21.8 | 0.46 | 19.02 | 29.31 | 158.04 | 30.46 | 208,326 | 20,036 | 5,438 |
| Santa Maria (SMA) | 98,300,335 | 49.0 | 0.22 | 25.68 | 34.61 | 153.76 | 39.62 | 358,593 | 21,738 | 4,614 |
| São Miguel (SMG) | 739,124,608 | 31.1 | 0.57 | 29.00 | 37.20 | **1066.49** | 53.57 | **5,654,189** | **360,692** | **44,253** |
| Terceira (TER) | 402,254,525 | 38.5 | 0.23 | 26.50 | 34.31 | 136.15 | 42.48 | 1,570,389 | 50,993 | 7,125 |

## Outlier assessment

- All islands have a long right tail: P99 is 19–29 m, but maxima reach 75–160 m,
  and São Miguel reaches 1066 m. Values above ~40 m are almost certainly
  artefacts (misclassified noise returns, DTM interpolation gaps, or merge
  artefacts), not trees.
- **São Miguel is the clear priority**: 360,692 pixels > 40 m (0.049 % of its
  valid pixels) and 44,253 > 50 m, including a compact 5-pixel spike of
  1060–1066 m near EPSG:5015 x ≈ 619,706, y ≈ 4,180,002
  (col 19083–19084, row 16582–16584), plus a cluster near x ≈ 644,800,
  y ≈ 4,180,130 reaching ~247 m.
- Other islands with notable high spikes: São Jorge (158 m), Santa Maria
  (154 m), Flores (130 m), Terceira (136 m), Corvo (107 m).
- The outliers are spatially clustered rather than uniformly scattered; see
  `height_outlier_blocks.csv` (256-row windows) and
  `height_extreme_pixels.csv` (all 487,792 pixels above 40 m, with UTM
  coordinates).

## Files

- `height_distribution_by_island.csv` — one summary row per island.
- `height_histogram.csv` — 0.25 m bins, one row per island/bin.
- `height_outlier_blocks.csv` — outlier counts per 256-row window.
- `height_extreme_pixels.csv` — pixel locations above 40 m.
- `height_distribution_report.html` — histograms and top outlier windows per island.
