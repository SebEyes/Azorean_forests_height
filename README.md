# Azorean forest vegetation height

A Julia workflow for deriving **1 m vegetation-height rasters** from the Portuguese Directorate-General for Territory (DGT) ACORES-LAZ classified LiDAR point-cloud collection. The pipeline retrieves LAZ tiles from the DGT catalogue, estimates terrain and vegetation surfaces, writes GeoTIFF canopy-height products, mosaics them by island, and provides raster quality-control summaries.

For each 1 m cell, the processor estimates a digital terrain model (DTM) from class 2 ground returns, fills nearby ground gaps, takes the maximum elevation of vegetation returns in classes 3, 4, and 5 (DSMv), and calculates `max(0, DSMv - DTM)`. Cells with no vegetation are zero; cells without ground support inside the fill radius are NoData. Other point classes are excluded. The source raster CRS is read from the LAZ header where available; use `--epsg` when it is absent.

> **Code acknowledgment:** The code in this repository was written by **OPENCODE (Model Big pickle)**. The author carried out validation tests to confirm that it was operating correctly.

## Repository layout

```text
.
├── Project.toml                    # Direct Julia project dependencies
├── Manifest.toml                   # Resolved Julia environment (Julia 1.13.3)
├── code/
│   ├── 00-julia_env_initialisation.jl
│   ├── download_azores_laz.jl      # DGT catalogue search and LAZ download
│   ├── laz_to_vegetation_height.jl # One LAZ tile -> vegetation-height GeoTIFF
│   ├── compute_archipelago.jl      # Download, process, and mosaic by island
│   ├── raster_quality_check.jl     # Distribution/outlier QC reports
│   └── filter_raster_p99.jl        # Mask values above each island mosaic's estimated P99
├── data/vegetation_height/qc/      # Committed QC outputs
└── doc/Data_source.pdf             # Data-source documentation
```

LAZ point clouds and GeoTIFF rasters (`*.laz`, `*.tif`) are excluded by `.gitignore` due to their size; obtain or generate them as described below. The tracked `data/vegetation_height/qc/` folder contains QC tables, an HTML report, and related diagnostic outputs.

## Requirements

- Julia **1.13.3**.
- Julia packages: `ArchGDAL`, `HTTP`, `JSON3`, and `PointClouds` (plus the Julia standard library `Dates`). `Manifest.toml` locks the resolved package versions.
- For DGT downloads, access to the DGT ACORES-LAZ catalogue and a valid DGT login (including an OTP if enabled), or an authenticated browser `connect.sid` cookie.
- Storage for large point-cloud inputs and raster outputs. The pipeline removes temporary LAZ files by default, but generated rasters can be large.

From the repository root, instantiate the project environment:

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

The script `code/00-julia_env_initialisation.jl` also activates and instantiates the project.

## Usage

Run commands from the repository root. Each script provides `--help` for its full option list.

### 1. Inspect the DGT catalogue

No login is needed to list catalogue items:

```sh
julia code/download_azores_laz.jl --list-only
julia code/download_azores_laz.jl --island COR --limit 2 --dry-run
```

### 2. Download LAZ tiles (optional standalone step)

Provide credentials as command-line options or through environment variables. Keep secrets out of shell history and never commit them.

```sh
export DGT_USERNAME='your-account'
export DGT_PASSWORD='your-password'
# If required by the account:
export DGT_OTP='current-one-time-code'

julia -t 4 code/download_azores_laz.jl \
  --island COR --limit 2 --outdir data/azores_laz
```

Alternatively, the downloader accepts `--cookie VALUE` for a browser session cookie. Do not share or commit a session cookie. To download the full catalogue, omit `--island` and `--limit`. The catalogue downloader's default output directory is `data/azores_laz`.

### 3. Process one LAZ tile

```sh
julia code/laz_to_vegetation_height.jl INPUT.laz \
  -o data/vegetation_height/tiles/OUTPUT_vegheight.tif
```

The processor defaults to 1 m resolution, ground class 2, vegetation classes 3/4/5, and chunked point processing. See `--help` for options including `--fill-distance`, `--chunk-size`, `--threads`, `--nodata`, and `--epsg`.

### 4. Run the end-to-end island/archipelago pipeline

The pipeline searches the catalogue, authenticates, downloads and processes tiles, removes each downloaded LAZ after processing, and merges available tile rasters into one GeoTIFF per island. It may take substantial time, network traffic, CPU, and disk space.

```sh
export DGT_USERNAME='your-account'
export DGT_PASSWORD='your-password'
# If required:
export DGT_OTP='current-one-time-code'

julia --threads=4 code/compute_archipelago.jl \
  --island COR --limit 2 \
  --tiles-dir data/vegetation_height/tiles \
  --outdir data/vegetation_height
```

Omit `--island COR --limit 2` to target the whole archipelago. Useful controls include `--resume` (reuse existing tile rasters), `--delete-tiles` (remove per-tile rasters after mosaic creation), `--keep-laz` (retain downloaded LAZ files for debugging), and `--workdir DIR`. By default, per-tile products go to `data/vegetation_height/tiles`, island mosaics to `data/vegetation_height`, and temporary LAZ files to a system temporary directory. An island mosaic is named `<CODE>_vegetation_height.tif` (e.g. `COR_vegetation_height.tif`).

### 5. Run raster quality checks

After generating island mosaics:

```sh
julia --threads=6 code/raster_quality_check.jl \
  --raster-dir data/vegetation_height \
  --outdir data/vegetation_height/qc
```

The QC script scans raster windows and writes `height_distribution_by_island.csv` (island summaries), `height_histogram.csv`, `height_outlier_blocks.csv`, `height_extreme_pixels.csv`, and the self-contained `height_distribution_report.html`. Its defaults and optional island/threshold settings are documented by `--help`.

### 6. Filter island mosaics above their estimated P99 (optional)

To create filtered copies after mosaicking (and optionally after QC), run:

```sh
julia --threads=4 code/filter_raster_p99.jl \
  --raster-dir data/vegetation_height \
  --outdir data/vegetation_height
```

By default, the script finds merged inputs named `<CODE>_vegetation_height.tif`, computes a separate P99 threshold for each island from all valid pixels, and writes `<CODE>_vegetation_height_p99.tif`. Values strictly above the threshold are changed to NoData; the original mosaics are left unchanged. `p99_filter_summary.csv` records the threshold, valid-pixel count, NoData/non-finite count, masked-pixel count, and percentage masked for each island. The threshold is a histogram-based estimate (default bin width `0.01` m); values at or above the default `--max-height 120` m share an overflow bin, so it is not an exact raw-pixel quantile. Raster scanning/writing is windowed (default 256 rows) to limit memory use.

Use `--island COR` to filter only one island, `--bin-width M` or `--max-height M` to change the histogram settings, `--block-rows N` to change the window height, and `--suffix STR` to change the output suffix. Defaults are `--raster-dir data/vegetation_height` and `--outdir` equal to the input directory. For example, to keep filtered outputs separate:

```sh
julia code/filter_raster_p99.jl \
  --raster-dir data/vegetation_height \
  --outdir data/vegetation_height/qc \
  --island COR
```

Run `julia code/filter_raster_p99.jl --help` for all options. The script activates and instantiates the repository Julia environment at startup.

## Data source and interpretation

The primary point-cloud source for the Julia pipeline is the Portuguese Directorate-General for Territory (DGT) ACORES-LAZ collection, accessed through the DGT catalogue (`https://cdd.dgterritorio.gov.pt/dgt-fe/mapa?collection=ACORES-LAZ`); the downloader enumerates this collection and requires DGT authentication to download its LAZ assets. Consult [`doc`](doc/) for official source information and the scripts' help/comments for algorithm details. **The generated vegetation-height surface is an estimate based on classified returns and the configured ground-gap filling, not a direct field measurement**. P99-filtered rasters are further processed copies that mask high values; they do not replace the source mosaics. Review the QC products and validate against the intended study use before interpreting extreme values.

## Code maintenance

For code maintenance, contact [Sébastien Lhoumeau](https://github.com/SebEyes).