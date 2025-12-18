# Changelog

## [3.1.1](https://github.com/snakemake-workflows/rna-seq-star-deseq2/compare/v3.1.0...v3.1.1) (2025-12-18)


### Bug Fixes

* get sra download to work and test it ([#102](https://github.com/snakemake-workflows/rna-seq-star-deseq2/issues/102)) ([a9854f6](https://github.com/snakemake-workflows/rna-seq-star-deseq2/commit/a9854f62574e5879323bd200e6e4b4272a975ce9))

## [3.1.0](https://github.com/snakemake-workflows/rna-seq-star-deseq2/compare/v3.0.1...v3.1.0) (2025-09-10)


### Features

* add new (logging) output to star_align ([#97](https://github.com/snakemake-workflows/rna-seq-star-deseq2/issues/97)) ([d365a2c](https://github.com/snakemake-workflows/rna-seq-star-deseq2/commit/d365a2c1859f69d014054fbf0f32b7e8b4726df6))

## [3.0.1](https://github.com/snakemake-workflows/rna-seq-star-deseq2/compare/v3.0.0...v3.0.1) (2025-09-05)


### Bug Fixes

* consistently provide extra params for star rules ([#95](https://github.com/snakemake-workflows/rna-seq-star-deseq2/issues/95)) ([b975537](https://github.com/snakemake-workflows/rna-seq-star-deseq2/commit/b975537fe32e696691ec3d72e785fb7524b7d252))

## [3.0.0](https://github.com/snakemake-workflows/rna-seq-star-deseq2/compare/v2.1.2...v3.0.0) (2025-09-04)


### ⚠ BREAKING CHANGES

* use fastp instead of cutadapt for adapter trimming ([#93](https://github.com/snakemake-workflows/rna-seq-star-deseq2/issues/93))
* update all snakemake wrapper and conda environment tool versions

### Features

* use fastp instead of cutadapt for adapter trimming ([#93](https://github.com/snakemake-workflows/rna-seq-star-deseq2/issues/93)) ([3c86b1a](https://github.com/snakemake-workflows/rna-seq-star-deseq2/commit/3c86b1a3124bf30a084d489acd1d9f85e55f1060))

## [2.1.2](https://github.com/snakemake-workflows/rna-seq-star-deseq2/compare/v2.1.1...v2.1.2) (2024-06-05)


### Bug Fixes

* use derived input for star_index ([#81](https://github.com/snakemake-workflows/rna-seq-star-deseq2/issues/81)) ([87fffe6](https://github.com/snakemake-workflows/rna-seq-star-deseq2/commit/87fffe6a1beaa86e95c3564061d2720cc73308c7))

## [2.1.1](https://github.com/snakemake-workflows/rna-seq-star-deseq2/compare/v2.1.0...v2.1.1) (2024-03-25)


### Bug Fixes

* release-please branch to `master` and set permissions ([#79](https://github.com/snakemake-workflows/rna-seq-star-deseq2/issues/79)) ([4b781cf](https://github.com/snakemake-workflows/rna-seq-star-deseq2/commit/4b781cfa14fb5474108594fbaefa0ac8519f19dc))
* remove unused ftp RemoteProvider and require recent snakemake 8 ([#76](https://github.com/snakemake-workflows/rna-seq-star-deseq2/issues/76)) ([0f18be7](https://github.com/snakemake-workflows/rna-seq-star-deseq2/commit/0f18be7618a8dfb998455edf1da89b7cfb2e1301))


### Performance Improvements

* update all wrapper to latest v3.5.3 ([#78](https://github.com/snakemake-workflows/rna-seq-star-deseq2/issues/78)) ([bc9ab71](https://github.com/snakemake-workflows/rna-seq-star-deseq2/commit/bc9ab713f7c11b04bae296a27970aceeb12ab1ae))
