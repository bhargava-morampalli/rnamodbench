# nf-core/rnamodifications: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0-dev] - Unreleased

### Added

- Initial pipeline release
- Input validation with nf-schema
- Read mapping to 16S and 23S rRNA references using minimap2
- QC statistics using samtools flagstat, depth, and NanoPlot
- Signal-level data preparation:
  - Read ID extraction from mapped FASTQs
  - FAST5 subsetting using ont_fast5_api
  - Multi-read to single-read FAST5 conversion
  - Tombo resquiggling
  - F5C indexing
- Signal processing with F5C eventalign
- RNA modification detection using three tools:
  - Tombo level_sample_compare
  - Yanocomp (GMM-based detection)
  - Xpore (differential modification analysis)
- Support for native vs IVT comparison with 3 replicates each
- Comprehensive test suite with nf-test (23 module tests)
- Container support via Docker, Singularity, and Wave
- nf-core 3.5.x standards compliance

### Changed

- N/A

### Fixed

- N/A

### Deprecated

- N/A

### Removed

- N/A
