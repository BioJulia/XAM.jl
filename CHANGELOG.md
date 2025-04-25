# Changelog

All notable changes to XAM.jl will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]
- Added conversion capability for SAM to BAM ([#63](https://github.com/BioJulia/XAM.jl/pull/63)).

## [0.4.1]

### Fixed
- Corrected `isproperpair` method ([#78](https://github.com/BioJulia/XAM.jl/pull/78)).

## [0.4.0]

### Added
- Added BAM.Reader index support for BAI object ([#56](https://github.com/BioJulia/XAM.jl/pull/56)).
- Added doi badge.
- Added test to ensure EOF_BLOCK gets written.
- Added `isreversestrand`.
- Added `isfirstsegment`.
- Added `islastsegment`.

### Changed

- Subtype from XAMReader and XAMWriter from common abstract types.
- Subtype from XAMRecord.
- Unified flag queries.
- Improved Slack link.
- Updated to use [Automa](https://github.com/BioJulia/Automa.jl) v1 ([#65](https://github.com/BioJulia/XAM.jl/pull/65)).
- Pointed the Unit Tests badge at the develop branch.
- Pluralised flag.
- Renamed `ismateunmapped` to `isnextunmapped`.
- Renamed `isreverse` to `isreversecomplemented`.
- Renamed `isforward` to `isforwardstrand`.
- `ispositivestrand` aliases `isforwardstrand`.
- `isnegativestrand` aliases `isreversestrand`.
- Renamed `ismatereverse` to `isnextreversecomplemented`.
- `isread1` aliases `isfirstsegment`.
- `isread2` aliases `islastsegment`.

### Fixed
- Updated hts-files.md ([#62](https://github.com/BioJulia/XAM.jl/pull/62)).
- Corrected the behaviour of `isprimaryalignment` with `isprimary`.

### Removed
- Moved the functionality of `isprimary` into `isprimaryalignment`.


## [0.3.1]

### Changed

- Upgraded to BioAlignments v3 ([#55](https://github.com/BioJulia/XAM.jl/pull/55)).


## [0.3.0] - 2022-10-10

## Added

- Crosschecks for SAM and BAM ([#29](https://github.com/BioJulia/XAM.jl/pull/29)).
- Improved documentation for flags ([#43](https://github.com/BioJulia/XAM.jl/pull/43)).

### Changed

- `BAM.quality` performance improved ([#21](https://github.com/BioJulia/XAM.jl/issues/21)).
- Updated BioAlignments to v2.2 and BioSequences to v3 ([#48](https://github.com/BioJulia/XAM.jl/pull/48)).

### Fixed

- `BAM.Record` layout now matches the BAM specs ([#26](https://github.com/BioJulia/XAM.jl/pull/26)).

[Unreleased]: https://github.com/BioJulia/XAM.jl/compare/v0.4.0...HEAD
[0.4.0]: https://github.com/BioJulia/XAM.jl/compare/v0.3.1...0.4.0
[0.3.1]: https://github.com/BioJulia/XAM.jl/compare/v0.3.0...v0.3.1
[0.3.0]: https://github.com/BioJulia/XAM.jl/compare/v0.2.8...v0.3.0
