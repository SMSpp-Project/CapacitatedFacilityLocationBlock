# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

### Changed

### Fixed

## [0.3.0] - 2026-09-12

### Added

- the two "Benders friendly" formulations of the problem, one with slack arcs
  and one with feasibility cuts

### Changed

- the converter of the text instances to netCDF is `cfl2nc4`, rather than
  `txt2nc4`

- the Solution versions of `is_feasible()` and `is_optimal()` are
  `is_sol_feasible()` and `is_sol_optimal()`, as in the core

- the version of the module is the git tag of its repository, or the
  VERSION.txt of a release tarball, and the shared library carries it: its
  SONAME is major.minor while the major is 0, and it is installed with an
  RPATH relative to itself, so that an installed tree keeps working wherever
  it is moved

### Removed

- ScenarioReductionSolver and its tester, superseded by the standalone
  ScenarioReductionSolver module

## [0.2.0] - 2025-12-12

### Added

- tester for ScenarioReductionSolver

- [huge] ScenarioReductionSolver implementing a bunch of
  fast heuristics for the specific version of the problem
  used in scenario reduction applications

- possibility of a constraint on maximum number of facilities
  (for scenario reduction applications)

### Changed 

- major data handling upgrade: instances are not included
  in the repo for space/time efficiency but a script to
  download them is provided

- adapted to new standard organization of makefiles

### Fixed 


## [0.1.1] - 2024-02-27

### Changed

- makefiles and Cmake files updated to new global SMS++ scheme

- documentation updated accordingly

## [0.1.0] - 2022-06-28

- First test release.


[Unreleased]: https://gitlab.com/smspp/capacitatedfacilitylocationblock/-/compare/0.3.0...develop
[0.3.0]: https://gitlab.com/smspp/capacitatedfacilitylocationblock/-/compare/0.2.0...0.3.0
[0.2.0]: https://gitlab.com/smspp/capacitatedfacilitylocationblock/-/compare/0.1.1...0.2.0
[0.1.1]: https://gitlab.com/smspp/capacitatedfacilitylocationblock/-/compare/0.1.0...0.1.1
[0.1.0]: https://gitlab.com/smspp/capacitatedfacilitylocationblock/-/tags/0.1.0
