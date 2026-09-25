# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- `set_structure()`: a structure Configuration that asks for the knapsack
  formulation constructs the `BinaryKnapsackBlock` of the facilities when
  the `BlockConfig` is applied, before any abstract representation, so that
  whoever reads the tree of sub-Block sees them, e.g., a
  `LagrangianDualSolver` that decomposes the Block recursively; the
  knapsacks are loaded again when the abstract representation is generated
  if the data have changed since, and `load()` keeps the structure

- in the knapsack formulation the Block has an `FRealObjective` of its own,
  with no Variable, besides those of the knapsacks, so that it can be the
  Block of a `LagBFunction`, e.g., as the scenario of a
  `TwoStageStochasticBlock` decomposed by a `LagrangianDualSolver`

### Changed

- the data archive is downloaded by version: `DATA_VERSION` in CMakeLists.txt
  names the version of the Package Registry to read, and the archive and the
  marker of its extraction carry it in their name, so that a tree holding
  an older extraction (the cache of the CI, or a clone extracted before)
  downloads and extracts again instead of running on the old data;
  data/upload-txt publishes the archive under that version

- the makefile asks for `-O3 -DNDEBUG` and nothing else, the macro of the
  patch for `boost::any` on macOS having no reason to be there since there is
  no `boost::any` left in the core

- whoever links the module keeps it: the classes of a module register
  themselves in the factory from a static initialiser, and a linker that
  drops what looks unused takes the registration away with it, so the target
  now tells whoever links it to keep the symbol that forces the module in,
  and on ELF, where naming the symbol is not enough, the library as a whole

- `chg_facility_costs()`, `chg_transportation_costs()`,
  `chg_facility_capacities()` and `chg_customer_demands()` take their data
  as a `std::span< const double >`, whose length they check against the
  Range or the Subset instead of reading past the end, and are registered in
  the methods factory in that form too; the forms taking an iterator stay,
  and defer to the span ones. What they pass on to the BinaryKnapsackBlock
  and the MCFBlock inside is a span as well

### Fixed

- on macOS a program linking the module lost the classes the module
  registers in the factories when the linker dropped the library, as it
  does under `-dead_strip_dylibs`, which conda sets: the target now asks the
  linker for the symbol that forces the module in (`-u`), which ld64,
  unlike the ELF linker, counts as a use of the library

- the step that fetches the data archive of this module says what went wrong
  when it goes wrong: the download is checked, an archive that did not arrive
  is removed instead of being left on disk for the build to take for the real
  one, and the message names the URL. A server that answers with an error page
  used to leave a file of a few bytes there, which made the next build fail
  while extracting it, with the message of `tar` and no mention of the
  download

- a facility fixed open or closed in the data (`FacilityFix`) fixes, in
  the knapsack formulation, the opening of the facility, i.e., the last item
  of its knapsack, instead of the item with the index of the facility,
  i.e., the assignment of a customer

- `close_facilities()`, in both its forms, counted the facilities that were
  already fixed where it acts on the free ones, so that it returned without
  doing anything whenever all of them were free, the state every instance
  starts in

- the Subset form of `"CapacitatedFacilityLocationBlock::close_facilities"`
  was registered in the methods factory on `open_facilities()`, so that
  closing facilities by name opened them

- the Subset form of `chg_transportation_costs()` was not registered in the
  methods factory, while the documentation said it was

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

## [0.1.1] - 2024-02-27

### Changed

- makefiles and Cmake files updated to new global SMS++ scheme

- documentation updated accordingly

## [0.1.0] - 2022-06-28

### Added

- First test release.

[Unreleased]: https://gitlab.com/smspp/capacitatedfacilitylocationblock/-/compare/0.3.0...develop
[0.3.0]: https://gitlab.com/smspp/capacitatedfacilitylocationblock/-/compare/0.2.0...0.3.0
[0.2.0]: https://gitlab.com/smspp/capacitatedfacilitylocationblock/-/compare/0.1.1...0.2.0
[0.1.1]: https://gitlab.com/smspp/capacitatedfacilitylocationblock/-/compare/0.1.0...0.1.1
[0.1.0]: https://gitlab.com/smspp/capacitatedfacilitylocationblock/-/tags/0.1.0
