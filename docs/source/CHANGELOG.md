# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [unreleased]

### Fixed

- Hidden `RuntimeError` in `moldrug.fitness` module.
- Bug printing number of generations in final info.

### Removed

- Unused code in `moldrug.home` module.

### Added

- Handling vina RuntimeError and keeping track for debug. This feature is used to identify what is the error. In the future will be removed.
- Two new fitness functions: `moldrug.fitness.CostOnlyVina` and `moldrug.fitness.CostMultiReceptorsOnlyVina`. They only use the information of vina scoring function. See the [docs](https://moldrug.readthedocs.io/en/latest/source/modules/fitness.html) for more info about it.

## [0.1.0] - 2022-07-25

### Fixed

- Minor code cleaning.
- Bug during the import of user custom cost function.

### Added

- `outdir` option for the command line.
- User custom desirability.

## [0.0.4] - 2022-07-21

### Fixed

- Minor compatibility issue with Python 3.8 (issue [#4](https://github.com/ale94mleon/MolDrug/issues/4)).
- Problem with the user custom cost function supplied on the command line.
- `Local` class compatible with the command line.
- Minor code cleaning.
- Better code covered during testing

[unreleased]: https://github.com/ale94mleon/MolDrug/compare/0.1.0...HEAD
[0.1.0]: https://github.com/ale94mleon/MolDrug/compare/0.0.4...0.1.0
[0.0.4]: https://github.com/ale94mleon/MolDrug/compare/0.0.3...0.0.4
