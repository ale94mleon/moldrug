# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [unreleased]

### Fixed
-   Hidden `RuntimeError` in fitness module

## [0.1.0] - 2022-07-25

### Fixed

- Minor cleaning code.
- Bug during the import of user custom cost function.

### Added

- `outdir` option for the command line.
- User custom desirability.

## [0.0.4] - 2022-07-21

### Fixed

- Minor compatibility issue with Python 3.8 (issue [#4](https://github.com/ale94mleon/MolDrug/issues/4)).
- Problem with the user custom cost function supplied on the command line.
- `Local` class compatible with the command line.
- Minor cleaning code.
- Better code covered during testing

[unreleased]: https://github.com/ale94mleon/MolDrug/compare/0.0.1...HEAD
[0.1.0]: https://github.com/ale94mleon/MolDrug/compare/0.0.4...0.0.1
[0.0.4]: https://github.com/ale94mleon/MolDrug/compare/0.0.3...0.0.4
