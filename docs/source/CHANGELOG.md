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
- 3D conformation in the `mol` attribute of the Individual during initialization.
- The use of `grow_mol` in the initialization of the populating when `get_similar = True`. Now the population is initialized with `mutate_mol` and the same set of crem parameters used during the searching.
- The automatic addition of Hs in the case where `min_size` and/or `max_size` were equal to zero. Now if your intention is work with the hydrogens, you must provided a SMILES with the explicit Hs. In the future the input could be just a RDKit mol.

### Added

- Handling vina RuntimeError and keeping track for debug. This feature is used to identify what is the error. In the future will be removed.
- Two new fitness functions: `moldrug.fitness.CostOnlyVina` and `moldrug.fitness.CostMultiReceptorsOnlyVina`. They only use the information of vina scoring function. See the [docs](https://moldrug.readthedocs.io/en/latest/source/modules/fitness.html) for more info about it.
- Tracking of atom indexes during generations in order to use `protected_ids` and `replace_ids` options of `mutate_mol` function of [CReM](https://github.com/DrrDom/crem). Before it was not possible; the use of these features generate undesired solutions because the indexes are not static over generations. Even so, there are still some problems for symmetric molecules. We are working on it.

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
