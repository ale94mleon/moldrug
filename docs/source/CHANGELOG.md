# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [unreleased]

### Fixed

- Refactored changes for for meeko-0.5.0.
- Update docs.

## [3.4.0] - 2023.03.10

### Added

- Continuation option to the command line.
- `moldrug.cli.CommandLineHelper` class to work with the parameters passed through the command line.
- `checkpoint` option to `moldrug.utils.GA`.
- [MolDrug-Dashboard](https://ale94mleon-moldrug-streamlitstreamlit-app-nltunu.streamlit.app) add-on. This is not included on the package itself, but could be used online or locally. In the case of locally you must check [Streamlit](https://streamlit.io/), the [requirements.txt](https://github.com/ale94mleon/MolDrug/blob/main/streamlit/requirements.txt) and the [app script](https://github.com/ale94mleon/MolDrug/blob/main/streamlit/streamlit_app.py).
- `retunr_mol` option to `utils.to_dataframe`

### Changed

- The warnings are not printed for: `moldrug.fitness._vinadock` and `moldrug.constraintconf.generate_conformers`. Now, only at the end of a MolDrug simulation a note will be printed if the `error.tar.gz` file is created.
- `moldrug.utils.run` does not print extra info if the command fails. It only raises the corresponded `RuntimeError`.
- `moldrug.fitness.__vinadock` by `moldrug.fitness._vinadock`.
- Remove conformers that clash with the protein in case of score_only, for local_only vina will handle the clash.

### Fixed

- Small bug during initialization of population with multiple `seed_mol`. Now `seed_mol` with the same amount of elements as `popsize` are not submitted to mutations, only evaluation.
- Problems with default parameters definition on the command line. Parameters with default values by the `type` of run defined in the configuration file are not needed to redefined any more; `moldrug.cli` will guess those.

## [3.3.0] - 2022.12.21

### Changed

- `moldrug.utils.__tar_errors__` by `moldrug.utils.tar_errors`.
- Default value of `moldrug.utils.tar_errors` is `error` instead of `.error`.
- `moldrug.constraintconf.generate_conformers` outputs warnings and errors to `error` instead of `.error`.
- `moldrug.fitness._vinadock` outputs warnings and errors to `error` instead of `.error`.

### Fixed

- The use of `moldrug.utils.tar_errors` inside of `moldrug.utils.Local` and `moldrug.utils.GA`.
- Clean code.

## [3.2.5] - 2022.12.20

### Fixed

- Improve docs.
- Cleaning of temporal files in `/tmp` directory. Now it is created temporal directories in the working directory with the pattern: `.costfunc_MolDrug_XXXXXXXX`.
- Cleaning errors. Now all warnings and errors are saved in .error directory and at the end they are compress to `error.tar.gz`.

## [3.2.2] - 2022.12.12

### Fixed

- Bug: The initial individual was not printed properly.

### Removed

- Redundant code in `moldrug.utils.GA`

## [3.2.0] - 2022.11.28

### Fixed

- Bug: The output error name when constraint fails has a idx prefix. E.g. `33_conf_27_error.pbz2` now is: `idx_33_conf_27_error.pbz2`. Now it is easy to delete all of this files at the end of the simulation if they are not needed. (on the last version the naming was not changing)

### Removed

- `moldrug.fitness.is_inside_box`
- Bug: `constraint_check_inside_box` option for the cost functions of `moldrug.fitness`

## [3.1.0] - 2022.11.28

### Added

- `moldrug.constraintconf.gen_aligned_conf`
- In case that `moldrug.constraintconf.generate_conformers` fails with `rdkit.Chem.AllChem.ConstrainedEmbed` it will try with `moldrug.constraintconf.gen_aligned_conf`.
- `moldrug.fitness.is_inside_box`
- `constraint_check_inside_box` arguments to the cost functions of `moldrug.fitness`. If the coordinates of the constraint conformation are outside the box; use always local_only, by default False

### Changed

- The output error name when constraint fails has a idx prefix. E.g. 33_conf_27_error.pbz2 now is: idx_33_conf_27_error.pbz2. Now it is easy to delete all of this files at the end of the simulation if they are not needed.

### Fixed

- Clean code.
- Improve docs.

## [3.0.3] - 2022.11.26

### Added

- Warning in case `moldrug.utils.Local` or `moldrug.utils.GA` are called with a different **MolDrug** as they were initialized.

### Changed

- Convert to absolute path `receptor_pdbqt_path` and `vina_executable` (in case that it points to a file) inside of `moldruf.fitness.__vinadock`.

### Fixed

- Bug for hydrogens coordinates when constrain docking was used.
- Improve docs.

## [3.0.1] - 2022.11.24

### Fixed

- Cleaning code.
- Sort the initial population based on the cost attribute when it is saved on disk.
- Improve docs.

## [3.0.0] - 2022.11.23

### Changed

- Name of `moldrug.fitness.get_mol_cost` to `moldrug.fitness.__get_mol_cost` function.
- The class `moldrug.utils.GA` does not have any more the method `roulette_wheel_selection`; now is part a function that could be called from `moldrug.utils`
- `max` for `max_conf` in `moldrug.constraintconf.constraintconf()` function.
- Entrance point constraintconf was changed to  constraintconf_moldrug and now it is link to `moldrug.cli.__constraintconf_cmd` instead `moldrug.constrainconf.constraintconf_cmd`.
- Name of the function `moldrug.fitness.vinadock` now is `moldrug.fitness.__vinadock`.
- Name of the function `moldrug.cli.moldrug_cmd` now is `moldrug.cli.__moldrug_cmd`.

### Fixed

- Cleaning the code.
- If `vina_executable` is provided (to any cost function) and it represents a path. It will be try to convert to absolute path. Previously relative path to the executable were not understood properly.
- Improve docs.

### Added

- `ad4map` in all the cost functions of the `moldrug.fitness` module. This parameters specify the path where the ad4 map files are. To use this feature you must have the AutoDcok Vina v1.2.3 of above. Now you can use the force fields of AD4 inside of Vina. Future release will extend the integration with this versions.
- `moldrug.utils.to_dataframe`. THis function was previously isolated as a method of the class `moldrug.utils.GA`; now it could also be called as a function.
- `kept_gens` attribute to the `Individual`s inside of `moldrug.utils.GA`. This is a set that contains the generations for which the Individual was conserved.
- `acceptance` attribute to `moldrug.utils.GA`. This is a dictionary that has as keyword the generation ID, and as values a dictionary with keywords: `accepted` (number of generated molecules accepted on the current generation) and `generated` (number of total molecules generated)
- Print `Accepted rate= accepted /  generated` during running.
- Add hydrogens before create pdbt file with meeko when constrain docking i used.
- `seed_mol` of `moldrug.utils.GA` now could be a list (or iterable in a general way) of RDKit molecules. This feature could be used to combine several MolDrug runs and create a final runs with this combined population.
- `seed_mol` from the command line could be: a valid SMILES, a list of valid SMILES or a list of path to the `_pop.pbz2` binary files. In the last case all the populations will be combined and sorted based on the cost attribute. If the result population is less that `popsize` new structures will be generated to complete the initial population. The individuals of this initial population will be reinitialized and the cost function will be calculated.

## [2.1.12] - 2022.09.29

### Added

- `score_only` bool parameter to `moldrug.fitness.get_mol_cost`.
- Print starting date when MolDrug is called from the command line.

### Removed

- Type Hints `int` for attribute `idx` on`moldrug.utils.Individual`.

### Changed

- If vina fails inside `moldrug.fitness.vinadock`; give as pdbqt the string "VinaFailed".
- If the molecule has a molecular weight highest than `wt_cutoff` when `moldrug.fitness.CostOnlyVina` (or `moldrug.fitness.CostMultiReceptorsOnlyVina`) is called; the pdbqt attribute of the returned Individual will be the string "TooHeavy" (or the list of strings List["TooHeavy"])

## [2.1.7] - 2022.09.02

### Fixed

- Bug on `moldrug.fitness.vinadock` during searching of MCS between `Individual.mol` and `constraint_ref`. Before was needed to manually specified the atom ids of the `seed_mol` that match `constraint_ref`, now it is not needed any more.
- Bug during handling exception in `moldrug.constraintconf.generate_conformers`.
- Bug during handling exception `moldrug.constraintconf.generate_conformers` in `moldrug.fitness.vinadock`.
- Bug(s) when constraint docking is used on different vina versions. The output of vina is not the same and therefore `moldrug.fitness.vinadock` failed.

### Changed

- In case `constraint = True` in `moldrug.fitness.vinadock`, `ref_smi` will be the the MCF between `individual.mol` amd `constraint_ref` instead the SMILES string of `constraint_ref` when `moldrug.constrainconf.generate_conformers` is internally called.

### Added

- `moldrug.fitness.get_mol_cost` function.
- Attribute `genID` to the generated individuals during a `moldrug.utils.GA` run.

## [2.1.0] - 2022.08.30

### Fixed

- Bug during the calculation of probabilities when costs are larger numbers.
- Expose hidden error if some Exception ocurred during parallel run.

### Added

- `moldrug.constrainconf` module.
- Raise `ValueError` if `ref_smi` is invalid in `moldrug.utils.constrainconf.generate_conformers`.

### Changed

- In case `constraint = True` in `moldrug.fitness.vinadock`, `ref_smi` will be the SMILES string of `constraint_ref` when `moldrug.constrainconf.generate_conformers` is internally called. This is in order to avoid error when `moldrug.utils.constrainconf.generate_conformers` tries to guess `ref_smi` based on MCS and fails, [see this RDKit bug](https://github.com/rdkit/rdkit/issues/5518). The work around for constraint docking is explained here: [Constraint Docking](https://moldrug.readthedocs.io/en/latest/notebooks/constraint_docking.html).
- `moldrug.fitness.generate_conformers` does not fail. In case of Exception it returns the same `mol` without conformers and write the error in a log file into the working directory.
- The attribute name `bestcost` by `best_cost` of `moldrug.utils.GA`.
- The functions `duplicate_conformers`, `get_mcs`, `generate_conformers`, `constraintconf` and `constraintconf_cmd` and the class `ProteinLigandClashFilter` were moved from `moldrug.fitness` module to `moldrug.constrainconf` module.
- Entrance point constraintconf now it is link to `moldrug.constrainconf.constraintconf_cmd` instead `moldrug.fitness.constraintconf_cmd`.

## [2.0.0] - 2022.08.25

### Added

- The functions `duplicate_conformers`, `get_mcs`, `generate_conformers`, `constraintconf` and `constraintconf_cmd` and the class `ProteinLigandClashFilter`. The code was borrowed from [Pat Walters](https://github.com/PatWalters/fragment_expansion/blob/master/rdkit_eval/rd_gen_restricted_confs.py). It is used if constraint docking is needed.
- `constraintconf` can be called from the command line.
- `moldrug.fitness.vinadock()` a simple wrapper around vina. This function will be used for all the implemented cost functions inside of the module `moldrug.fitness`. It could be used for constraint docking.
- `moldrug.data.constraintref`. This module is used for testing in case constraint docking is needed. It has two MolBlock strings: `r_6lu7` and `r_x0161`. That could be easily converted in RDKit molecules.

    ```python
    from rdkit import Chem
    from moldrug.data import constraintref
    mol = Chem.MolFromMolBlock(constraintref.r_x0161)
    ```

    This molecule is needed for the keyword argument `constraint_ref` of the functions of the `moldrug.fitness` module in case of constraint docking is used.
- Constraint docking capability in all implemented cost functions of the module `moldrug.fitness`.
- `moldrug.data.receptor_pdb`. This module is similar to `moldrug.data.receptor_pdbqt` but in pdb format.
- Documentation and tutorials.

### Changed

- `moldrug.utils.make_sdf` only will create the sdf file based on the `pdbqt` attribute. If `pdbqt` is a list, it will work as previous version works with `pdbqts` attribute.
- Name of the module `moldrug.data.receptors` to `moldrug.data.receptor_pdbqt`.
- Name of keyword argument `receptor_path` to `receptor_pdbqt_path` on the cost functions: `moldrug.fitness.Cost` and `moldrug.fitness.CostOnlyVina`.
- Name of keyword arguments `receptor_path`, `vina_score_types`, `boxcenters` and `boxsizes` to `receptor_pdbqt_path`, `vina_score_type`, `boxcenter` and `boxsize` respectively  on the cost functions: `moldrug.fitness.CostMultiReceptors` and `moldrug.fitness.CostMultiReceptorsOnlyVina`.
- `smiles` attribute in `moldrug.utils.Individual` now it is always without explicit Hs, despite if the mol attribute has them.

## [1.1.0] - 2022.08.23

### Changed

- `moldrug.utils.Individual` now is a hashable object.
- `moldrug.utils.GA.SawIndividuals` now is a `set` instead of a `list`
- `moldrug.utils.update_reactant_zone` sets the keywords `matchValences` and `ringMatchesRingOnly` to True on rdFMCS.FindMCS. This prevent undesired effects. E.g:

    ```python
    from moldrug.utils import update_reactant_zone
    from rdkit import Chem
    mol1 = Chem.MolFromSmiles('c1ccccc1')
    mol2 = Chem.MolFromSmiles('CCCN')
    update_reactant_zone(parent=mol1,offspring=mol2 parent_replace_ids = [0,2])
    ```

    Before the results was: `([0, 2, 3], [])`. Now it is: `([0, 1, 2, 3], [])`. The behavior was because the `ringMatchesRingOnly` is set to `False` by default inside RDKit.

### Removed

- `moldrug.utils.timeit`. No longer needed.

### Added

- Documentation.

## [1.0.2] - 2022.08.08

### Removed

- `Popen` option in `moldrug.utils.run`.

### Changed

- `RuntimeError` by `warnings.warn` when vina run fails and save every error as `idx_error.pbz2`. Where idx is the index of the failed individual.
- Print format when `mutate` fails inside of `moldrug.utils.GA`

### Added

- Print MolDrug's version when the command line is used.

## [1.0.0] - 2022.07.30

### Fixed

- Hidden `RuntimeError` in `moldrug.fitness` module.
- Bug printing number of generations in final info.

### Removed

- Unused code in `moldrug.home` module.
- 3D conformation in the `mol` attribute of the Individual during initialization.
- The use of `grow_mol` in the initialization of the populating when `get_similar = True`. Now the population is initialized with `mutate_mol` and the same set of crem parameters used during the searching.
- The automatic addition of Hs in the case where `min_size` and/or `max_size` were equal to zero. Now if your intention is work with the hydrogens, you must provided a SMILES with the explicit Hs. In the future the input could be just a RDKit mol. Now you must specify if you would like to add explicit Hs to the molecule withe keyword `AddHs`; default is False and is used for both `moldrug.utils.GA` and `moldrug.utils.Local`.

### Added

- Handling vina RuntimeError and keeping track for debug. This feature is used to identify what is the error. In the future will be removed.
- Two new fitness functions: `moldrug.fitness.CostOnlyVina` and `moldrug.fitness.CostMultiReceptorsOnlyVina`. They only use the information of vina scoring function. See the [docs](https://moldrug.readthedocs.io/en/latest/source/modules/fitness.html) for more info about it.
- Tracking of atom indexes during generations in order to use `protected_ids` and `replace_ids` options of `mutate_mol` function of [CReM](https://github.com/DrrDom/crem). Before it was not possible; the use of these features generate undesired solutions because the indexes are not static over generations. Even so, there are still some problems for symmetric molecules. We are working on it.

### Changed

- The whole MolDrug works base on the RDKit mol instead of the SMILES string:
    1. `moldrug.utils.Individual` is now initialized `mol` instead `smiles`. Now the SMILES string is generated internally, it still used as identifying for the instance.
    2. `moldrug.utils.GA` changed `smiles` for `seed_mol` and it is not needed the `mol` variable any more.
    3. `moldrug.utils.Local` changed `mol` for `seed_mol` in the initialization variables..
    4. `moldrug.utils.confgen` changed smiles for `mol` variable.


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

[unreleased]: https://github.com/ale94mleon/MolDrug/compare/3.4.0...HEAD
[3.4.0]: https://github.com/ale94mleon/MolDrug/compare/3.3.0...3.4.0
[3.3.0]: https://github.com/ale94mleon/MolDrug/compare/3.2.5...3.3.0
[3.2.5]: https://github.com/ale94mleon/MolDrug/compare/3.2.2...3.2.5
[3.2.2]: https://github.com/ale94mleon/MolDrug/compare/3.2.0...3.2.2
[3.2.0]: https://github.com/ale94mleon/MolDrug/compare/3.1.0...3.2.0
[3.1.0]: https://github.com/ale94mleon/MolDrug/compare/3.0.3...3.1.0
[3.0.3]: https://github.com/ale94mleon/MolDrug/compare/3.0.1...3.0.3
[3.0.1]: https://github.com/ale94mleon/MolDrug/compare/3.0.0...3.0.1
[3.0.0]: https://github.com/ale94mleon/MolDrug/compare/2.1.12...3.0.0
[2.1.12]: https://github.com/ale94mleon/MolDrug/compare/2.1.7...2.1.12
[2.1.7]: https://github.com/ale94mleon/MolDrug/compare/2.1.0...2.1.7
[2.1.0]: https://github.com/ale94mleon/MolDrug/compare/2.0.0...2.1.0
[2.0.0]: https://github.com/ale94mleon/MolDrug/compare/1.1.0...2.0.0
[1.1.0]: https://github.com/ale94mleon/MolDrug/compare/1.0.2...1.1.0
[1.0.2]: https://github.com/ale94mleon/MolDrug/compare/1.0.0...1.0.2
[1.0.0]: https://github.com/ale94mleon/MolDrug/compare/0.1.0...1.0.0
[0.1.0]: https://github.com/ale94mleon/MolDrug/compare/0.0.4...0.1.0
[0.0.4]: https://github.com/ale94mleon/MolDrug/compare/0.0.3...0.0.4
