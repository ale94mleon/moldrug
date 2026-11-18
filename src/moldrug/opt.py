#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Optimization classes: the genetic algorithm GA and the local search Local."""

import datetime
import logging
import os
import random
import time
from copy import deepcopy
from typing import Callable, Dict, Iterable, List, Optional, Union

import numpy as np
from crem.crem import grow_mol, mutate_mol
from rdkit import Chem, RDLogger

from moldrug import __version__
from moldrug.runner import Runner, RunnerMode
from moldrug.utils import (Individual, _make_kwargs_copy, compressed_pickle,
                           format_duration, full_pickle, get_similar_mols,
                           is_iter, make_sdf, roulette_wheel_selection, softmax,
                           tar_errors, to_dataframe, update_reactant_zone)

logger = logging.getLogger(__name__)

#: How many replacements CReM is asked to sample per mutation, ``None`` meaning "enumerate
#: all of them". Capping it saves the work of building mutants that get discarded, the more
#: so the larger the database, but it is not only a speed setting: CReM gives each fragment
#: of the molecule an equal quota, however many
#: replacements the database holds for it, so sampling re-weights which part of the molecule
#: gets modified, away from the groups the database knows well. See the CReM sampling page
#: of the documentation for the measurements behind this default.
DEFAULT_MAX_REPLACEMENTS = None

RDLogger.DisableLog('rdApp.*')


class Local:
    """Generate close solutions to the seed molecule.

    It use :meth:`crem.crem.grow_mol`.

    Attributes
    ----------
    randomseed : Union[None, int]
        The seed of :attr:`rng`.
    rng : random.Random
        The generator every random decision of this class is drawn from. It is seeded with
        ``randomseed`` and pickled with the class, so a run restored from a checkpoint
        carries on where it left off.
    __moldrug_version__ : str
        The molDrug version.
    costfunc : object
        The cost function set by the user.
    crem_db_path : str
        Path to the CReM data base.
    costfunc_kwargs : dict
        The keyword arguments of the costfunc.
    grow_crem_kwargs : dict
        The keyword arguments to pass to :meth:`crem.crem.grow_mol`.
    timing : dict
        With keys ``generation`` and ``evaluation``, how many seconds have been spent
        building molecules and calculating their cost. Both accumulate over calls, while
        the summary at the end of a call reports that call alone. They do not add up to
        the total time of a call, since cleaning the working directory falls outside both.
    AddHs : bool
        In case explicit hydrogens should be added.
    pop : list[:meth:`moldrug.Individuals`]
        The final population sorted by cost.
    """

    def __init__(self, seed_mol: Chem.rdchem.Mol, crem_db_path: str, costfunc: object, grow_crem_kwargs: Dict = None,
                 costfunc_kwargs: Dict = None, AddHs: bool = False, randomseed: Union[None, int] = None,
                 deffnm: str = 'local') -> None:
        """Creator.

        Parameters
        ----------
        seed_mol : Chem.rdchem.Mol
            The seed molecule from which the population will be generated.
        crem_db_path : str
            The pathway to the CReM data base.
        costfunc : object
            The cost function to work with (any from :mod:`moldrug.fitness` or a valid user defined).
        grow_crem_kwargs : Dict, optional
            The keywords of the grow_mol function of CReM, by default None
        costfunc_kwargs : Dict, optional
            The keyword arguments of the selected cost function, by default None
        AddHs : bool, optional
            If True the explicit hyrgones will be added, by default False
        randomseed : Union[None, int], optional
           Set a random seed for reproducibility, by default None
        deffnm : str
            Just a place holder for compatibility with the CLI.

        Raises
        ------
        Exception
            In case that some problem occured during the creation of the Individula from the seed_mol
        ValueError
            In case of incorrect definition of grow_crem_kwargs and/or costfunc_kwargs.
            They must be None or a dict instance.
        """
        self.randomseed = randomseed
        self.rng = random.Random(randomseed)  # nosec

        self.__moldrug_version = __version__

        if grow_crem_kwargs is None:
            grow_crem_kwargs = dict()
        elif not isinstance(grow_crem_kwargs, dict):
            raise ValueError(f'grow_crem_kwargs must be None or a dict instance. {grow_crem_kwargs} was provided')

        if costfunc_kwargs is None:
            costfunc_kwargs = dict()
        elif not isinstance(costfunc_kwargs, dict):
            raise ValueError(f'grow_crem_kwargs must be None or a dict instance. {costfunc_kwargs} was provided')

        if AddHs:
            self._seed_mol = Chem.AddHs(seed_mol)
        else:
            self._seed_mol = seed_mol

        self.InitIndividual = Individual(self._seed_mol, idx=0, randomseed=self.randomseed)
        if not self.InitIndividual.pdbqt:
            raise Exception("For some reason, it was not possible to create for the class Individula "
                            "a pdbqt from the seed_smiles. Consider to check the validity of the SMILES string!")
        if os.path.isfile(crem_db_path):
            self.crem_db_path = os.path.abspath(crem_db_path)
        else:
            raise FileNotFoundError(f"{crem_db_path = } does not exists or is not accesible")

        self.grow_crem_kwargs = grow_crem_kwargs
        self.costfunc = costfunc
        self.costfunc_kwargs = costfunc_kwargs
        self.timing = {'generation': 0.0, 'evaluation': 0.0}
        self.pop = [self.InitIndividual]

    def __call__(self, njobs: int = 1, pick: int = None, runner: Optional[Runner] = None):
        """Call definition.

        Parameters
        ----------
        njobs : int, optional
            The number of jobs for parallelization, the module multiprocessing will be used, by default 1
        pick : int, optional
            How many molecules take from the generated throgh the grow_mol CReM operation,
            by default None which means all generated.
        runner: Runner, optional
            Providing this parameter instead of njobs allows execution on a cluster with Dask. Few other modes
            of executions are also available, but useful mostly for debugging and profiling.
        """
        if njobs > 1 and runner is not None:
            raise ValueError("Both njobs > 1 and runner have been specified. Please use only one of the parameters.")

        if runner is None:
            runner = Runner(RunnerMode.MULTIPROCESSING, process_count=njobs)

        # Check version of moldrug
        if self.__moldrug_version != __version__:
            logger.error(f"{self.__class__.__name__} was initilized with moldrug-{self.__moldrug_version} "
                         f"but was called with moldrug-{__version__}")
        ts = time.time()
        tg = time.perf_counter()
        self.grow_crem_kwargs.update({'return_mol': True})
        new_mols = list(grow_mol(self._seed_mol, self.crem_db_path, **self.grow_crem_kwargs))
        if pick:
            self.rng.shuffle(new_mols)
            new_mols = new_mols[:pick]
            new_mols = [item[1] for item in new_mols]

        idx0 = len(self.pop)
        for i, mol in enumerate(new_mols):
            individual = Individual(mol, idx=idx0 + i, randomseed=self.randomseed)
            if individual.pdbqt:
                self.pop.append(individual)
        generation = time.perf_counter() - tg
        self.timing['generation'] += generation

        # Calculating cost of each individual
        # Creating the arguments
        te = time.perf_counter()
        args_list = []
        # Make a copy of the self.costfunc_kwargs
        kwargs_copy, costfunc_jobs_tmp_dir = _make_kwargs_copy(self.costfunc, self.costfunc_kwargs)

        for individual in self.pop:
            args_list.append((individual, kwargs_copy))

        logger.info('Calculating cost function...')
        self.pop = runner.run(self.__costfunc__, args_list)

        # Clean directory
        costfunc_jobs_tmp_dir.cleanup()
        evaluation = time.perf_counter() - te
        self.timing['evaluation'] += evaluation
        # Tar errors
        tar_errors('error')

        # Printing how long was the simulation, and where the time went. The two phases do
        # not add up to the total: cleaning the working directory is outside both.
        total = time.time() - ts
        logger.info(f"Total time: {format_duration(total)}.\n"
                    f"Generating molecules: {format_duration(generation)} "
                    f"({100 * generation / total:.1f} %).\n"
                    f"Evaluating the cost function: {format_duration(evaluation)} "
                    f"({100 * evaluation / total:.1f} %).\n"
                    f"Finished at {datetime.datetime.now().strftime('%c')}.\n")

    def __costfunc__(self, args_list):
        """Evaluate the cost function on one individual, as the runner expects it."""
        Individual, kwargs = args_list
        # This is just to use the progress bar on pool.imap
        return self.costfunc(Individual, **kwargs)

    def pickle(self, title: str, compress: bool = False):
        """Pickle the whole Local class.

        Parameters
        ----------
        title : str
            Name of the object which will be compleated with the correposnding
            extension depending if compress is set to True or False.
        compress : bool, optional
            Use compression, by default False. If True :meth:`moldrug.compressed_pickle` will be used;
            if not :meth:`moldrug.full_pickle` will be used instead.
        """
        cls = self.__class__
        result = cls.__new__(cls)
        result.__dict__.update(self.__dict__)
        if compress:
            compressed_pickle(title, result)
        else:
            full_pickle(title, result)

    def to_dataframe(self, return_mol: bool = False):
        """Create a DataFrame from self.pop.

        Returns
        -------
        pandas.DataFrame
            The DataFrame
        """
        return to_dataframe(self.pop, return_mol=return_mol)


#######################
# Optimazer functions
#######################

class GA:
    """An implementation of a genetic algorithm to search in the chemical space.

    Attributes
    ----------
    randomseed : Union[None, int]
        The seed of :attr:`rng`.
    rng : random.Random
        The generator every random decision of this class is drawn from. It is seeded with
        ``randomseed`` and pickled with the class, so a run restored from a checkpoint
        carries on where it left off.
    __moldrug_version__ : str
        The moldrug version.
    costfunc : object
        The cost function set by the user.
    crem_db_path : str
        Path to the CReM data base.
    maxiter : int
        Maximum number of iteratinos to perform.
    popsize : int
        Population size.
    beta : float
        Selection pressure.
    costfunc_kwargs : dict
        The keyword arguments of the costfunc.
    costfunc_ncores : int
        The number of cores to use for costfunc.
    nc : int
        Number of childs of offsprints ``= round(pc * popsize)``
    get_similar : bool
        Bias the search upon similar molecules. If True :meth:`modrug.get_similar_mols` is used after the
        mutation with CReM instead random choice.
    mutate_crem_kwargs : dict
        The keyword arguments to pass to :meth:`crem.crem.mutate_mol`.
    save_pop_every_gen : int
        Frequency to save the pickle file o fthe population during the optimazation.
    checkpoint : bool
        Safe chekpoint file, this help to restart a simualation.
    deffnm : str
        Prefix for the genereated files.
    NumCalls : int
        How many times the ``__call__`` method has been called.
    NumGens : int
        he number of generations performed by the class. Subsequent ``__call__``
        executions update this number acordennly.
    SawIndividuals : set[:meth:`moldrug.Individuals`]
        All the Individulas saw during the optimizations.
    acceptance : dict
        A dictionary with key the Generation id and as value another dictionary
        with keys ``accepeted`` and ``generated`` with the number of accepted and genereated
        individuals on the generation respectively.
    timing : dict
        A dictionary with key the Generation id and as value another dictionary with keys
        ``generation`` and ``evaluation``, holding how many seconds that generation spent
        building molecules and calculating their cost. Generation 0 is the initialization
        of the population. The two do not add up to the total time of a call: saving,
        sorting and checkpointing are outside both.
    AddHs : bool
        In case explicit hydrogens should be added for all genreated molecules.
    _seed_mol : list[Chem.rdchem.Mol]
        The list of seed molecules.
    InitIndividual : :meth:`moldrug.Individuals`
        The initial individual based on _seed_mol.
    pop : list[:meth:`moldrug.Individuals`]
        The final population sorted by cost.
    best_cost : list[float]
        The list of best cost for each generations.
    avg_cost : list[float]
        The list of average cost for each generations.

    TODO:

        * Extend to other generators:

            - mutate_crem_kwargs = None and some other keyword that gets the generator
              function, in this case the mutate method will be overwritten with the one the
              user provided. This function will take an Individual and return a new
              offspring. To be more compatible and not create issues, a good idea will be
              that this function accepts a self as first argument, and internally it will
              use the self of the GA class

    """

    def __init__(self, seed_mol: Union[Chem.rdchem.Mol, Iterable[Chem.rdchem.Mol]],
                 costfunc: Callable, costfunc_kwargs: Dict, crem_db_path: str, maxiter: int = 10, popsize: int = 20,
                 beta: float = 0.001, pc: float = 1, get_similar: bool = False, mutate_crem_kwargs: Union[None, Dict] = None,
                 save_pop_every_gen: int = 0, checkpoint: bool = False, deffnm: str = 'ga',
                 AddHs: bool = False, randomseed: Union[None, int] = None) -> None:
        """Create the GA.

        Parameters
        ----------
        seed_mol : Union[Chem.rdchem.Mol, Iterable[Chem.rdchem.Mol]]
            The seed molecule submitted to genetic algorithm optimization on the chemical space. Could be only one RDKit
            molecule or more than one specified in an Iterable object.
        costfunc : Callable
            The cost function to work with (any from :mod:`moldrug.fitness` or a valid user defined).
        costfunc_kwargs : Dict
            The keyword arguments of the selected cost function
        crem_db_path : str
            Path to the CReM data base.
        maxiter : int, optional
            Maximum number of iteration (or generation), by default 10.
        popsize : int, optional
            Population size, by default 20.
        beta : float, optional
            Selection pressure. Higher values means that the best individual
            are going to be sumitted for mutations more frquently, by default 0.001.
        pc : float, optional
            Proportion of children, by default 1
        get_similar : bool, optional
            If True the searching will be bias to similar molecules, by default False
        mutate_crem_kwargs : Union[None, Dict], optional
            Parameters for mutate_mol of CReM, by default {}. ``max_replacements`` is worth
            knowing about: it caps how many mutants CReM builds per mutation, which is much
            faster on a large database, but it also changes how the chemical space is
            explored. See the CReM sampling page of the documentation.
        save_pop_every_gen : int, optional
            Frequency to save the population, by default 0
        checkpoint : bool, optional
            If True the whole class will be saved as cpt with the frequency of save_pop_every_gen.
            This means that if save_pop_every_gen = 0 and checkpoint = True, no checkpoint will be
            output, by default False
        deffnm : str, optional
            Default prefix name for all generated files, by default 'ga'
        AddHs : bool, optional
           If True the explicit hydrogens will be added, by default False
        randomseed : Union[None, int], optional
           Set a random seed for reproducibility, by default None
        Raises
        ------
        TypeError
            In case that seed_mol is a wrong input.
        ValueError
            In case of incorrect definition of mutate_crem_kwargs. It must be None or a dict instance.
        ValueError
            In case of crem_db_path deos not exist.
        """
        self.randomseed = randomseed
        self.rng = random.Random(randomseed)  # nosec

        self.__moldrug_version__ = __version__
        if mutate_crem_kwargs is None:
            mutate_crem_kwargs = dict()
        elif not isinstance(mutate_crem_kwargs, dict):
            raise ValueError(f'mutate_crem_kwargs must be None or a dict instance. {mutate_crem_kwargs} was provided')

        self.costfunc = costfunc
        if os.path.exists(crem_db_path):
            self.crem_db_path = os.path.abspath(crem_db_path)
        else:
            raise FileNotFoundError(f"{crem_db_path = } does not exists or is not accesible")

        self.maxiter = maxiter
        self.popsize = popsize
        self.beta = beta
        self.costfunc_kwargs = costfunc_kwargs
        if 'ncores' in costfunc_kwargs:
            self.costfunc_ncores = self.costfunc_kwargs['ncores']
        else:
            self.costfunc_ncores = 1

        self.nc = round(pc * popsize)
        self.get_similar = get_similar
        # Internally update with default values, TODO, maybe I should remove this
        # and the users should handled CReM parameters by themself.
        # I think that this will avoid confusion.
        self.mutate_crem_kwargs = {
            'radius': 3,
            'min_size': 1,
            'max_size': 8,
            'min_inc': -5,
            'max_inc': 3,
            'ncores': 1,
            'max_replacements': DEFAULT_MAX_REPLACEMENTS,
        }
        self.mutate_crem_kwargs.update(mutate_crem_kwargs)
        # We need the molecule back, so this overrides any user definition of the keyword
        self.mutate_crem_kwargs['return_mol'] = True

        # Saving parameters
        self.save_pop_every_gen = save_pop_every_gen
        self.checkpoint = checkpoint
        self.deffnm = deffnm

        # Tracking parameters
        self.NumCalls = 0
        self.NumGens = 0
        self.SawIndividuals = set()
        self.acceptance = dict()
        self.timing = dict()

        # work with the seed molecule or population
        self.AddHs = AddHs

        # Convert to list the seed_mol in case that it is not
        if not is_iter(seed_mol) and isinstance(seed_mol, Chem.rdchem.Mol):
            self._seed_mol = [seed_mol]
        elif all([isinstance(mol, Chem.rdchem.Mol) for mol in seed_mol]):
            self._seed_mol = seed_mol
        else:
            raise TypeError("seed_mol is not Chem.rdchem.Mol neither a Iterable[Chem.rdchem.Mol]")

        if self.AddHs:
            self._seed_mol = [Chem.AddHs(mol) for mol in self._seed_mol]
        # if 'protected_ids' in self.mutate_crem_kwargs or 'replace_ids' in self.mutate_crem_kwargs:
        #     _ = [atom.SetIntProp('label_moldrug', atom.GetIdx()) for atom in seed_mol.GetAtoms()]

        # Create the first Individual
        self.InitIndividual = Individual(self._seed_mol[0], idx=0, randomseed=self.randomseed)
        self.pop = []

    def __call__(self, njobs: int = 1, runner: Optional[Runner] = None):
        """Call definition.

        Parameters
        ----------
        njobs : int, optional
            The number of jobs for parallelization, the module multiprocessing will be used, by default 1,
        runner: Runner, optional
            Providing this parameter instead of njobs allows execution on a cluster with Dask. Few other modes
            of executions are also available, but useful mostly for debugging and profiling.

        Raises
        ------
        ValueError
            If both njobs > 1 and runner were provided.
        RuntimeError
            Error during the initialization of the population.
        """
        if njobs > 1 and runner is not None:
            raise ValueError("Both njobs > 1 and runner have been specified. Please use only one of the parameters.")

        if runner is None:
            runner = Runner(RunnerMode.MULTIPROCESSING, process_count=njobs)

        ts = time.time()
        # Earlier calls have already filled self.timing, and the summary at the end of this
        # one reports this call alone.
        spent_before = {phase: sum(gen[phase] for gen in self.timing.values())
                        for phase in ('generation', 'evaluation')}
        # Counting the calls
        self.NumCalls += 1

        # Check version of moldrug
        if self.__moldrug_version__ != __version__:
            logger.error(f"{self.__class__.__name__} was initialized with moldrug-{self.__moldrug_version__} "
                         f"but was called with moldrug-{__version__}")

        # In case that the population exists there is not need to initialize.
        if len(self.pop) == 0:
            self._initialize_population(runner)

        # Saving tracking variables, the first population, outside the if to take into account second calls
        # with different population provided by the user.
        self.SawIndividuals.update(self.pop)

        # Saving population in disk if it was required
        if self.save_pop_every_gen:
            self._save_population(sorted(self.pop))

        # Main Loop
        # Another control variable. In case that the __call__ method is used more than ones.
        number_of_previous_generations = len(self.best_cost)
        for it in range(self.maxiter):
            # Saving Number of Generations
            self.NumGens += 1

            popc = self._generate_offspring()

            if popc:  # Only if there are new members
                # Add idx label to each individual
                NumbOfSawIndividuals = len(self.SawIndividuals)
                for (i, individual) in enumerate(popc):
                    individual.idx = i + NumbOfSawIndividuals
                popc = self._evaluate_costs(
                    popc, runner,
                    f'Evaluating generation {self.NumGens} / {self.maxiter + number_of_previous_generations}:')

            self._update_population(popc)

            # Saving population in disk if it was required
            if self.save_pop_every_gen:
                # Save every save_pop_every_gen and always the last population
                if self.NumGens % self.save_pop_every_gen == 0 or it + 1 == self.maxiter:
                    self._save_population(self.pop)

            # Show Iteration Information
            logger.info(f"Generation {self.NumGens}: Best Individual: {self.pop[0]}")
            logger.info(f"Acceptance rate: {self.acceptance[self.NumGens]['accepted']} / {self.acceptance[self.NumGens]['generated']}\n")

        self._log_summary(ts, spent_before)

    def _crem_kwargs(self):
        """Build the CReM keywords for one call, with a fresh sampling seed.

        When ``max_replacements`` is set, CReM samples the replacements with a generator of
        its own, ``random.Random(seed)``, which does not read :attr:`rng`. Leaving ``seed``
        unset would therefore make a run with ``randomseed`` unreproducible, and pinning it
        to a constant would hand back the same molecule on every mutation of the same
        parent. Drawing the seed from :attr:`rng` gives both: it varies within a run and
        repeats across runs.

        A ``seed`` the user put in ``mutate_crem_kwargs`` is left alone.

        Returns
        -------
        dict
            A copy of ``mutate_crem_kwargs``, possibly with ``seed`` added.
        """
        kwargs = self.mutate_crem_kwargs.copy()
        if kwargs.get('max_replacements') is not None and 'seed' not in kwargs:
            kwargs['seed'] = self.rng.randint(0, 2**31 - 1)
        return kwargs

    def _initialize_population(self, runner: Runner):
        """Build the first population and evaluate it.

        The seed molecules are kept and the rest of the population is generated from them
        with CReM, then the whole population is costed and the bookkeeping that the main
        loop relies on is set up.

        Parameters
        ----------
        runner : Runner
            The runner used to evaluate the cost function.

        Raises
        ------
        RuntimeError
            If the seed molecules did not generate any new molecule.
        """
        GenInitStructs = []
        # in case that the input has the popsize memebers there is not need to generate new structures
        if len(self._seed_mol) < self.popsize:
            for mol in self._seed_mol:
                tmp_GenInitStructs = list(mutate_mol(mol, self.crem_db_path, **self._crem_kwargs()))
                tmp_GenInitStructs = [mol for (_, mol) in tmp_GenInitStructs]
                GenInitStructs += tmp_GenInitStructs
            # Checking for possible scenarios
            if len(GenInitStructs) == 0:
                raise RuntimeError("Something really strange happened. The seed_mol did not "
                                   "generate any new molecule during the initialization of the population. "
                                   "Check the provided crem parameters!")
            if len(GenInitStructs) < (self.popsize - len(self._seed_mol)):
                logger.warning('The initial population has repeated elements')
                # temporal solution
                GenInitStructs += self.rng.choices(GenInitStructs,
                                                   k=self.popsize - len(GenInitStructs) - len(self._seed_mol))
            elif len(GenInitStructs) > (self.popsize - 1):
                # Selected random sample from the generation
                GenInitStructs = self.rng.sample(GenInitStructs, k=self.popsize - len(self._seed_mol))
            else:
                # Everything is ok!
                pass

        # Adding the inputs to the initial population
        for i, mol in enumerate(self._seed_mol):
            individual = Individual(mol, idx=i, randomseed=self.randomseed)
            if individual.pdbqt:
                self.pop.append(individual)

        # Completing the population with the generated structures
        for i, mol in enumerate(GenInitStructs):
            if self.AddHs:
                individual = Individual(Chem.AddHs(mol), idx=i + len(self._seed_mol), randomseed=self.randomseed)
            else:
                individual = Individual(mol, idx=i + len(self._seed_mol), randomseed=self.randomseed)
            if individual.pdbqt:
                self.pop.append(individual)

        # Make sure that the population do not have more than popsize members and it is without repeated elements.
        # That could happens if seed_mol has more molecules than popsize
        self.pop = sorted(set(self.pop), key=lambda x: x.idx)[:self.popsize]

        self.pop = self._evaluate_costs(
            self.pop, runner, f'Creating the first population with {len(self.pop)} members:')

        # Adding generation information
        for individual in self.pop:
            individual.genID = self.NumGens
            individual.kept_gens = set([self.NumGens])

        self.acceptance[self.NumGens] = {
            'accepted': len(self.pop[:]),
            'generated': len(self.pop[:])
        }

        # Get the same order population in case cost is the same. Sorted by idx and then by cost
        if self.randomseed:
            self.pop = sorted(self.pop, key=lambda x: x.idx)
        self.pop = sorted(self.pop)
        # Print some information of the initial population
        logger.info(f"Initial Population: Best Individual: {self.pop[0]}")
        logger.info(f"Acceptance rate: {self.acceptance[self.NumGens]['accepted']} / {self.acceptance[self.NumGens]['generated']}\n")
        # Updating the info of the first individual (parent)
        # to print at the end how well performed the method (cost function)
        # Because How the population was initialized and because we are using pool.imap (ordered).
        # The parent is the first Individual of self.pop.
        # We have to use deepcopy because Individual is a mutable object
        # Because above set were used, we have to sorter based on idx
        self.InitIndividual = deepcopy(
            min(
                sorted(self.pop, key=lambda x: x.idx)[:len(self._seed_mol)]
            )
        )
        # Best Cost of Iterations
        self.best_cost = []
        self.avg_cost = []

    def _generate_offspring(self):
        """Select parents and mutate them into the offspring of the current generation.

        Only offspring that are new, not already in this generation, and have a valid pdbqt
        are kept, so fewer than ``nc`` individuals may come back.

        Returns
        -------
        list[Individual]
            The offspring, without their cost calculated yet.
        """
        ts = time.perf_counter()
        # Probabilities Selections
        probs = softmax((-self.beta * np.array(self.pop)).astype('float64'))
        if any(np.isnan(probs)):
            probs = np.nan_to_num(probs)
        # TODO: This cycle should run in this way only if no user generetor was provided
        # with and if, else statment I could correct, and then the genereator functions is completlly up to the user,
        # then I do not need to worry in how the selection is made,
        # In this case self.nc will not have any validity unless the user use it with its evaluator
        # the checking of SawIndivduals must be done after the user funcrion return the popc
        # The other that I need to change is that if the if it is a new genereator the genereation of the initil population is different
        # the other is that checking for redundancy may be complicated in the case, that molecules are, for example peptides,
        # in this case other identifier like the aa sequnce should be ued intead. For that the user may need a different Individual instance
        # a one more efficient, there are a lot of if here :`-)
        popc = []
        for _ in range(self.nc):
            # Perform Roulette Wheel Selection
            parent = self.pop[roulette_wheel_selection(probs, self.rng)]

            # Perform Mutation (this mutation is some kind of crossover but with CReM library)
            children = self.mutate(parent)

            # Save offspring population
            # I will save only those offsprings that were not seen and that have a correct pdbqt file
            if children not in self.SawIndividuals and children not in popc and children.pdbqt:
                children.genID = self.NumGens
                children.kept_gens = set()
                popc.append(children)
        self._record_time('generation', time.perf_counter() - ts)
        return popc

    def _evaluate_costs(self, individuals: List[Individual], runner: Runner, message: str):
        """Calculate the cost of every individual, in parallel if the runner allows it.

        Parameters
        ----------
        individuals : list[Individual]
            The individuals to evaluate.
        runner : Runner
            The runner that spreads the evaluations across workers.
        message : str
            Logged just before the evaluations start.

        Returns
        -------
        list[Individual]
            The evaluated individuals, in the order they were given.
        """
        ts = time.perf_counter()
        # Make a copy of the self.costfunc_kwargs
        kwargs_copy, costfunc_jobs_tmp_dir = _make_kwargs_copy(self.costfunc, self.costfunc_kwargs)

        # The problem here is that we are not being general for other possible Cost functions.
        args_list = [(individual, kwargs_copy) for individual in individuals]

        logger.info(message)
        evaluated = runner.run(self.__costfunc__, entries=args_list)

        # Clean directory
        costfunc_jobs_tmp_dir.cleanup()
        self._record_time('evaluation', time.perf_counter() - ts)
        return evaluated

    def _update_population(self, popc: List[Individual]):
        """Merge the offspring into the population and keep the best ``popsize`` of them.

        Also updates the acceptance rate, the cost history and the set of individuals seen
        so far.

        Parameters
        ----------
        popc : list[Individual]
            The evaluated offspring of the current generation.
        """
        # Merge, Sort and Select
        self.pop += popc
        if self.randomseed:
            self.pop = sorted(self.pop, key=lambda x: x.idx)
        self.pop = sorted(self.pop)
        self.pop = self.pop[:self.popsize]

        # Update the kept_gens attribute
        self.acceptance[self.NumGens] = {
            'accepted': 0,
            'generated': len(popc)
        }
        for individual in self.pop:
            if not individual.kept_gens:
                self.acceptance[self.NumGens]['accepted'] += 1
            individual.kept_gens.add(self.NumGens)

        # Store Best Cost
        self.best_cost.append(self.pop[0].cost)

        # Store Average cost
        self.avg_cost.append(np.mean(self.pop))

        # Saving tracking variables
        self.SawIndividuals.update(popc)

    def _save_population(self, population: List[Individual]):
        """Write the population to disk, and the whole class if checkpointing is on.

        Parameters
        ----------
        population : list[Individual]
            The population to save.
        """
        compressed_pickle(f"{self.deffnm}_pop", (self.NumGens, population))
        make_sdf(population, sdf_name=f"{self.deffnm}_pop")
        if self.checkpoint:
            compressed_pickle('cpt', self)

    def _record_time(self, phase: str, seconds: float):
        """Add ``seconds`` to one phase of the current generation of :attr:`timing`.

        Parameters
        ----------
        phase : str
            Either ``'generation'`` or ``'evaluation'``.
        seconds : float
            How long the phase took.
        """
        self.timing.setdefault(self.NumGens, {'generation': 0.0, 'evaluation': 0.0})
        self.timing[self.NumGens][phase] += seconds

    def _log_summary(self, ts: float, spent_before: Dict[str, float]):
        """Report how the simulation went, and how long it took.

        Parameters
        ----------
        ts : float
            The value of ``time.time()`` when the call started.
        spent_before : dict[str, float]
            How many seconds each phase of :attr:`timing` already held when the call
            started, so that a second call reports its own time and not the whole history.
        """
        # Printing summary information
        logger.info(f"\t\t{20*'=+'}\n")
        logger.info(f"The simulation finished successfully after {self.NumGens} generations with"
                    f"a population of {self.popsize} individuals. "
                    f"A total number of {len(self.SawIndividuals)} Individuals were seen during the simulation.")
        logger.info(f"Initial Individual: {self.InitIndividual}")
        logger.info(f"Final Individual: {self.pop[0]}")
        logger.info(f"The cost function dropped in {self.InitIndividual - self.pop[0]} units.")
        logger.info(f"\t\t{20*'=+'}\n")

        # Tar errors
        tar_errors('error')

        # Printing how long was the simulation, and where the time went. The two phases do
        # not add up to the total: saving the population, sorting and the checkpoints are
        # outside both of them.
        total = time.time() - ts
        generation = sum(gen['generation'] for gen in self.timing.values()) - spent_before['generation']
        evaluation = sum(gen['evaluation'] for gen in self.timing.values()) - spent_before['evaluation']
        logger.info(f"Total time ({self.maxiter} generations): {format_duration(total)}.\n"
                    f"Generating molecules: {format_duration(generation)} "
                    f"({100 * generation / total:.1f} %).\n"
                    f"Evaluating the cost function: {format_duration(evaluation)} "
                    f"({100 * evaluation / total:.1f} %).\n"
                    f"Finished at {datetime.datetime.now().strftime('%c')}.\n")

    def __costfunc__(self, args_list):
        """Evaluate the cost function on one individual, as the runner expects it."""
        individual, kwargs = args_list
        # This is just to use the progress bar on pool.imap
        return self.costfunc(individual, **kwargs)

    def mutate(self, individual: Individual):
        """Genetic operators.

        Parameters
        ----------
        individual : Individual
            The individual to mutate.

        Returns
        -------
        Individual
            A new Individual.
        """
        # Here is were I have to check if replace_ids or protected_ids where provided.
        mutate_crem_kwargs_to_work_with = self._crem_kwargs()
        if 'replace_ids' in self.mutate_crem_kwargs and 'protected_ids' in self.mutate_crem_kwargs:
            mutate_crem_kwargs_to_work_with['replace_ids'], mutate_crem_kwargs_to_work_with['protected_ids'] = update_reactant_zone(
                self.InitIndividual.mol, individual.mol, parent_replace_ids=self.mutate_crem_kwargs['replace_ids'],
                parent_protected_ids=self.mutate_crem_kwargs['protected_ids'])
        elif 'replace_ids' in self.mutate_crem_kwargs:
            mutate_crem_kwargs_to_work_with['replace_ids'], _ = update_reactant_zone(
                self.InitIndividual.mol, individual.mol, parent_replace_ids=self.mutate_crem_kwargs['replace_ids'])
        elif 'protected_ids' in self.mutate_crem_kwargs:
            _, mutate_crem_kwargs_to_work_with['protected_ids'] = update_reactant_zone(
                self.InitIndividual.mol, individual.mol,
                parent_protected_ids=self.mutate_crem_kwargs['protected_ids'])

        try:
            mutants = list(mutate_mol(individual.mol, self.crem_db_path, **mutate_crem_kwargs_to_work_with))
            if not mutants:
                raise RuntimeError('CReM did not generate any mutant, check the provided CReM parameters')
            # Bias the searching to similar molecules
            if self.get_similar:
                mol = get_similar_mols(mols=[mol for _, mol in mutants],
                                       ref_mol=self.InitIndividual.mol, pick=1, rng=self.rng, beta=0.01)[0]
            else:
                _, mol = self.rng.choice(mutants)
        except Exception as error:
            # Report what went wrong. CReM takes **kwargs, so a mistake in
            # mutate_crem_kwargs arrives here as a ValueError naming the unknown column,
            # which is indistinguishable from a failed mutation without the message.
            logger.warning(f'The mutation on {individual} did not work '
                           f'({type(error).__name__}: {error}), it will be returned the same individual')
            mol = individual.mol
        if self.AddHs:
            mol = Chem.AddHs(mol)
        return Individual(mol, randomseed=self.randomseed)

    def pickle(self, title: str, compress: bool = False):
        """Pickle the whole GA class.

        Parameters
        ----------
        title : str
            Name of the object which will be completed with the corresponding
            extension depending if compress is set to True or False.
        compress : bool, optional
            Use compression, by default False. If True :meth:`moldrug.compressed_pickle` will be used;
            if not :meth:`moldrug.full_pickle` will be used instead.
        """
        cls = self.__class__
        result = cls.__new__(cls)
        result.__dict__.update(self.__dict__)
        if compress:
            compressed_pickle(title, result)
        else:
            full_pickle(title, result)

    def to_dataframe(self, return_mol: bool = False):
        """Create a DataFrame from self.SawIndividuals.

        Returns
        -------
        pandas.DataFrame
            The DataFrame
        """
        return to_dataframe(self.SawIndividuals, return_mol=return_mol)
