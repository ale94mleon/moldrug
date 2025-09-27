#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import bz2
import collections.abc
import os
import random
import shutil
import subprocess
import tempfile
from copy import deepcopy
from inspect import signature
from typing import List, Union

import dill as pickle
import numpy as np
import pandas as pd
from meeko import (MoleculePreparation, PDBQTMolecule, PDBQTWriterLegacy,
                   RDKitMolCreate)
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, DataStructs, Descriptors, Lipinski, rdFMCS

from moldrug.logging_utils import LogLevel, log

RDLogger.DisableLog('rdApp.*')
# # in order to pickle the isotope properties of the molecule
# Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.AllProps)


class Individual:
    """
    Base class to work with GA, Local and all the fitness functions.
    Individual is a mutable object. Only the attribute smiles it is not mutable and is used
    for hash. Therefore this class is hashable based on the smiles attribute.
    This one is also used for '==' comparison
    If two Individuals has the same smiles not matter if the rest of the elements are different,
    they will be considered the same.
    The cost attribute is used for arithmetic operations.
    It also admit copy and deepcopy operations.
    Known issue, in case that we would like to use a numpy array of individuals.
    It is needed to change the dtype of the generated arrays

    Attributes
    ----------
    mol:  Chem.rdchem.Mol
        The molecule object
    idx: Union[int, str]
        The identifier
    pdbqt: str
        A pdbqt string representation of the molecule, used for docking with Vina. It is generated
        during the initialization of the class
    smiles: str (property)
        The SMILES representation of the mol attribute without explicit hydrogens,
        this attribute (property) is immutable.
    cost: float
        This attribute is used to interact with the fitness functions of :mod:`moldrug.fitness`

    Example
    -------
    .. ipython:: python

        from moldrug import utils, fitness
        import numpy as np
        from copy import copy, deepcopy
        from rdkit import Chem
        i1 = utils.Individual(mol = Chem.MolFromSmiles('CC'), idx = 1, cost = 5)
        i2 = utils.Individual(mol = Chem.MolFromSmiles('CC'), idx = 2, cost = 4)
        i3 = utils.Individual(mol = Chem.MolFromSmiles('CCC'), idx = 3, cost = 4)
        # Show the '==' operation
        print(i1 == i2, i1 == i3)
        # Show that Individual is a hashable object based on the smiles
        print(set([i1,i2,i3]))
        # Show arithmetic operations
        print(i1+i2)
        # How to work with numpy
        array = np.array([i1,i2, i3])
        array_2 = (array*2).astype('float64')
        print(array_2)
        # Show copy
        print(copy(i3), deepcopy(i3))
    """
    def __init__(self, mol: Chem.rdchem.Mol, idx: Union[int, str] = 0, pdbqt: str = None,
                 cost: float = np.inf, randomseed: Union[int, None] = None) -> None:
        """This is the constructor of the class.

        Parameters
        ----------
        mol : Chem.rdchem.Mol, optional
            A valid RDKit molecule.
        idx :Union[int str], optional
            An identification, by default 0
        pdbqt : str, optional
            A valid pdbqt string. If it is not provided it will be generated from mol through utils.confgen
            and the mol attribute will be update with the 3D model, by default None
        cost : float, optional
            This attribute is used to perform operations between Individuals and
            should be used for the cost functions, by default np.inf
        randomseed : Union[None, int], optional
            Provide a seed for the random number generator so that the "same" coordinates can be obtained
            for the attribute pdbqt on multiple runs. If None, the RNG will not be seeded, by default None
        """
        self.mol = mol

        if not pdbqt:
            try:
                self.pdbqt = confgen(self.mol, randomseed=randomseed)
            except Exception:
                self.pdbqt = None
        else:
            self.pdbqt = pdbqt

        self.cost = cost
        self.idx = idx

    @property
    def smiles(self):
        return Chem.MolToSmiles(Chem.RemoveHs(self.mol))

    def __repr__(self):
        return f"{self.__class__.__name__}(idx = {self.idx}, smiles = {self.smiles}, cost = {self.cost})"

    def __hash__(self):
        return hash(self.smiles)

    def __eq__(self, other: object) -> bool:
        if self.__class__ is other.__class__:
            return self.smiles == other.smiles  # self.cost == other.cost and
        else:
            return False  # self.smiles == other

    def __gt__(self, other: object) -> bool:
        return self.cost > other.cost

    def __ge__(self, other: object) -> bool:
        return self.cost >= other.cost

    def __lt__(self, other: object) -> bool:
        return self.cost < other.cost

    def __le__(self, other: object) -> bool:
        return self.cost <= other.cost

    def __neg__(self) -> float:
        return -self.cost

    def __abs__(self) -> float:
        return abs(self.cost)

    def __add__(self, other: object):
        return self.cost + other

    def __radd__(self, other: object):
        return other + self.cost

    def __sub__(self, other: object):
        return self.cost - other

    def __rsub__(self, other: object):
        return other - self.cost

    def __mul__(self, other: object):
        return self.cost * other

    def __rmul__(self, other: object):
        return other * self.cost

    def __truediv__(self, other: object):
        return self.cost / other

    def __rtruediv__(self, other: object):
        return other / self.cost

    def __floordiv__(self, other: object):
        return self.cost // other

    def __rfloordiv__(self, other: object):
        return other // self.cost

    def __mod__(self, other: object):
        return self.cost % other

    def __rmod__(self, other: object):
        return other % self.cost

    def __divmod__(self, other: object):
        return (self.cost // other, self.cost % other)

    def __rdivmod__(self, other: object):
        return (other // self.cost, other % self.cost)

    def __pow__(self, other: object):
        return self.cost ** other

    def __rpow__(self, other: object):
        return other ** self.cost

    def exp(self):
        return np.exp(self.cost)

    def __copy__(self):
        cls = self.__class__
        result = cls.__new__(cls)
        result.__dict__.update(self.__dict__)
        return result

    def __deepcopy__(self, memo):
        cls = self.__class__
        result = cls.__new__(cls)
        memo[id(self)] = result
        for k, v in self.__dict__.items():
            setattr(result, k, deepcopy(v, memo))
        return result


def make_sdf(individuals: List[Individual], sdf_name: str = 'out'):
    """This function create a sdf file from a list of Individuals based on their pdbqt attribute
    This assume that the cost function update the pdbqt attribute after the docking with the conformations obtained
    In the case of multiple receptor the attribute should be a list of valid pdbqt strings.
    Here will export several sdf depending how many pdbqt string are in the pdbqt attribute.

    Parameters
    ----------
    individuals : list[Individual]
        A list of individuals
    sdf_name : str, optional
        The name for the output file. Could be a ``path + sdf_name``.
        The sdf extension will be added by the function, by default 'out'

    Example
    -------
    .. ipython:: python

        import tempfile, os
        from moldrug import utils
        from rdkit import Chem
        # Create some temporal dir
        tmp_path = tempfile.TemporaryDirectory()
        # Creating two individuals
        I1 = utils.Individual(Chem.MolFromSmiles('CCCCl'))
        I2 = utils.Individual(Chem.MolFromSmiles('CCOCCCF'))
        # Creating the pdbqt attribute as a list with the pdbqt attribute (this is just a silly example)
        I1.pdbqt = [I1.pdbqt, I1.pdbqt]
        I2.pdbqt = [I2.pdbqt, I2.pdbqt]
        utils.make_sdf([I1, I2], sdf_name = os.path.join(tmp_path.name, 'out'))
        # Two files were created
        # In the other hand, if the attribute pdbqt is not a list, only one file is going to be created
        # Set pdbqt to the original value
        I1.pdbqt = I1.pdbqt[0]
        I2.pdbqt = I2.pdbqt[0]
        utils.make_sdf([I1, I2], sdf_name = os.path.join(tmp_path.name, 'out'))
        # Only one file will be created if the pdbqt has not len in some of
        # the individuals or they presents different lens as well.
        # In this case the pdbqts will be completely ignored and pdbqt attribute
        # will be used for the construction of the sdf file
        I1.pdbqt = [I1.pdbqt, I1.pdbqt, I1.pdbqt]
        I2.pdbqt = [I2.pdbqt, I2.pdbqt]
        utils.make_sdf([I1, I2], sdf_name = os.path.join(tmp_path.name, 'out'))
    """
    pdbqt_tmp = tempfile.NamedTemporaryFile(suffix='.pdbqt')

    # Check for the attribute pdbqt in all passed individuals and that all of them have the same number of pdbqt
    check = True
    NumbOfpdbqt = set()
    for individual in individuals:
        if isinstance(individual.pdbqt, List):
            NumbOfpdbqt.add(len(individual.pdbqt))
        else:
            check = False
            break
    if len(NumbOfpdbqt) > 1:
        check = False

    if check is True:
        for i in range(list(NumbOfpdbqt)[0]):
            with Chem.SDWriter(f"{sdf_name}_{i+1}.sdf") as w:
                for individual in individuals:
                    with open(pdbqt_tmp.name, 'w') as f:
                        f.write(individual.pdbqt[i])
                    try:
                        pdbqt_mol = PDBQTMolecule.from_file(pdbqt_tmp.name, skip_typing=True)
                        mol = RDKitMolCreate.from_pdbqt_mol(pdbqt_mol)[0]
                        mol.SetProp("_Name",
                                    f"idx :: {individual.idx}, smiles :: {individual.smiles}, "
                                    f"cost :: {individual.cost}")
                        w.write(mol)
                    except Exception:
                        # Should be that the pdbqt is not valid
                        log(
                            f"{individual} does not have a valid pdbqt: {individual.pdbqt}.",
                            LogLevel.ERROR
                        )
            log(f"File {sdf_name}_{i+1}.sdf was created!")
    else:
        with Chem.SDWriter(f"{sdf_name}.sdf") as w:
            for individual in individuals:
                with open(pdbqt_tmp.name, 'w') as f:
                    if len(NumbOfpdbqt) == 0:
                        f.write(individual.pdbqt)
                    else:
                        f.write(individual.pdbqt[0])
                try:
                    pdbqt_mol = PDBQTMolecule.from_file(pdbqt_tmp.name, skip_typing=True)
                    mol = RDKitMolCreate.from_pdbqt_mol(pdbqt_mol)[0]
                    mol.SetProp("_Name",
                                f"idx :: {individual.idx}, smiles :: {individual.smiles}, "
                                f"cost :: {individual.cost}")
                    w.write(mol)
                except Exception:
                    # Should be that the pdbqt is not valid
                    log(
                        f"{individual} does not have a valid pdbqt: {individual.pdbqt}.",
                        LogLevel.ERROR
                    )
        log(f"File {sdf_name}.sdf was createad!")


def _make_kwargs_copy(costfunc, costfunc_kwargs,):
    """Make a copy of the self.costfunc_kwargs.
    It creates a temporal directory.

    Returns
    -------
    dict
        A copy of self.costfunc_kwargs with wd changed if needed
    """
    kwargs_copy = costfunc_kwargs.copy()
    costfunc_jobs_tmp_dir = tempfile.TemporaryDirectory(prefix='.costfunc_moldrug_', dir='.')
    if 'wd' in signature(costfunc).parameters:
        kwargs_copy['wd'] = costfunc_jobs_tmp_dir.name
    return kwargs_copy, costfunc_jobs_tmp_dir


def tar_errors(error_path: str = 'error'):
    """Clean errors in the working directory.
    Convert to error.tar.gz the error_path
    and delete the directory.

    Parameters
    ----------
    error_path : str
        Where the errors are storged.
    """
    if os.path.isdir(error_path):
        if os.listdir(error_path):
            shutil.make_archive('error', 'gztar', error_path)
            
            log(f"\t\t{20*'=+'}")
            log("Check the running warnings and erorrs in error.tar.gz file!", LogLevel.WARNING)
            log(f"\t\t{20*'=+'}")
        shutil.rmtree(error_path)

######################
# Selection functions
######################


def roulette_wheel_selection(p: List[float]):
    """Function to select the offsprings based on their fitness.

    Parameters
    ----------
    p : list[float]
        Probabilities

    Returns
    -------
    int
        The selected index
    """
    c = np.cumsum(p)
    r = sum(p) * random.random() # noqa
    ind = np.argwhere(r <= c)
    return ind[0][0]


def to_dataframe(individuals: List[Individual], return_mol: bool = False) -> pd.DataFrame:
    """Convert a list of individuals to a DataFrame

    Parameters
    ----------
    individuals : List[Individual]
        The list of individuals
    return_mol : bool, optional
        If True the attribute mol will bot be return, by default False

    Returns
    -------
    pd.DataFrame
        The DataFrame
    """
    list_of_dictionaries = []
    for individual in individuals:
        dictionary = individual.__dict__.copy()
        if not return_mol:
            del dictionary['mol']
        list_of_dictionaries.append(dictionary)
    return pd.DataFrame(list_of_dictionaries)

################################################
# Here are some important functions to work with
################################################


def run(command: str, shell: bool = True, executable: str = '/bin/bash'):
    """This function is just a useful wrapper around subprocess.run

    Parameters
    ----------
    command : str
        Any command to execute.
    shell : bool, optional
        keyword of ``subprocess.Popen`` and ``subprocess.Popen``, by default True
    executable : str, optional
        keyword of ``subprocess.Popen`` and ``subprocess.Popen``, by default '/bin/bash'

    Returns
    -------
    object
        The processes returned by Run.

    Raises
    ------
    RuntimeError
        In case of non-zero exit status on the provided command.
    """

    process = subprocess.run(command, shell=shell, executable=executable,
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    returncode = process.returncode
    if returncode != 0:
        # print(f'Command {command} returned non-zero exit status {returncode}')
        raise RuntimeError(process.stderr)
    return process


def confgen(mol: Chem.rdchem.Mol, return_mol: bool = False, randomseed: Union[int, None] = None):
    """Create a 3D model from a smiles and return a pdbqt string and, a mol if ``return_mol = True``.

    Parameters
    ----------
    mol : Chem.rdchem.Mol
        A valid RDKit molecule.
    return_mol : bool, optional
        If true the function will also return the ``rdkit.Chem.rdchem.Mol``, by default False
    randomseed : Union[None, int], optional
       Provide a seed for the random number generator so that the same coordinates
       can be obtained for a molecule on multiple runs. If None, the RNG will not be seeded, by default None

    Returns
    -------
    tuple or str
        If ``return_mol = True`` it will return a tuple ``(str[pdbqt], Chem.rdchem.Mol)``,
        if not only a ``str`` that represents the pdbqt.
    """
    mol = Chem.AddHs(mol)
    if randomseed is None:
        randomSeed = -1
    else:
        randomSeed = randomseed

    AllChem.EmbedMolecule(mol, randomSeed=randomSeed)
    # The optimization introduce some sort of non-reproducible results.
    # For that reason is not used when randomseed is set
    if not randomseed:
        AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
    preparator = MoleculePreparation()
    mol_setups = preparator.prepare(mol)
    pdbqt_string = PDBQTWriterLegacy.write_string(mol_setups[0])[0]
    if return_mol:
        return (pdbqt_string, mol)
    else:
        return pdbqt_string


def update_reactant_zone(parent: Chem.rdchem.Mol, offspring: Chem.rdchem.Mol,
                         parent_replace_ids: List[int] = None, parent_protected_ids: List[int] = None):
    """This function will find the difference between offspring and parent
    based on the Maximum Common Substructure (MCS).
    This difference will be consider offspring_replace_ids.
    Because after a reaction the indexes of the product could change respect to the reactant,
    the parent_replace_ids could change.
    The function will map the index of the parent to the offspring based on MCS. If on those indexes some of the
    parent_replace_ids are still present, they will be updated based on the offspring and also added
    to offspring_replace_ids. Similarly will be done for the parent_protected_ids.

    Parameters
    ----------
    parent : Chem.rdchem.Mol
        The original molecule from where offspring was generated
    offspring : Chem.rdchem.Mol
        A derivative of parent
    parent_replace_ids : List[int], optional
        A list of replaceable indexes in the parent, by default None
    parent_protected_ids : List[int], optional
        A list of protected indexes in the parent, by default None
    Returns
    -------
    tuple[list[int]]
        The function returns a tuple composed by two list of integers.
        The first list is offspring_replace_ids and the second one offspring_protected_ids.
    """

    # Finding Maximum Common Substructure (MCS) and getting the SMARTS
    mcs = rdFMCS.FindMCS([parent, offspring], matchValences=True, ringMatchesRingOnly=True)
    mcs_mol = Chem.MolFromSmarts(mcs.smartsString)

    # Get the index corresponding to the MCS for both parent and offspring
    # The ordering of the indices corresponds to the atom ordering in GetSubstructMatch.
    # Therefore we have a mapping between the two molecules
    match_parent = parent.GetSubstructMatch(mcs_mol)
    match_offspring = offspring.GetSubstructMatch(mcs_mol)
    mapping = dict(zip(match_parent, match_offspring))

    offspring_replace_ids = []
    for atom in offspring.GetAtoms():
        if atom.GetIdx() not in match_offspring:
            offspring_replace_ids.append(atom.GetIdx())

    # Because after a reaction the index could change we must update the original index.
    # If the try raise exception is because the parent_replace_ids
    # are not present in the MCS what means that they were already submitted to a reaction
    # Therefore we dont need to update them neither add to the offspring_replace_ids
    if parent_replace_ids:
        for idx in parent_replace_ids:
            try:
                offspring_replace_ids.append(mapping[idx])
            except KeyError:
                pass

    offspring_protected_ids = []
    if parent_protected_ids:
        for idx in parent_protected_ids:
            try:
                offspring_protected_ids.append(mapping[idx])
            except KeyError:
                pass
    return offspring_replace_ids, offspring_protected_ids


def get_sim(ms: List[Chem.rdchem.Mol], ref_fps: List):
    """Get the molecules with higher similarity to each member of ref_fps.

    Parameters
    ----------
    ms : list[Chem.rdchem.Mol]
        List of molecules
    ref_fps : list[AllChem.GetMorganFingerprintAsBitVect(mol, 2)]
        A list of reference fingerprints

    Returns
    -------
    list[Chem.rdchem.Mol]
        A list of molecules with the higher similarity with their corresponded ref_fps value.
    """
    output = []
    fps1 = [AllChem.GetMorganFingerprintAsBitVect(m, 2) for m in ms]
    for fp in fps1:
        v = DataStructs.BulkTanimotoSimilarity(fp, ref_fps)
        i = np.argmax(v)
        output.append([v[i], i])
    return output


def get_similar_mols(mols: List, ref_mol: Chem.rdchem.Mol, pick: int, beta: float = 0.01):
    """Pick the similar molecules from mols respect to ref_mol using a roulette wheel selection strategy.

    Parameters
    ----------
    mols : list
        The list of molecules from where to pick molecules.
    ref_mol : Chem.rdchem.Mol
        The reference molecule
    pick : int
        Number of molecules to pick from mols
    beta : float, optional
        Selection threshold, by default 0.01

    Returns
    -------
    list
        A list of picked molecules.
    """
    if pick >= len(mols):
        return mols
    else:
        ref_fp = AllChem.GetMorganFingerprintAsBitVect(ref_mol, 2)
        fps = [AllChem.GetMorganFingerprintAsBitVect(mol, 2) for mol in mols]
        similarities = np.array(DataStructs.BulkTanimotoSimilarity(ref_fp, fps))
        probs = np.exp(beta * similarities) / np.sum(np.exp(beta * similarities))
        cumsum = np.cumsum(probs)
        indexes = set()
        while len(indexes) != pick:
            r = sum(probs) * random.random() # noqa
            indexes.add(np.argwhere(r <= cumsum)[0][0])
        return [mols[index] for index in indexes]


def lipinski_filter(mol: Chem.rdchem.Mol, maxviolation: int = 2):
    """Implementation of Lipinski filter.

    Parameters
    ----------
    mol : Chem.rdchem.Mol
        An RDKit molecule.
    maxviolation : int, optional
        Maximum number of violations. Above this value the function return False, by default 2

    Returns
    -------
    bool
        True if the molecule present less than maxviolation violations; otherwise False.
    """
    filter_prop = {
        'NumHAcceptors': {'method': Lipinski.NumHAcceptors, 'cutoff': 10},
        'NumHDonors': {'method': Lipinski.NumHDonors, 'cutoff': 5},
        'wt': {'method': Descriptors.MolWt, 'cutoff': 500},
        'MLogP': {'method': Descriptors.MolLogP, 'cutoff': 5},
        'NumRotatableBonds': {'method': Lipinski.NumRotatableBonds, 'cutoff': 10},
        'TPSA': {'method': Chem.MolSurf.TPSA, 'cutoff': 140},
    }
    cont = 0
    for prop in filter_prop:

        if filter_prop[prop]['method'](mol) > filter_prop[prop]['cutoff']:
            cont += 1
        if cont >= maxviolation:
            return False
    return True


def lipinski_profile(mol: Chem.rdchem.Mol):
    """See: https://www.rdkit.org/docs/source/rdkit.Chem.Lipinski.html?highlight=lipinski#module-rdkit.Chem.Lipinski

    Parameters
    ----------
    mol : Chem.rdchem.Mol
        An RDKit molecule.

    Returns
    -------
    dict
        A dictionary with molecular properties.
    """
    properties = {
        'NumHAcceptors': {'method': Lipinski.NumHAcceptors, 'cutoff': 10},
        'NumHDonors': {'method': Lipinski.NumHDonors, 'cutoff': 5},
        'wt': {'method': Descriptors.MolWt, 'cutoff': 500},
        'MLogP': {'method': Descriptors.MolLogP, 'cutoff': 5},
        'NumRotatableBonds': {'method': Lipinski.NumRotatableBonds, 'cutoff': 10},
        'TPSA': {'method': Chem.MolSurf.TPSA, 'cutoff': range(0, 140)},

        'FractionCSP3': {'method': Lipinski.FractionCSP3, 'cutoff': None},
        'HeavyAtomCount': {'method': Lipinski.HeavyAtomCount, 'cutoff': None},
        'NHOHCount': {'method': Lipinski.NHOHCount, 'cutoff': None},
        'NOCount': {'method': Lipinski.NOCount, 'cutoff': None},
        'NumAliphaticCarbocycles': {'method': Lipinski.NumAliphaticCarbocycles, 'cutoff': None},
        'NumAliphaticHeterocycles': {'method': Lipinski.NumAliphaticHeterocycles, 'cutoff': None},
        'NumAliphaticRings': {'method': Lipinski.NumAliphaticRings, 'cutoff': None},
        'NumAromaticCarbocycles': {'method': Lipinski.NumAromaticCarbocycles, 'cutoff': None},
        'NumAromaticHeterocycles': {'method': Lipinski.NumAromaticHeterocycles, 'cutoff': None},
        'NumAromaticRings': {'method': Lipinski.NumAromaticRings, 'cutoff': None},
        'NumHeteroatoms': {'method': Lipinski.NumHeteroatoms, 'cutoff': None},
        'NumSaturatedCarbocycles': {'method': Lipinski.NumSaturatedCarbocycles, 'cutoff': None},
        'NumSaturatedHeterocycles': {'method': Lipinski.NumSaturatedHeterocycles, 'cutoff': None},
        'NumSaturatedRings': {'method': Lipinski.NumSaturatedRings, 'cutoff': None},
        'RingCount': {'method': Lipinski.RingCount, 'cutoff': None},
    }
    profile = {}
    for prop in properties:
        profile[prop] = properties[prop]['method'](mol)
        # print(f"{prop}: {properties[prop]['method'](mol)}. cutoff: {properties[prop]['cutoff']}")
    return profile


def LargerTheBest(Value: float, LowerLimit: float, Target: float, r: float = 1) -> float:
    """Desirability function used when larger values are the targets. If Value is higher
    or equal than the target it will return 1; if it is lower than LowerLimit it will return 0;
    else a number between 0 and 1.
    You can also check:
    doi:10.1016/j.chemolab.2011.04.004
    https://www.youtube.com/watch?v=quz4NW0uIYw&list=PL6ebkIZFT4xXiVdpOeKR4o_sKLSY0aQf_&index=3

    Parameters
    ----------
    Value : float
        Value to test.
    LowerLimit : float
        Lower value accepted. Lower than this one will return 0.
    Target : float
        The target value. On this value (or higher) the function takes 1 as value.
    r : float, optional
        This is the exponent of the interpolation. Could be used to control the interpolation, by default 1

    Returns
    -------
    float
        A number between 0 and 1. Been 1 the desireable value to get.
    """
    if Value < LowerLimit:
        return 0.0
    elif Value <= Target:
        return ((Value - LowerLimit) / (Target - LowerLimit))**r
    else:
        return 1.0


def SmallerTheBest(Value: float, Target: float, UpperLimit: float, r: float = 1) -> float:
    """Desirability function used when lower values are the targets. If Value is lower or
    equal than the target it will return 1; if it is higher than UpperLimit it will return 0;
    else a number between 0 and 1.

    Parameters
    ----------
    Value : float
        Value to test.
    Target : float
        The target value. On this value (or lower) the function takes 1 as value.
    UpperLimit : float
        Upper value accepted. Higher than this one will return 0.
    r : float, optional
        This is the exponent of the interpolation. Could be used to control the interpolation, by default 1

    Returns
    -------
    float
        A number between 0 and 1. Been 1 the desireable value to get.
    """
    if Value < Target:
        return 1.0
    elif Value <= UpperLimit:
        return ((UpperLimit - Value) / (UpperLimit - Target))**r
    else:
        return 0.0


def NominalTheBest(Value: float, LowerLimit: float, Target: float,
                   UpperLimit: float, r1: float = 1, r2: float = 1) -> float:
    """Desirability function used when a target value is desired. If Value is
    lower or equal than the LowerLimit it will return 0; as well values higher or equal than
    UpperLimit; else a number between 0 and 1.

    Parameters
    ----------
    Value : float
        Value to test.
    LowerLimit : float
        Lower value accepted. Lower than this one will return 0.
    Target : float
        The target value. On this value the function takes 1 as value.
    UpperLimit : float
        Upper value accepted. Higher than this one will return 0.
    r1 : float, optional
        This is the exponent of the interpolation from LowerLimit to Target.
        Could be used to control the interpolation, by default 1
    r2 : float, optional
        This is the exponent of the interpolation from Target to UpperLimit.
        Could be used to control the interpolation, by default 1

    Returns
    -------
    float
        A number between 0 and 1. Been 1 the desireable value to get.
    """
    if Value < LowerLimit:
        return 0.0
    elif Value <= Target:
        return ((Value - LowerLimit) / (Target - LowerLimit))**r1
    elif Value <= UpperLimit:
        return ((UpperLimit - Value) / (UpperLimit - Target))**r2
    else:
        return 0.0


def DerringerSuichDesirability():
    """A warper around the implemented desirability functions

    Returns
    -------
    dict
        A dict with key name of the desirability and value the corresponded function
    """
    my_dict = {
        'LargerTheBest': LargerTheBest,
        'SmallerTheBest': SmallerTheBest,
        'NominalTheBest': NominalTheBest
    }
    return my_dict


# Saving data
def full_pickle(title: str, data: object):
    """Normal pickle.

    Parameters
    ----------
    title : str
        Name of the file without extension, .pkl will be added by default.
    data : object
        Any serializable python object
    """
    with open(f'{title}.pkl', 'wb') as pkl:
        pickle.dump(data, pkl)


def loosen(file: str):
    """Unpickle a pickled object.

    Parameters
    ----------
    file : str
        The path to the file who store the pickle object.

    Returns
    -------
    object
        The python object.
    """
    with open(file, 'rb') as pkl:
        data = pickle.load(pkl)
    return data


def compressed_pickle(title: str, data: object):
    """Compress Python object. First cPickle it and then bz2.BZ2File compressed it.

    Parameters
    ----------
    title : str
         Name of the file without extensions, .pbz2 will be added by default
    data : object
        Any serializable python object
    """
    with bz2.BZ2File(f'{title}.pbz2', 'w') as f:
        pickle.dump(data, f)


def decompress_pickle(file: str):
    """Decompress CPickle objects compressed first with bz2 formats

    Parameters
    ----------
    file : str
         This is the cPickle files compressed with bz2.BZ2File. (as a convention with extension .pbz2, but not needed)

    Returns
    -------
    object
        The python object.
    """
    data = bz2.BZ2File(file, 'rb')
    data = pickle.load(data)
    return data


def is_iter(obj):
    """Check if obj is iterable

    Parameters
    ----------
    obj : Any
        Any python object

    Returns
    -------
    bool
        Tru if obj iterable, False if not
    """
    try:
        for _ in obj:
            return True
    except TypeError:
        return False


def import_sascorer():
    """Function to import sascorer from RDConfig.RDContribDir of RDKit

    Returns
    -------
    module
        The sascorer module ready to use.
    """
    # In order to import sascorer from RDConfig.RDContribDir
    import importlib.util as importlib_util

    from rdkit.Chem import RDConfig
    spec = importlib_util.spec_from_file_location(
        'sascorer', os.path.join(RDConfig.RDContribDir, 'SA_Score', 'sascorer.py'))
    sascorer = importlib_util.module_from_spec(spec)
    spec.loader.exec_module(sascorer)
    return sascorer


def deep_update(target_dict: dict, update_dict: dict) -> dict:
    """Recursively update a dictionary with the key-value pairs from another dictionary.
    Inpired on https://stackoverflow.com/questions/3232943/update-value-of-a-nested-dictionary-of-varying-depth

    Parameters
    ----------
    target_dict : dict
        The dictionary to be updated
    update_dict : dict
        The dictionary providing the updates

    Example
    -------
    .. ipython:: python

        from moldrug.utils import deep_update
        target = {'a': 1, 'b': {'c': 2, 'd': 3}}
        updates = {'b': {'c': 4, 'e': 5}, 'f': 6}
        result = deep_update(target, updates)
        print(result)
        # Output: {'a': 1, 'b': {'c': 4, 'd': 3, 'e': 5}, 'f': 6}

    Returns
    -------
    dict
        The updated dictionary
    """
    for key, value in update_dict.items():
        if isinstance(value, collections.abc.Mapping):
            # Recursive update for nested dictionaries
            target_dict[key] = deep_update(target_dict.get(key, {}), value)
        else:
            # Update the value if it's not a dictionary
            target_dict[key] = value
    return target_dict


def softmax(x, axis=None):
    x_max = np.amax(x, axis=axis, keepdims=True)
    exp_x_shifted = np.exp(x - x_max)
    return exp_x_shifted / np.sum(exp_x_shifted, axis=axis, keepdims=True)
###########################################################################
#   Classes to work with Vina
###########################################################################


class Atom:
    """This is a simple class to wrap a pdbqt Atom.
    It is based on https://userguide.mdanalysis.org/stable/formats/reference/pdbqt.html#writing-out.
    """
    def __init__(self, line):
        self.lineType = "ATOM"
        self.serial = int(line[6:11])
        self.name = line[12:16].strip()
        self.altLoc = line[16]
        self.resName = line[17:21].strip()
        self.chainID = line[21]
        self.resSeq = int(line[22:26])
        self.iCode = line[26]
        self.x = float(line[30:38])
        self.y = float(line[38:46])
        self.z = float(line[46:54])
        self.occupancy = line[54:60].strip()
        self.tempFactor = line[60:66].strip()
        self.partialChrg = line[66:76].strip()
        self.atomType = line[78:80].strip()

    def __getitem__(self, key):
        return self.__dict__[key]


class CHUNK_VINA_OUT:
    """This class will be used by VINA_OUT in order to read the pdbqt ouput of a vina docking results.
    """
    def __init__(self, chunk):
        self.chunk = chunk
        self.atoms = []
        self.run = None
        self.freeEnergy = None
        self.RMSD1 = None
        self.RMSD2 = None
        self.parse()

    def parse(self):
        for line in self.chunk:
            if line.startswith("MODEL"):
                self.run = int(line[5:])
            elif line.startswith("REMARK VINA RESULT:"):
                (self.freeEnergy, self.RMSD1, self.RMSD2) = [float(number) for number in line.split(":")[-1].split()]
            elif line.startswith("ATOM"):
                self.atoms.append(Atom(line))

    def get_atoms(self):
        """Return a list of all atoms.

        If to_dict is True, each atom is represented as a dictionary.
        Otherwise, a list of Atom objects is returned."""
        return [x.__dict__ for x in self.atoms]

    def write(self, name=None):
        if name:
            with open(name, "w") as f:
                f.writelines(self.chunk)
        else:
            with open(f"Run_{self.run}.pdbqt", "w") as f:
                f.writelines(self.chunk)


class VINA_OUT:
    """
    Vina class to handle vina output. Think about use meeko in the future!
    """
    def __init__(self, file):
        self.file = file

        self.chunks = []
        self.parse()

    def parse(self):
        with open(self.file, "r") as input_file:
            lines = input_file.readlines()
        i = 0
        while i < len(lines):

            if lines[i].startswith("MODEL"):
                j = i
                tmp_chunk = []
                while (not lines[j].startswith("ENDMDL")) and (j < len(lines)):
                    tmp_chunk.append(lines[j])
                    j += 1
                    i += 1
                tmp_chunk.append("ENDMDL\n")
                self.chunks.append(CHUNK_VINA_OUT(tmp_chunk))

            i += 1

    def BestEnergy(self, write=False):
        min_chunk = min(self.chunks, key=lambda x: x.freeEnergy)
        if write:
            min_chunk.write("best_energy.pdbqt")
        return min_chunk



if __name__ == '__main__':
    pass
