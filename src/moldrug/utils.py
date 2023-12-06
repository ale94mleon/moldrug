#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import bz2
import collections.abc
import datetime
import multiprocessing as mp
import os
import random
import shutil
import subprocess
import tempfile
import time
from copy import deepcopy
from inspect import signature
from typing import Dict, Iterable, List, Union
from warnings import warn

import dill as pickle
import numpy as np
import pandas as pd
import tqdm
from crem.crem import grow_mol, mutate_mol
from meeko import (MoleculePreparation, PDBQTMolecule, PDBQTWriterLegacy,
                   RDKitMolCreate)
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, DataStructs, Descriptors, Lipinski, rdFMCS

from moldrug import __version__

RDLogger.DisableLog('rdApp.*')
# # in order to pickle the isotope properties of the molecule
# Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.AllProps)

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

#################################
# Classes to work with moldrug
#################################


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
                        print(f"{individual} does not have a valid pdbqt: {individual.pdbqt}.")
            print(f" File {sdf_name}_{i+1}.sdf was created!")
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
                    print(f"{individual} does not have a valid pdbqt: {individual.pdbqt}.")
        print(f"File {sdf_name}.sdf was createad!")


def _make_kwargs_copy(costfunc, costfunc_kwargs,):
    """Make a copy of the self.costfunc_kwargs.
    It creates a temporal directory.

    Returns
    -------
    dict
        A copy of self.costfunc_kwargs with wd changed if needed
    """
    kwargs_copy = costfunc_kwargs.copy()
    costfunc_jobs_tmp_dir = tempfile.TemporaryDirectory(prefix='.costfunc_MolDrug_', dir='.')
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
            print(f"\n{50*'=+'}")
            print("Note: Check the running warnings and erorrs in error.tar.gz file!")
            print(f"{50*'=+'}\n")
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


class Local:
    """This class is used to genereate close solutions to the seed molecule.
    It use :meth:`crem.crem.grow_mol`.

    Attributes
    ----------
    randomseed : Union[None, int]
        The random seed to use with random module.
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
    AddHs : bool
        In case explicit hydrogens should be added.
    pop : list[:meth:`moldrug.utils.Individuals`]
        The final population sorted by cost.
    """
    def __init__(self, seed_mol: Chem.rdchem.Mol, crem_db_path: str, costfunc: object, grow_crem_kwargs: Dict = None,
                 costfunc_kwargs: Dict = None, AddHs: bool = False, randomseed: Union[None, int] = None,
                 deffnm: str = 'local') -> None:
        """Creator

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

        Raises
        ------
        Exception
            In case that some problem occured during the creation of the Individula from the seed_mol
        ValueError
            In case of incorrect definition of grow_crem_kwargs and/or costfunc_kwargs.
            They must be None or a dict instance.
        """
        self.randomseed = randomseed
        if self.randomseed is not None:
            random.seed(randomseed)

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
        if os.path.exists(crem_db_path):
            self.crem_db_path = os.path.abspath(crem_db_path)
        else:
            raise FileNotFoundError(f"{crem_db_path = } does not exists or is not accesible")

        self.grow_crem_kwargs = grow_crem_kwargs
        self.costfunc = costfunc
        self.costfunc_kwargs = costfunc_kwargs
        self.pop = [self.InitIndividual]

    def __call__(self, njobs: int = 1, pick: int = None):
        """Call deffinition

        Parameters
        ----------
        njobs : int, optional
            The number of jobs for parallelization, the module multiprocessing will be used, by default 1
        pick : int, optional
            How many molecules take from the generated throgh the grow_mol CReM operation,
            by default None which means all generated.
        """
        # Check version of MolDrug
        if self.__moldrug_version != __version__:
            warn(f"{self.__class__.__name__} was initilized with moldrug-{self.__moldrug_version} "
                 f"but was called with moldrug-{__version__}")
        self.grow_crem_kwargs.update({'return_mol': True})
        new_mols = list(grow_mol(self._seed_mol, self.crem_db_path, **self.grow_crem_kwargs))
        if pick:
            random.shuffle(new_mols)
            new_mols = new_mols[:pick]
            new_mols = [item[1] for item in new_mols]

        idx0 = len(self.pop)
        for i, mol in enumerate(new_mols):
            individual = Individual(mol, idx=idx0 + i, randomseed=self.randomseed)
            if individual.pdbqt:
                self.pop.append(individual)

        # Calculating cost of each individual
        # Creating the arguments
        args_list = []
        # Make a copy of the self.costfunc_kwargs
        kwargs_copy, costfunc_jobs_tmp_dir = _make_kwargs_copy(self.costfunc, self.costfunc_kwargs)

        for individual in self.pop:
            args_list.append((individual, kwargs_copy))

        print('Calculating cost function...')
        pool = mp.Pool(njobs)
        self.pop = [individual for individual in tqdm.tqdm(pool.imap(self.__costfunc__, args_list), total=len(args_list))]
        pool.close()

        # Clean directory
        costfunc_jobs_tmp_dir.cleanup()
        # Tar errors
        tar_errors('error')

        # Printing how long was the simulation
        print(f"Finished at {datetime.datetime.now().strftime('%c')}.\n")

    def __costfunc__(self, args_list):
        Individual, kwargs = args_list
        # This is just to use the progress bar on pool.imap
        return self.costfunc(Individual, **kwargs)

    def pickle(self, title: str, compress: bool = False):
        """Method to pickle the whole Local class

        Parameters
        ----------
        title : str
            Name of the object which will be compleated with the correposnding
            extension depending if compress is set to True or False.
        compress : bool, optional
            Use compression, by default False. If True :meth:`moldrug.utils.compressed_pickle` will be used;
            if not :meth:`moldrug.utils.full_pickle` will be used instead.
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
        The random seed to use with random module.
    __moldrug_version__ : str
        The molDrug version.
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
        Bias the search upon similar molecules. If True :meth:`modrug.utils.get_similar_mols` is used after the
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
    SawIndividuals : set[:meth:`moldrug.utils.Individuals`]
        All the Individulas saw during the optimizations.
    acceptance : dict
        A dictionary with key the Generation id and as value another dictionary
        with keys ``accepeted`` and ``generated`` with the number of accepted and genereated
        individuals on the generation respectively.
    AddHs : bool
        In case explicit hydrogens should be added.
    _seed_mol : list[Chem.rdchem.Mol]
        The list of seed molecules.
    InitIndividual : :meth:`moldrug.utils.Individuals`
        The initial individual based on _seed_mol.
    pop : list[:meth:`moldrug.utils.Individuals`]
        The final population sorted by cost.
    best_cost : list[float]
        The list of best cost for each generations.
    avg_cost : list[float]
        The list of average cost for each generations.
    """
    def __init__(self, seed_mol: Union[Chem.rdchem.Mol, Iterable[Chem.rdchem.Mol]],
                 costfunc: object, costfunc_kwargs: Dict, crem_db_path: str, maxiter: int = 10, popsize: int = 20,
                 beta: float = 0.001, pc: float = 1, get_similar: bool = False, mutate_crem_kwargs: Union[None, Dict] = None,
                 save_pop_every_gen: int = 0, checkpoint: bool = False, deffnm: str = 'ga',
                 AddHs: bool = False, randomseed: Union[None, int] = None) -> None:
        """Constructor

        Parameters
        ----------
        seed_mol : Union[Chem.rdchem.Mol, Iterable[Chem.rdchem.Mol]]
            The seed molecule submitted to genetic algorithm optimization on the chemical space. Could be only one RDKit
            molecule or more than one specified in an Iterable object.
        costfunc : object
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
            Selection pressure, by default 0.001
        pc : float, optional
            Proportion of children, by default 1
        get_similar : bool, optional
            If True the searching will be bias to similar molecules, by default False
        mutate_crem_kwargs : Union[None, Dict], optional
            Parameters for mutate_mol of CReM, by default {}
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
        if self.randomseed is not None:
            random.seed(randomseed)

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
        self.mutate_crem_kwargs = {
            'radius': 3,
            'min_size': 1,
            'max_size': 8,
            'min_inc': -5,
            'max_inc': 3,
            'ncores': 1,
        }
        self.mutate_crem_kwargs.update(mutate_crem_kwargs)

        # Saving parameters
        self.save_pop_every_gen = save_pop_every_gen
        self.checkpoint = checkpoint
        self.deffnm = deffnm

        # Tracking parameters
        self.NumCalls = 0
        self.NumGens = 0
        self.SawIndividuals = set()
        self.acceptance = dict()

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
        #     _ = [atom.SetIntProp('label_MolDrug', atom.GetIdx()) for atom in seed_mol.GetAtoms()]

        # Create the first Individual
        self.InitIndividual = Individual(self._seed_mol[0], idx=0, randomseed=self.randomseed)
        self.pop = []

    def __call__(self, njobs: int = 1):
        """Call definition

        Parameters
        ----------
        njobs : int, optional
            The number of jobs for parallelization, the module multiprocessing will be used, by default 1,

        Raises
        ------
        RuntimeError
            Error during the initialization of the population.
        """
        ts = time.time()
        # Counting the calls
        self.NumCalls += 1

        # Check version of MolDrug
        if self.__moldrug_version__ != __version__:
            warn(f"{self.__class__.__name__} was initialized with moldrug-{self.__moldrug_version__} "
                 f"but was called with moldrug-{__version__}")

        # Here we will update if needed some parameters for
        # the crem operations that could change between different calls.
        # We need to return the molecule, so we override the possible user definition respect to this keyword
        self.mutate_crem_kwargs['return_mol'] = True

        # Initialize Population
        # In case that the populating exist there is not need to initialize.
        if len(self.pop) == 0:
            GenInitStructs = []
            # in case that the input has the popsize memebers there is not need to generate new structures
            if len(self._seed_mol) < self.popsize:
                for mol in self._seed_mol:
                    tmp_GenInitStructs = list(mutate_mol(mol, self.crem_db_path, **self.mutate_crem_kwargs))
                    tmp_GenInitStructs = [mol for (_, mol) in tmp_GenInitStructs]
                    GenInitStructs += tmp_GenInitStructs
                # Checking for possible scenarios
                if len(GenInitStructs) == 0:
                    raise RuntimeError("Something really strange happened. The seed_mol did not "
                                       "generate any new molecule during the initialization of the population. "
                                       "Check the provided crem parameters!")
                if len(GenInitStructs) < (self.popsize - len(self._seed_mol)):
                    print('The initial population has repeated elements')
                    # temporal solution
                    GenInitStructs += random.choices(GenInitStructs,
                                                     k=self.popsize - len(GenInitStructs) - len(self._seed_mol))
                elif len(GenInitStructs) > (self.popsize - 1):
                    # Selected random sample from the generation
                    GenInitStructs = random.sample(GenInitStructs, k=self.popsize - len(self._seed_mol))
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

            # Calculating cost of each individual
            # Creating the arguments
            args_list = []
            # Make a copy of the self.costfunc_kwargs
            # Make a copy of the self.costfunc_kwargs
            kwargs_copy, costfunc_jobs_tmp_dir = _make_kwargs_copy(self.costfunc, self.costfunc_kwargs)

            for individual in self.pop:
                args_list.append((individual, kwargs_copy))

            print(f'\n\nCreating the first population with {len(self.pop)} members:')
            try:
                pool = mp.Pool(njobs)
                self.pop = [individual for individual in tqdm.tqdm(pool.imap(self.__costfunc__, args_list), total=len(args_list))]
                pool.close()
            except Exception as e1:
                warn("Parallelization did not work. Trying with serial...")
                try:
                    self.pop = [self.__costfunc__(args) for args in tqdm.tqdm(args_list, total=len(args_list))]
                except Exception as e2:
                    raise RuntimeError("Serial did not work either. Here are the ucurred exceptions:\n"
                                       f"=========Parellel=========:\n {e1}\n"
                                       f"==========Serial==========:\n {e2}")
            # Clean directory
            costfunc_jobs_tmp_dir.cleanup()

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
            print(f"Initial Population: Best Individual: {self.pop[0]}")
            print(f"Accepted rate: {self.acceptance[self.NumGens]['accepted']} / {self.acceptance[self.NumGens]['generated']}\n")
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

        # Saving tracking variables, the first population, outside the if to take into account second calls
        # with different population provided by the user.
        self.SawIndividuals.update(self.pop)

        # Saving population in disk if it was required
        if self.save_pop_every_gen:
            compressed_pickle(f"{self.deffnm}_pop", (self.NumGens, sorted(self.pop)))
            make_sdf(sorted(self.pop), sdf_name=f"{self.deffnm}_pop")
            if self.checkpoint:
                compressed_pickle('cpt', self)

        # Main Loop
        # Another control variable. In case that the __call__ method is used more than ones.
        number_of_previous_generations = len(self.best_cost)
        for it in range(self.maxiter):
            # Saving Number of Generations
            self.NumGens += 1

            # Probabilities Selections
            probs = softmax((-self.beta * np.array(self.pop)).astype('float64'))
            if any(np.isnan(probs)):
                probs = np.nan_to_num(probs)

            popc = []
            for _ in range(self.nc):
                # Perform Roulette Wheel Selection
                parent = self.pop[roulette_wheel_selection(probs)]

                # Perform Mutation (this mutation is some kind of crossover but with CReM library)
                children = self.mutate(parent)

                # Save offspring population
                # I will save only those offsprings that were not seen and that have a correct pdbqt file
                if children not in self.SawIndividuals and children.pdbqt:
                    children.genID = self.NumGens
                    children.kept_gens = set()
                    popc.append(children)

            if popc:  # Only if there are new members
                # Calculating cost of each offspring individual (Doing Docking)

                # Creating the arguments
                args_list = []
                # Make a copy of the self.costfunc_kwargs
                kwargs_copy, costfunc_jobs_tmp_dir = _make_kwargs_copy(self.costfunc, self.costfunc_kwargs)

                NumbOfSawIndividuals = len(self.SawIndividuals)
                for (i, individual) in enumerate(popc):
                    # Add idx label to each individual
                    individual.idx = i + NumbOfSawIndividuals
                    # The problem here is that we are not being general for other possible Cost functions.
                    args_list.append((individual, kwargs_copy))
                print(f'Evaluating generation {self.NumGens} / {self.maxiter + number_of_previous_generations}:')

                # Calculating cost fucntion in parallel
                try:
                    pool = mp.Pool(njobs)
                    popc = [individual for individual in tqdm.tqdm(pool.imap(self.__costfunc__, args_list), total=len(args_list))]
                    pool.close()
                except Exception as e1:
                    warn("Parallelization did not work. Trying with serial...")
                    try:
                        popc = [self.__costfunc__(args) for args in tqdm.tqdm(args_list, total=len(args_list))]
                    except Exception as e2:
                        raise RuntimeError("Serial did not work either. Here are the ucurred exceptions:\n"
                                           f"=========Parellel=========:\n {e1}\n"
                                           f"==========Serial==========:\n {e2}")

                # Clean directory
                costfunc_jobs_tmp_dir.cleanup()

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

            # Saving population in disk if it was required
            if self.save_pop_every_gen:
                # Save every save_pop_every_gen and always the last population
                if self.NumGens % self.save_pop_every_gen == 0 or it + 1 == self.maxiter:
                    compressed_pickle(f"{self.deffnm}_pop", (self.NumGens, self.pop))
                    make_sdf(self.pop, sdf_name=f"{self.deffnm}_pop")
                    if self.checkpoint:
                        compressed_pickle('cpt', self)

            # Show Iteration Information
            print(f"Generation {self.NumGens}: Best Individual: {self.pop[0]}.")
            print(f"Accepted rate: {self.acceptance[self.NumGens]['accepted']} / {self.acceptance[self.NumGens]['generated']}\n")

        # Printing summary information
        print(f"\n{50*'=+'}\n")
        print(f"The simulation finished successfully after {self.NumGens} generations with"
              f"a population of {self.popsize} individuals. "
              f"A total number of {len(self.SawIndividuals)} Individuals were seen during the simulation.")
        print(f"Initial Individual: {self.InitIndividual}")
        print(f"Final Individual: {self.pop[0]}")
        print(f"The cost function dropped in {self.InitIndividual - self.pop[0]} units.")
        print(f"\n{50*'=+'}\n")

        # Tar errors
        tar_errors('error')

        # Printing how long was the simulation
        print(f"Total time ({self.maxiter} generations): {time.time() - ts:>5.2f} (s).\n"
              f"Finished at {datetime.datetime.now().strftime('%c')}.\n")

    def __costfunc__(self, args_list):
        individual, kwargs = args_list
        # This is just to use the progress bar on pool.imap
        return self.costfunc(individual, **kwargs)

    def mutate(self, individual: Individual):
        """Genetic operators

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
        mutate_crem_kwargs_to_work_with = self.mutate_crem_kwargs.copy()
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
            # Bias the searching to similar molecules
            if self.get_similar:
                mol = get_similar_mols(mols=[mol for _, mol in mutants],
                                       ref_mol=self.InitIndividual.mol, pick=1, beta=0.01)[0]
            else:
                _, mol = random.choice(mutants)  # nosec
        except Exception:
            print(f'Note: The mutation on {individual} did not work, it will be returned the same individual')
            mol = individual.mol
        if self.AddHs:
            mol = Chem.AddHs(mol)
        return Individual(mol, randomseed=self.randomseed)

    def pickle(self, title: str, compress: bool = False):
        """Method to pickle the whole GA class

        Parameters
        ----------
        title : str
            Name of the object which will be completed with the corresponding
            extension depending if compress is set to True or False.
        compress : bool, optional
            Use compression, by default False. If True :meth:`moldrug.utils.compressed_pickle` will be used;
            if not :meth:`moldrug.utils.full_pickle` will be used instead.
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


if __name__ == '__main__':
    pass
