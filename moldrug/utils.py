#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs, Lipinski, Descriptors
from meeko import MoleculePreparation, PDBQTMolecule
from crem.crem import mutate_mol, grow_mol

from copy import deepcopy
from inspect import getfullargspec
import multiprocessing as mp
import tempfile, subprocess, random, time, shutil, tqdm, bz2, pickle, _pickle as cPickle, numpy as np, pandas as pd
from typing import List, Dict

from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

######################################################################################################################################
#                                   Here are some important functions to work with                                                   #
######################################################################################################################################


def timeit(method: object):
    """Calculate the time of a process.
    Useful as decorator of functions

    Parameters
    ----------
    method : object
        A python function
    """
    def timed(*args, **kw):
        ts = time.time()
        result = method(*args, **kw)
        te = time.time()
        if 'log_time' in kw:
            name = kw.get('log_name', method.__name__.upper())
            kw['log_time'][name] = int((te - ts) * 1000)
        else:
            print('%r  %2.2f ms' % (method.__name__, (te - ts) * 1000))
        return result
    return timed


def run(command: str, shell: bool = True, executable: str = '/bin/bash', Popen: bool = False):
    """This function is just a useful wrapper around subprocess.Popen, subprocess.run

    Parameters
    ----------
    command : str
        Any command to execute.
    shell : bool, optional
        keyword of ``subprocess.Popen`` and ``subprocess.Popen``, by default True
    executable : str, optional
        keyword of ``subprocess.Popen`` and ``subprocess.Popen``, by default '/bin/bash'
    Popen : bool, optional
        If True it will launch popen if not Run, by default False

    Returns
    -------
    object
        The processes returned by Popen or Run.

    Raises
    ------
    RuntimeError
        In case of non-zero exit status on the provided command.
    """
    if Popen:
        # In this case you could access the pid as: run.pid
        process = subprocess.Popen(command, shell=shell, executable=executable, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    else:
        process = subprocess.run(command, shell=shell, executable=executable, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        returncode = process.returncode
        if returncode != 0:
            print(f'Command {command} returned non-zero exit status {returncode}')
            raise RuntimeError(process.stderr)
    return process


def confgen(smiles: str, return_mol: bool = False):
    """Create a 3D model from a smiles and return a pdbqt string and, a mol if ``return_mol = True``.

    Parameters
    ----------
    smiles : str
        A valid SMILES string.
    return_mol : bool, optional
        If true the function will also return the ``rdkit.Chem.rdchem.Mol``, by default False

    Returns
    -------
    tuple or str
        If ``return_mol = True`` it will return a tuple ``(str[pdbqt], Chem.rdchem.Mol)``, if not only a ``str`` that represents the pdbqt.
    """
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(mol)
    AllChem.MMFFOptimizeMolecule(mol)
    preparator = MoleculePreparation()
    preparator.prepare(mol)
    pdbqt_string = preparator.write_pdbqt_string()
    if return_mol:
        return (pdbqt_string, mol)
    else:
        return pdbqt_string


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
        probs = np.exp(beta*similarities) / np.sum(np.exp(beta*similarities))
        cumsum = np.cumsum(probs)
        indexes = set()
        while len(indexes) != pick:
            r = sum(probs)*np.random.rand()
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
    filter = {
        'NumHAcceptors': {'method': Lipinski.NumHAcceptors, 'cutoff': 10},
        'NumHDonors': {'method': Lipinski.NumHDonors, 'cutoff': 5},
        'wt': {'method': Descriptors.MolWt, 'cutoff': 500},
        'MLogP': {'method': Descriptors.MolLogP, 'cutoff': 5},
        'NumRotatableBonds': {'method': Lipinski.NumRotatableBonds, 'cutoff': 10},
        'TPSA': {'method': Chem.MolSurf.TPSA, 'cutoff': 140},
    }
    cont = 0
    for property in filter:

        if filter[property]['method'](mol) > filter[property]['cutoff']:
            cont += 1
        if cont >= maxviolation:
            return False
    return True

def lipinski_profile(mol:Chem.rdchem.Mol):
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
        'NumHAcceptors': {'method':Lipinski.NumHAcceptors,'cutoff':10},
        'NumHDonors': {'method':Lipinski.NumHDonors,'cutoff':5},
        'wt': {'method':Descriptors.MolWt,'cutoff':500},
        'MLogP': {'method':Descriptors.MolLogP,'cutoff':5},
        'NumRotatableBonds': {'method':Lipinski.NumRotatableBonds,'cutoff':10},
        'TPSA': {'method':Chem.MolSurf.TPSA,'cutoff':range(0,140)},

        'FractionCSP3': {'method':Lipinski.FractionCSP3,'cutoff':None},
        'HeavyAtomCount': {'method':Lipinski.HeavyAtomCount,'cutoff':None},
        'NHOHCount': {'method':Lipinski.NHOHCount,'cutoff':None},
        'NOCount': {'method':Lipinski.NOCount,'cutoff':None},
        'NumAliphaticCarbocycles': {'method':Lipinski.NumAliphaticCarbocycles,'cutoff':None},
        'NumAliphaticHeterocycles': {'method':Lipinski.NumAliphaticHeterocycles,'cutoff':None},
        'NumAliphaticRings': {'method':Lipinski.NumAliphaticRings,'cutoff':None},
        'NumAromaticCarbocycles': {'method':Lipinski.NumAromaticCarbocycles,'cutoff':None},
        'NumAromaticHeterocycles': {'method':Lipinski.NumAromaticHeterocycles,'cutoff':None},
        'NumAromaticRings': {'method':Lipinski.NumAromaticRings,'cutoff':None},
        'NumHeteroatoms': {'method':Lipinski.NumHeteroatoms,'cutoff':None},
        'NumSaturatedCarbocycles': {'method':Lipinski.NumSaturatedCarbocycles,'cutoff':None},
        'NumSaturatedHeterocycles': {'method':Lipinski.NumSaturatedHeterocycles,'cutoff':None},
        'NumSaturatedRings': {'method':Lipinski.NumSaturatedRings,'cutoff':None},
        'RingCount': {'method':Lipinski.RingCount,'cutoff':None},
    }
    profile = {}
    for property in properties:
        profile[property] = properties[property]['method'](mol)
        #print(f"{property}: {properties[property]['method'](mol)}. cutoff: {properties[property]['cutoff']}")
    return profile

#Desirability. doi:10.1016/j.chemolab.2011.04.004; https://www.youtube.com/watch?v=quz4NW0uIYw&list=PL6ebkIZFT4xXiVdpOeKR4o_sKLSY0aQf_&index=3
def LargerTheBest(Value:float, LowerLimit:float, Target:float, r:float = 1)  -> float:
    """Desirability function used when larger values are the targets. If Value is higher or equal than the target it will return 1; if it is lower than LowerLimit it will return 0; else a number between 0 and 1.

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
        return ((Value -LowerLimit)/(Target-LowerLimit))**r
    else:
        return 1.0

def SmallerTheBest(Value:float, Target:float, UpperLimit:float, r:float = 1) -> float:
    """Desirability function used when lower values are the targets. If Value is lower or equal than the target it will return 1; if it is higher than UpperLimit it will return 0; else a number between 0 and 1.

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
        return ((UpperLimit-Value)/(UpperLimit-Target))**r
    else:
        return 0.0

def NominalTheBest(Value:float, LowerLimit:float, Target:float, UpperLimit:float, r1:float = 1, r2:float = 1) -> float:
    """Desirability function used when a target value is desired. If Value is lower or equal than the LowerLimit it will return 0; as well values higher or equal than  UpperLimit; else a number between 0 and 1.

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
        This is the exponent of the interpolation from LowerLimit to Target. Could be used to control the interpolation, by default 1
    r2 : float, optional
        This is the exponent of the interpolation from Target to UpperLimit. Could be used to control the interpolation, by default 1

    Returns
    -------
    float
        A number between 0 and 1. Been 1 the desireable value to get.
    """
    if Value < LowerLimit:
        return 0.0
    elif Value <= Target:
        return ((Value -LowerLimit)/(Target-LowerLimit))**r1
    elif Value <= UpperLimit:
        return ((UpperLimit-Value)/(UpperLimit-Target))**r2
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
def full_pickle(title:str, data:object):
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

def loosen(file:str):
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

def compressed_pickle(title:str, data:object):
    """Compress python object. First cPickle it and then bz2.BZ2File compressed it.

    Parameters
    ----------
    title : str
         Name of the file without extensions, .pbz2 will be added by default
    data : object
        Any serializable python object
    """
    with bz2.BZ2File(f'{title}.pbz2', 'w') as f:
        cPickle.dump(data, f)

def decompress_pickle(file:str):
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
    data = cPickle.load(data)
    return data

def import_sascorer():
    """Function to import sascorer from RDConfig.RDContribDir of RDKit

    Returns
    -------
    module
        The sascorer module ready to use.
    """
    # In order to import sascorer from RDConfig.RDContribDir
    from rdkit.Chem import RDConfig
    import os, importlib.util as importlib_util
    spec=importlib_util.spec_from_file_location('sascorer', os.path.join(RDConfig.RDContribDir, 'SA_Score', 'sascorer.py'))
    sascorer = importlib_util.module_from_spec(spec)
    spec.loader.exec_module(sascorer)
    return sascorer
######################################################################################################################################




######################################################################################################################################
#                                               Classes to work with Vina                                                            #
######################################################################################################################################
class Atom:
    """This is a simple class to wrap a pdbqt Atom. It is based on https://userguide.mdanalysis.org/stable/formats/reference/pdbqt.html#writing-out.
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

    def write(self, name = None):
        if name:
            with open(name,"w") as f:
                f.writelines(self.chunk)
        else:
            with open(f"Run_{self.run}.pdbqt","w") as f:
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

    def BestEnergy(self, write = False):
        min_chunk = min(self.chunks, key= lambda x: x.freeEnergy)
        if write: min_chunk.write("best_energy.pdbqt")
        return min_chunk

######################################################################################################################################




######################################################################################################################################
#                                             Classes to work with moldrug                                                                   #
######################################################################################################################################

class Individual:
    """
    Base class to work with GA, Local and all the fitness functions.
    Individual is a mutable object. It uses the smiles string for '=='
    operator and the cost attribute for arithmetic operations.
    Known issue, in case that we would like to use a numpy array of individuals. It is needed to change the ditype of the generated arrays
    In order to use other operations, or cast to a list
    array = np.array([c1,c2])
    array_2 = (array*2).astype('float64')
    It also admit copy and deepcopy operations
    """
    def __init__(self,smiles:str = None, mol:Chem.rdchem.Mol = None, idx:int = 0, pdbqt:str = None, cost:float = np.inf) -> None:
        """This is the constructor of the class.

        Parameters
        ----------
        smiles : str, optional
            A valid RDKit SMILES, by default None
        mol : Chem.rdchem.Mol, optional
            A valid RDKit molecule. If not provided it will be generated from smiles, by default None
        idx : int, optional
            An identification, by default 0
        pdbqt : str, optional
            A valid pdbqt string. If it is not provided it will be generated from mol through utils.confgen and the mol attribute will be update with the 3D model, by default None
        cost : float, optional
            This attribute is used to perform operations between Individuals and should be used for the cost functions, by default np.inf
        """
        self.smiles = smiles

        if not mol:
            try:
                self.mol = Chem.MolFromSmiles(smiles)
            except Exception:
                self.mol = None
        else:
            self.mol = mol

        self.idx = idx

        if not pdbqt:
            try:
                self.pdbqt, mol3D = confgen(smiles, return_mol=True)
                if not mol:
                    self.mol = Chem.RemoveHs(mol3D)
            except Exception:
                self.pdbqt = None
        else:
            self.pdbqt = pdbqt

        self.cost = cost

    def __repr__(self):
        return f"{self.__class__.__name__}(idx = {self.idx}, smiles = {self.smiles}, cost = {self.cost})"

    def __eq__(self, other: object) -> bool:
        if self.__class__ is other.__class__:
            return self.smiles == other.smiles #self.cost == other.cost and
        else:
            return False # self.smiles == other

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
        return  other  + self.cost

    def __sub__(self, other: object):
        return self.cost - other
    def __rsub__(self, other: object):
        return  other  - self.cost

    def __mul__(self, other: object):
        return self.cost * other
    def __rmul__(self, other: object):
        return  other  * self.cost

    def __truediv__(self, other: object):
        return self.cost / other
    def __rtruediv__(self, other: object):
        return  other  / self.cost

    def __floordiv__(self, other: object):
        return self.cost // other
    def __rfloordiv__(self, other: object):
        return  other  // self.cost

    def __mod__(self, other: object):
        return self.cost % other
    def __rmod__(self, other: object):
        return  other  % self.cost

    def __divmod__(self, other: object):
        return (self.cost // other, self.cost % other)
    def __rdivmod__(self, other: object):
        return  (other  // self.cost, other  % self.cost)

    def __pow__(self, other: object):
        return self.cost ** other
    def __rpow__(self, other: object):
        return  other  ** self.cost

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

def make_sdf(individuals:List[Individual], sdf_name = 'out'):
    """This function create a sdf file from a list of Individuals based on their pdbqt (or pdbqts) attribute
    This assume that the cost function update the pdbqt attribute after the docking with the conformations obtained
    In the case of multiple receptor a new attribute named pdbqts should been added and it is only a list of pdbqt valid string.
    Here will export several sdf depending how many pdbqt string are in the pdbqts attribute.

    Parameters
    ----------
    individuals : list[Individual]
        A list of individuals
    sdf_name : str, optional
        The name for the output file. Could be a ``path + sdf_name``. The sdf extension will be added by the function, by default 'out'

    Example
    -------
    .. ipython:: python

        import tempfile, os
        from moldrug import utils
        # Create some temporal dir
        tmp_path = tempfile.TemporaryDirectory()
        # Creating two individuals
        I1 = utils.Individual('CCCCl')
        I2 = utils.Individual('CCOCCCF')
        # Creating the pdbqts attribute with the pdbqt attribute (this is just a silly example)
        I1.pdbqts = [I1.pdbqt, I1.pdbqt]
        I2.pdbqts = [I2.pdbqt, I2.pdbqt]
        utils.make_sdf([I1, I2], sdf_name = os.path.join(tmp_path.name, 'out'))
        # Two files were created
        # In the other hand, if the attribute pdbqts is not present, only one file is going to be created
        # Delete the pdbqts attribute
        delattr(I1, 'pdbqts')
        delattr(I2, 'pdbqts')
        utils.make_sdf([I1, I2], sdf_name = os.path.join(tmp_path.name, 'out'))
        # Only one file will be created if the pdbqts has not len in some of the individuals or they presents different lens as well. In this case the pdbqts will be completely ignored and pdbqt attribute will be used for the construction of the sdf file
        I1.pdbqts = [I1.pdbqt, I1.pdbqt, I1.pdbqt]
        I2.pdbqts = [I2.pdbqt, I2.pdbqt]
        utils.make_sdf([I1, I2], sdf_name = os.path.join(tmp_path.name, 'out'))
    """
    pdbqt_tmp = tempfile.NamedTemporaryFile(suffix='.pdbqt')

    # Check for the attribute pdbqts in all passed individuals and that all of them have the same number of pdbqt
    check = True
    NumbOfpdbqt = set()
    for individual in individuals:
        if 'pdbqts' in dir(individual):
            NumbOfpdbqt.add(len(individual.pdbqts))
        else:
            check = False
            break
    if len(NumbOfpdbqt) == 0 or len(NumbOfpdbqt) > 1:
        check = False

    if check == True:
        for i in range(list(NumbOfpdbqt)[0]):
            with Chem.SDWriter(f"{sdf_name}_{i+1}.sdf") as w:
                for individual in individuals:
                    with open(pdbqt_tmp.name, 'w') as f:
                        f.write(individual.pdbqts[i])
                    pdbqt_mol = PDBQTMolecule.from_file(pdbqt_tmp.name, skip_typing=True)
                    mol = pdbqt_mol.export_rdkit_mol()
                    mol.SetProp("_Name",f"idx :: {individual.idx}, smiles :: {individual.smiles}, cost :: {individual.cost}")
                    w.write(mol)
            print(f" File {sdf_name}_{i+1}.sdf was createad!")
    else:
        with Chem.SDWriter(f"{sdf_name}.sdf") as w:
            for individual in individuals:
                with open(pdbqt_tmp.name, 'w') as f:
                    f.write(individual.pdbqt)
                pdbqt_mol = PDBQTMolecule.from_file(pdbqt_tmp.name, skip_typing=True)
                mol = pdbqt_mol.export_rdkit_mol()
                mol.SetProp("_Name",f"idx :: {individual.idx}, smiles :: {individual.smiles}, cost :: {individual.cost}")
                w.write(mol)
        print(f"File {sdf_name}.sdf was createad!")

class Local:
    """For local search
    """
    def __init__(self, mol:Chem.rdchem.Mol, crem_db_path:str, costfunc:object, grow_crem_kwargs:Dict = {}, costfunc_kwargs:Dict = {}) -> None:
        self.mol = mol
        self.InitIndividual = Individual(Chem.MolToSmiles(self.mol), self.mol, idx = 0)
        if not self.InitIndividual.pdbqt:
            raise Exception(f"For some reason, it was not possible to create the class Individula was not able to create a pdbqt from the seed_smiles. Consider to check the validity of the SMILES string!")
        self.crem_db_path = crem_db_path
        self.grow_crem_kwargs = grow_crem_kwargs
        self.costfunc = costfunc
        self.costfunc_kwargs = costfunc_kwargs
        self.pop = [self.InitIndividual]



    def __call__(self, njobs:int = 1, pick:int = None):
        self.grow_crem_kwargs.update({'return_mol':True})
        new_mols = list(grow_mol(
            self.mol,
            self.crem_db_path,
            **self.grow_crem_kwargs
            ))
        if pick:
            random.shuffle(new_mols)
            new_mols = new_mols[:pick]
            new_mols = [Chem.RemoveHs(item[1]) for item in new_mols]

        idx0 = len(self.pop)
        for i, mol in enumerate(new_mols):
            individual = Individual(Chem.MolToSmiles(mol), mol, idx = idx0 + i)
            if individual.pdbqt:
                self.pop.append(individual)

        # Calculating cost of each individual
        # Creating the arguments
        args_list = []
        # Make a copy of the self.costfunc_kwargs
        kwargs_copy = self.costfunc_kwargs.copy()
        if 'wd' in getfullargspec(self.costfunc).args:
            costfunc_jobs = tempfile.TemporaryDirectory(prefix='costfunc')
            kwargs_copy['wd'] = costfunc_jobs.name

        for individual in self.pop:
            args_list.append((individual, kwargs_copy))

        print(f'Calculating cost function...')
        pool = mp.Pool(njobs)
        self.pop = [individual for individual in tqdm.tqdm(pool.imap(self.__costfunc__, args_list), total=len(args_list))]
        pool.close()

        if 'wd' in getfullargspec(self.costfunc).args: shutil.rmtree(costfunc_jobs.name)

    def __costfunc__(self, args_list):
        Individual, kwargs = args_list
        #This is just to use the progress bar on pool.imap
        return self.costfunc(Individual, **kwargs)

    def pickle(self,title, compress = False):
        cls = self.__class__
        result = cls.__new__(cls)
        result.__dict__.update(self.__dict__)
        if compress:
            compressed_pickle(title, result)
        else:
            full_pickle(title, result)

    def to_dataframe(self):
        list_of_dictionaries = []
        for individual in self.pop:
            dictionary = individual.__dict__.copy()
            del dictionary['mol']
            list_of_dictionaries.append(dictionary)
        return pd.DataFrame(list_of_dictionaries)

class GA:
    """An implementation of genetic algorithm to search in the chemical space.

    Attributes
    ----------
    -   pop:
    -   SawIndividuals:
    -   NumCalls:
    -   NumGens:

    """
    def __init__(self, seed_smiles:str, costfunc:object, costfunc_kwargs:Dict, crem_db_path:str, maxiter:int, popsize:int, beta:float = 0.001, pc:float =1, get_similar:bool = False, mutate_crem_kwargs:Dict = {}, save_pop_every_gen:int = 0, deffnm:str = 'ga') -> None:
        self.InitIndividual = Individual(seed_smiles, idx=0)
        if not self.InitIndividual.pdbqt:
            raise Exception(f"For some reason, it was not possible to create the class Individula was not able to create a pdbqt from the seed_smiles. Consider to check the validity of the SMILES string!")
        self.costfunc = costfunc
        self.crem_db_path = crem_db_path
        self.pop = [self.InitIndividual]

        self.maxiter = maxiter
        self.popsize = popsize
        self.beta = beta
        self.costfunc_kwargs = costfunc_kwargs
        if 'ncores' in costfunc_kwargs:
            self.costfunc_ncores = self.costfunc_kwargs['ncores']
        else:
            self.costfunc_ncores = 1

        self.nc = round(pc*popsize)
        self.get_similar = get_similar
        self.mutate_crem_kwargs = {
            'radius':3,
            'min_size':1,
            'max_size':8,
            'min_inc':-5,
            'max_inc':3,
            'ncores':1,
        }
        self.mutate_crem_kwargs.update(mutate_crem_kwargs)

        # Saving parameters
        self.save_pop_every_gen = save_pop_every_gen
        self.deffnm = deffnm

        # Tracking parameters
        self.NumCalls = 0
        self.NumGens = 0
        self.SawIndividuals = []

    @timeit
    def __call__(self, njobs:int = 1):
        # Counting the calls
        self.NumCalls += 1

        # Here we will update if needed some parameters for the crem operations that could change between differents calls.
        # We need to return the molecule, so we override the possible user definition respect to this keyword
        self.mutate_crem_kwargs['return_mol'] = True

        # Guess if we need to add explicit hydrogens in the mutate operation
        self.AddExplicitHs = False
        if 'min_size' in self.mutate_crem_kwargs:
            if self.mutate_crem_kwargs['min_size'] == 0:
                self.AddExplicitHs = True
        elif 'max_size' in self.mutate_crem_kwargs:
            if self.mutate_crem_kwargs['max_size'] == 0:
                self.AddExplicitHs = True

        # Initialize Population
        # In case that the populating exist there is not need to initialize.
        if len(self.pop) == 1:
            if self.get_similar:
                # Bias the searching to similar molecules
                GenInitStructs = list(
                    grow_mol(
                        Chem.AddHs(self.InitIndividual.mol),
                        self.crem_db_path,
                        radius=3,
                        min_atoms=1, max_atoms = 4,
                        return_mol= True,
                        ncores = self.mutate_crem_kwargs['ncores']
                        )
                    )
                GenInitStructs = get_similar_mols(mols = [Chem.RemoveHs(item[1]) for item in GenInitStructs], ref_mol=self.InitIndividual.mol, pick=self.popsize, beta=0.01)

            else:
                GenInitStructs = list(
                    mutate_mol(
                        Chem.AddHs(self.InitIndividual.mol),
                        self.crem_db_path,
                        **self.mutate_crem_kwargs,
                        )
                    )
                GenInitStructs = [Chem.RemoveHs(mol) for (_, mol) in GenInitStructs]

            # Checking for possible scenarios
            if len(GenInitStructs) < (self.popsize - 1):
                print('The initial population has repeated elements')
                # temporal solution
                GenInitStructs +=  random.choices(GenInitStructs, k = self.popsize - len(GenInitStructs) -1)
            elif len(GenInitStructs) > (self.popsize - 1):
                #Selected random sample from the generation
                GenInitStructs = random.sample(GenInitStructs, k = self.popsize -1)
            else:
                # Everything is ok!
                pass

            for i, mol in enumerate(GenInitStructs):
                individual = Individual(Chem.MolToSmiles(mol), mol, idx = i + 1) # 0 is the InitIndividual
                if individual.pdbqt:
                    self.pop.append(individual)

            # Calculating cost of each individual
            # Creating the arguments
            args_list = []
            # Make a copy of the self.costfunc_kwargs
            kwargs_copy = self.costfunc_kwargs.copy()
            if 'wd' in getfullargspec(self.costfunc).args:
                costfunc_jobs = tempfile.TemporaryDirectory(prefix='costfunc')
                kwargs_copy['wd'] = costfunc_jobs.name

            for individual in self.pop:
                args_list.append((individual, kwargs_copy))

            print(f'\n\nCreating the first population with {self.popsize} members:')
            pool = mp.Pool(njobs)
            self.pop = [individual for individual in tqdm.tqdm(pool.imap(self.__costfunc__, args_list), total=len(args_list))]
            pool.close()

            if 'wd' in getfullargspec(self.costfunc).args: shutil.rmtree(costfunc_jobs.name)

            # Print some information of the initial population
            print(f"Initial Population: Best individual: {min(self.pop)}")
            # Updating the info of the first individual (parent) to print at the end how well performed the method (cost function)
            # Because How the population was initialized and because we are using pool.imap (ordered). The parent is the first Individual of self.pop.
            # We have to use deepcopy because Individual is a mutable object
            self.InitIndividual = deepcopy(self.pop[0])

            # Best Cost of Iterations
            self.bestcost = []
            self.avg_cost = []

        # Saving tracking variables, the first population, outside the if to take into account second calls with different population provided by the user.
        for individual in self.pop:
            if individual not in self.SawIndividuals:
                self.SawIndividuals.append(individual)

        # Saving population in disk if it was required
        if self.save_pop_every_gen:
            compressed_pickle(f"{self.deffnm}_pop", (self.NumGens,self.pop))
            make_sdf(self.pop, sdf_name=f"{self.deffnm}_pop")

        # Main Loop
        number_of_previous_generations = len(self.bestcost) # Another control variable. In case that the __call__ method is used more than ones.
        for iter in range(self.maxiter):
            # Saving Number of Generations
            self.NumGens += 1

            # Probabilities Selections
            factors = (-self.beta * np.array(self.pop)).astype('float64')
            probs = np.exp(factors) / np.exp(factors).sum()

            popc = []
            for _ in range(self.nc):
                # Perform Roulette Wheel Selection
                parent = self.pop[self.roulette_wheel_selection(probs)]

                # Perform Mutation (this mutation is some kind of crossover but with CReM library)
                children = self.mutate(parent)

                # Save offspring population
                # I will save only those offsprings that were not seen and that have a correct pdbqt file
                if children not in self.SawIndividuals and children.pdbqt: popc.append(children)

            if popc: # Only if there are new members
                # Calculating cost of each offspring individual (Doing Docking)
                pool = mp.Pool(njobs)
                # Creating the arguments
                args_list = []
                # Make a copy of the self.costfunc_kwargs
                kwargs_copy = self.costfunc_kwargs.copy()
                if 'wd' in getfullargspec(self.costfunc).args:
                    costfunc_jobs = tempfile.TemporaryDirectory(prefix='costfunc')
                    kwargs_copy['wd'] = costfunc_jobs.name

                NumbOfSawIndividuals = len(self.SawIndividuals)
                for (i, individual) in enumerate(popc):
                    # Add idx label to each individual
                    individual.idx = i + NumbOfSawIndividuals
                    # The problem here is that we are not being general for other possible Cost functions.
                    args_list.append((individual,kwargs_copy))
                print(f'Evaluating generation {self.NumGens} / {self.maxiter + number_of_previous_generations}:')

                # Calculating cost fucntion in parallel
                popc = [individual for individual in tqdm.tqdm(pool.imap(self.__costfunc__, args_list), total=len(args_list))]
                pool.close()
                if 'wd' in getfullargspec(self.costfunc).args: shutil.rmtree(costfunc_jobs.name)

            # Merge, Sort and Select
            self.pop += popc
            self.pop = sorted(self.pop)
            self.pop = self.pop[:self.popsize]

            # Store Best Cost
            self.bestcost.append(self.pop[0].cost)

            # Store Average cost
            self.avg_cost.append(np.mean(self.pop))

            # Saving tracking variables
            for individual in popc:
                #Tracking variables. There is not need to check if the individual is in self.SawIndividual. It was already checked above.
                self.SawIndividuals.append(individual)

            # Saving population in disk if it was required
            if self.save_pop_every_gen:
                # Save every save_pop_every_gen and always the last population
                if self.NumGens % self.save_pop_every_gen == 0 or iter + 1 == self.maxiter:
                    compressed_pickle(f"{self.deffnm}_pop", (self.NumGens, self.pop))
                    make_sdf(self.pop, sdf_name=f"{self.deffnm}_pop")

            # Show Iteration Information
            print(f"Generation {self.NumGens}: Best Individual: {self.pop[0]}.\n")

        # Printing summary information
        print(f"\n{50*'=+'}\n")
        print(f'The simulation finished successfully after {self.maxiter} generations with a population of {self.popsize} individuals. A total number of {len(self.SawIndividuals)} Individuals were seen during the simulation.')
        print(f"Initial Individual: {self.InitIndividual}")
        print(f"Final Individual: {self.pop[0]}")
        print(f"The cost fucntion droped in {self.InitIndividual - self.pop[0]} units.")
        print(f"\n{50*'=+'}\n")

    def __costfunc__(self, args_list):
        individual, kwargs = args_list
        #This is just to use the progress bar on pool.imap
        return self.costfunc(individual, **kwargs)

    def mutate(self, individual):
        try:
            if self.AddExplicitHs:
                mutants = list(mutate_mol(Chem.AddHs(individual.mol), self.crem_db_path, **self.mutate_crem_kwargs))
                # Remove the Hs and place a dummy 0 as first element in the tuple.
                mutants = [(0, Chem.RemoveHs(mutant)) for _, mutant in mutants]
            else:
                mutants = list(mutate_mol(individual.mol, self.crem_db_path, **self.mutate_crem_kwargs))

            # Bias the searching to similar molecules
            if self.get_similar:
                mol = get_similar_mols(mols = [mol for _, mol in mutants], ref_mol=self.InitIndividual.mol, pick=1, beta=0.01)[0]
                smiles = Chem.MolToSmiles(mol)
            else:
                _, mol = random.choice(mutants)
                smiles = Chem.MolToSmiles(mol)
        except Exception:
            print('The mutation did not work, we returned the same individual')
            smiles, mol = individual.smiles, individual.mol
        return Individual(smiles,mol)


    def roulette_wheel_selection(self, p):
        c = np.cumsum(p)
        r = sum(p)*np.random.rand()
        ind = np.argwhere(r <= c)
        return ind[0][0]

    def pickle(self, title, compress = False):
        cls = self.__class__
        result = cls.__new__(cls)
        result.__dict__.update(self.__dict__)
        if compress:
            compressed_pickle(title, result)
        else:
            full_pickle(title, result)

    def to_dataframe(self):
        list_of_dictionaries = []
        for individual in self.SawIndividuals:
            dictionary = individual.__dict__.copy()
            del dictionary['mol']
            list_of_dictionaries.append(dictionary)
        return pd.DataFrame(list_of_dictionaries)
######################################################################################################################################

if __name__ == '__main__':
    pass