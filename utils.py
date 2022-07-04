#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from rdkit import Chem
from rdkit.Chem import AllChem, rdmolops, DataStructs, Lipinski, Descriptors
from openbabel import openbabel as ob
from crem.crem import mutate_mol, grow_mol

from copy import deepcopy
from inspect import getfullargspec
import multiprocessing as mp
import tempfile, subprocess, os, random, time, shutil, tqdm, warnings, bz2, pickle, _pickle as cPickle, numpy as np, pandas as pd
from datetime import datetime
from sklearn.ensemble import RandomForestRegressor

from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*') 

# Alias for vina
vina_executable = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'vina')

######################################################################################################################################
#                                   Here are some important functions to work with                                                   #
######################################################################################################################################
def timeit(method):
    """Calculate the time of a process. Useful as decorator of functions

    Args:
        method (function): Any function.
    """
    def timed(*args, **kw):
        ts = time.time()
        result = method(*args, **kw)
        te = time.time()
        if 'log_time' in kw:
            name = kw.get('log_name', method.__name__.upper())
            kw['log_time'][name] = int((te - ts) * 1000)
        else:
            print('%r  %2.2f ms' % \
                  (method.__name__, (te - ts) * 1000))
        return result
    return timed

def run(command:str, shell:bool = True, executable:str = '/bin/bash', Popen:bool = False):
    """This function is just a useful wrapper around subprocess.Popen, subprocess.run

    Args:
        command (str): Any command to execute on Linux
        shell (bool, optional): keyword of Popen and Run. Defaults to True.
        executable (str, optional): keyword of Popen and Run. Defaults to '/bin/bash'.
        Popen (bool, optional): If True it will launch popen if not Run. Defaults to False.

    Raises:
        RuntimeError: In case of non-zero exit status.

    Returns:
        object: the processes returned by Popen or Run.
    """
    if Popen:
        #In this case you could acces the pid as: run.pid
        process = subprocess.Popen(command, shell = shell, executable = executable, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text = True)
    else:
        process = subprocess.run(command, shell = shell, executable = executable, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text = True)
        returncode = process.returncode
        if returncode != 0:
            print(f'Command {command} returned non-zero exit status {returncode}')
            raise RuntimeError(process.stderr)
    return process

def obconvert(inpath:str, outpath:str):
    """Convert a molecule using openbabel

    Args:
        input (str, path): input molecule.
        output (str, path): must have the extension of the molecule.
    """
    in_ext = os.path.basename(inpath).split('.')[-1]
    out_ext = os.path.basename(outpath).split('.')[-1]

    obConversion = ob.OBConversion()
    obConversion.SetInAndOutFormats(in_ext, out_ext)
    mol = ob.OBMol()
    obConversion.ReadFile(mol, inpath)   # Open Babel will uncompressed automatically
    #mol.AddHydrogens()
    obConversion.WriteFile(mol, outpath)

def confgen(smiles:str, outformat:str = "pdbqt"):
    """Create a 3D model from a smiles and return in the specified format.

    Args:
        smiles (str): a valid smiles code.
        outformat (str, optional): The output molecule extension. Defaults to "pdbqt".
    """
    molin = tempfile.NamedTemporaryFile(suffix='.mol')
    molout = tempfile.NamedTemporaryFile(suffix=f'.{outformat}')
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(mol)
    AllChem.MMFFOptimizeMolecule(mol)
    Chem.MolToMolFile(mol, molin.name)
    obconvert(molin.name, molout.name)
    with open(molout.name, 'r') as o:
        string = o.read()
    return string

def fragments(mol:Chem.rdchem.Mol):
    """Create the fragments of the molecule based on the single bonds

    Args:
        mol (Chem.rdchem.Mol): An RDKit molecule.

    Returns:
        Chem.rdmolops.GetMolFrags: The fragments. you could cast the result using list and get all the fragments.
    """
    break_point = int(random.choice(np.where(np.array([b.GetBondType() for b in mol.GetBonds()]) == Chem.rdchem.BondType.SINGLE)[0]))
    # Chem.FragmentOnBonds(mol, break_point) # Could be used to increase randomness give more possible fragments and select two of them
    with Chem.RWMol(mol) as rwmol:
        b = rwmol.GetBondWithIdx(break_point)
        rwmol.RemoveBond(b.GetBeginAtomIdx(), b.GetEndAtomIdx())
    return rdmolops.GetMolFrags(rwmol, asMols = True)

def rdkit_numpy_convert(fp):
    """Convert a list of binary fingerprint to numpy arrays.

    Args:
        fp (_type_): list of binary fingerprints

    Returns:
        numpy.array: An array with dimension (number of elements in the list, number of bits on the fingerprint).
    """
    # fp - list of binary fingerprints
    output = []
    for f in fp:
        arr = np.zeros((1,))
        DataStructs.ConvertToNumpyArray(f, arr)
        output.append(arr)
    return np.asarray(output)

def get_top(ms:list, model):
    """Get the molecule with higher value predicted by the model.

    Args:
        ms (list): list of molecules
        model (sklearn model): A model for with the training set was a set of numpy array based on the AllChem.GetMorganFingerprintAsBitVect 

    Returns:
        _type_: The molecule with the higher predicted value.
    """
    fps1 = [AllChem.GetMorganFingerprintAsBitVect(m, 2) for m in ms]
    x1 = rdkit_numpy_convert(fps1)
    pred = model.predict(x1)
    i = np.argmax(pred)
    return ms[i], pred[i]

def get_sim(ms:list, ref_fps):
    """Get the molecules with higher similarity to each member of ref_fps.

    Args:
        ms (list): list of molecules
        ref_fps (AllChem.GetMorganFingerprintAsBitVect): A list of reference fingerprints

    Returns:
        list: A list of molecules with the higher similarity with their corresponded ref_fps value.
    """
    # ms - list of molecules
    # ref_fps - list of fingerprints of reference molecules
    output = []
    fps1 = [AllChem.GetMorganFingerprintAsBitVect(m, 2) for m in ms]
    for fp in fps1:
        v = DataStructs.BulkTanimotoSimilarity(fp, ref_fps)
        i = np.argmax(v)
        output.append([v[i], i])
    return output

def get_similar_mols(mols:list, ref_mol:Chem.rdchem.Mol, pick:int, beta:float = 0.01):
    """Pick the similar molecules from mols respect to ref_mol using a roulette wheel selection strategy.

    Args:
        mols (list): the list of molecules from where pick molecules will be chosen
        ref_mol (Chem.rdchem.Mol): The reference molecule
        pick (int): Number of molecules to pick from mols
        beta (float, optional): Selection threshold. Defaults to 0.01.

    Returns:
        list: A list of picked molecules.
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

def lipinski_filter(mol:Chem.rdchem.Mol, maxviolation:int = 2):
    """Implementation of Lipinski filter.

    Args:
        mol (Chem.rdchem.Mol): An RDKit molecule.
        maxviolation (int, optional): Maximum number of violations to accept the molecule. Defaults to 2.

    Returns:
        bool: True if the molecule present less than maxviolation of Number of violation; otherwise False.
    """
    filter = {
        'NumHAcceptors': {'method':Lipinski.NumHAcceptors,'cutoff':10},
        'NumHDonors': {'method':Lipinski.NumHDonors,'cutoff':5},
        'wt': {'method':Descriptors.MolWt,'cutoff':500},
        'MLogP': {'method':Descriptors.MolLogP,'cutoff':5},
        'NumRotatableBonds': {'method':Lipinski.NumRotatableBonds,'cutoff':10},
        'TPSA': {'method':Chem.MolSurf.TPSA,'cutoff':140},
    }
    cont = 0
    for property in filter:
        
        if filter[property]['method'](mol) > filter[property]['cutoff']:
            cont += 1
        if cont >= maxviolation:
            return False
    return True

def lipinski_profile(mol:Chem.rdchem.Mol):
    """Get several drug-like properties.
    See: https://www.rdkit.org/docs/source/rdkit.Chem.Lipinski.html?highlight=lipinski#module-rdkit.Chem.Lipinski

    Args:
        mol (Chem.rdchem.Mol): An RDKit molecule.

    Returns:
        dict: A dictionary with the properties.
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

    Args:
        Value (float): Value to test.
        LowerLimit (float): Lower value accepted. Lower than this one will return 0.
        Target (float): The target value. On this value (or higher) the function takes 1 as value.
        r (float, optional): This is the exponent of the interpolation. Could be used to control the interpolation. Defaults to 1.

    Returns:
        float: A number between 0 and 1. Been 1 the desireable value to get.
    """
    if Value < LowerLimit:
        return 0.0
    elif LowerLimit <= Value <= Target:
        return ((Value -LowerLimit)/(Target-LowerLimit))**r
    else:
        return 1.0

def SmallerTheBest(Value:float, Target:float, UpperLimit:float, r:float = 1) -> float:
    """Desirability function used when lower values are the targets. If Value is lower or equal than the target it will return 1; if it is higher than UpperLimit it will return 0; else a number between 0 and 1.

    Args:
        Value (float): Value to test.
        Target (float): The target value. On this value (or lower) the function takes 1 as value.
        UpperLimit (float): Upper value accepted. Higher than this one will return 0.
        r (float, optional): This is the exponent of the interpolation. Could be used to control the interpolation. Defaults to 1.

    Returns:
        float: A number between 0 and 1. Been 1 the desireable value to get.
    """
    if Value < Target:
        return 1.0
    elif Target <= Value <= UpperLimit:
        return ((UpperLimit-Value)/(UpperLimit-Target))**r
    else:
        return 0.0

def NominalTheBest(Value:float, LowerLimit:float, Target:float, UpperLimit:float, r1:float = 1, r2:float = 1) -> float:
    """Desirability function used when a target value is desired. If Value is lower or equal than the LowerLimit it will return 0; as well values higher or equal than  UpperLimit; else a number between 0 and 1.

    Args:
        Value (float):  Value to test.
        LowerLimit (float): Lower value accepted. Lower than this one will return 0.
        Target (float): The target value. On this value the function takes 1 as value.
        UpperLimit (float): Upper value accepted. Higher than this one will return 0.
        r1 (float, optional): This is the exponent of the interpolation from LowerLimit to Target. Could be used to control the interpolation. Defaults to 1.
        r2 (float, optional): This is the exponent of the interpolation from Target to UpperLimit. Could be used to control the interpolation. Defaults to 1.

    Returns:
        float: A number between 0 and 1. Been 1 the desireable value to get.
    """
    if Value < LowerLimit:
        return 0.0
    elif LowerLimit <= Value <= Target:
        return ((Value -LowerLimit)/(Target-LowerLimit))**r1
    elif Target <= Value <= UpperLimit:
        return ((UpperLimit-Value)/(UpperLimit-Target))**r2
    else:
        return 0.0

# Saving data
def full_pickle(title:str, data:object):
    """Normal pickle.

    Args:
        title (str): name of the file without extension, .pkl will be added by default.
        data (object): Any serializable python object
    """
    with open(f'{title}.pkl', 'wb') as pkl:
        pickle.dump(data, pkl)

def loosen(file:str):
    """Unpickle a pickled object.

    Args:
        file (str): The path to the file who store the pickle object.

    Returns:
        object: The python object.
    """
    with open(file, 'rb') as pkl:
        data = pickle.load(pkl)
    return data

def compressed_pickle(title:str, data:object):
    """Compress python object. First cPickle it and then bz2.BZ2File compressed it.

    Args:
        title (str): Name of the file without extensions, .pbz2 will be added by default
        data (object): Any serializable python object
    """
    with bz2.BZ2File(f'{title}.pbz2', 'w') as f: 
        cPickle.dump(data, f)   

def decompress_pickle(file:str):
    """Decompress CPickle objects compressed first with bz2 formats

    Args:
        file (str): This is the cPickle files compressed with bz2.BZ2File. (as a convention with extension .pbz2, but not needed)

    Returns:
        object: The python object.
    """
    data = bz2.BZ2File(file, 'rb')
    data = cPickle.load(data)
    return data
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

class Hetatm:
    """This is a simple class to wrap a pdbqt Hetatm. It is based on https://userguide.mdanalysis.org/stable/formats/reference/pdbqt.html#writing-out.
    """
    def __init__(self, line):
        self.lineType = "HETATM"
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

class Remark:
    """For now usesless
    """
    def __init__(self, line):
        pass

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
            else:
                pass

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
    To acces the chunks you need to take into account that 
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


class ScoringPredictor:
    """This class is used to predict a scoring based on the information of the SMILES.
    """
    def __init__(self, smiles_scoring:dict, receptor:str, boxcenter:list, boxsize:list, exhaustiveness:int) -> None:
        # Esto puede servir para simplemente buscar en la base de datos por la moelcula, y si esta dar el valor exacto de Vina sin tenr que hacer el calculo
        # Y por supuesto para predecir
        # La prediccion en la parte del "crossover" y el si esta la moelcuela para el caclulo real
        self.lastupdate = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        self.smiles_scoring = smiles_scoring
        # This are control variables to save for future use
        self.receptor = receptor
        self.boxcenter = boxcenter
        self.boxsize = boxsize
        self.exhaustiveness = exhaustiveness
        self.bfp = rdkit_numpy_convert([AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(s), 2) for s in self.smiles_scoring])
        self.model = None

    
    def __call__(self, **keywords):
        # n_estimators=100, n_jobs=njobs*self.costfunc_keywords['vina_cpus'], random_state=42, oob_score=True
        self.model = RandomForestRegressor(**keywords)
        self.model.fit(self.bfp, list(self.smiles_scoring.values()))
    
    def update(self, new_smiles_scoring:dict, receptor:str, boxcenter:list, boxsize:list, exhaustiveness:int) -> None:
        # Control that the new introduced data belongs to the same model
        assert self.receptor == receptor, f"The original object was constructed with the receptor {self.receptor} and you introduced {receptor}. Consider to generate a new object."
        assert self.boxcenter == boxcenter, f"The original object was constructed with the boxcenter {self.boxcenter} and you introduced {boxcenter}. Consider to generate a new object."
        assert self.boxsize == boxsize, f"The original object was constructed with the boxsize {self.boxsize} and you introduced {boxsize}. Consider to generate a new object."
        assert self.exhaustiveness == exhaustiveness, f"The original object was constructed with the exhaustiveness {self.exhaustiveness} and you introduced {exhaustiveness}. Consider to generate a new object."
        
        # Update the date 
        self.lastupdate = datetime.now().strftime
        
        # Look for new structures
        smiles_scoring2use = dict()
        for sm_sc in new_smiles_scoring:
            if sm_sc not in self.smiles_scoring:
                smiles_scoring2use[sm_sc] = new_smiles_scoring[sm_sc]
        if smiles_scoring2use:
            self.smiles_scoring.update(smiles_scoring2use)

            self.bfp = np.vstack((self.bfp, rdkit_numpy_convert([AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(s), 2) for s in smiles_scoring2use])))
            self.model = None
            print(f"{len(smiles_scoring2use)} new structures will be incorporate to the model.")
        else:
            print("The introduced smiles are already in the data base. No need to update.")

    def predict(self, smiles:list):
        if self.model and smiles:
            if type(smiles) == str: smiles = [smiles]
            fps = [AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(s), 2) for s in smiles]
            return self.model.predict(rdkit_numpy_convert(fps))
        else:
            if not self.model:
                warnings.warn('There are not model on this object, please call it (__call__) to create it in order to get a prediction.  If you update the class, you must call the class again in order to create a new model based on the updated information, right now was set to None" Right now you just got None!')
                return None
            if not smiles:
                warnings.warn('You did not provide a smiles. You will just get None as prediction')
                return None
    
    def pickle(self,title, compress = False):
        cls = self.__class__
        result = cls.__new__(cls)
        result.__dict__.update(self.__dict__)
        if compress:
            compressed_pickle(title, result)
        else:
            full_pickle(title, result) 
######################################################################################################################################




######################################################################################################################################
#                                             Classes to work lead                                                                   #
######################################################################################################################################
class Individual:
    """
    Base class to work with GA, Local and all the fitness fucntions.
    Individual is a mutable object. It uses the smiles string for '=='
    operator and the cost attribute for arithmetic operations.
    Known issue, in case that we would like to use a numpy array of individuals. It is needed to change the ditype of the generatead arrays
    In order to use other operations, or cast to a list
    array = np.array([c1,c2])
    array_2 = array*2
    array_2 = array_2.astype('float64')
    """
    def __init__(self,smiles:str = None, mol:Chem.rdchem.Mol = None, idx:int = 0, pdbqt = None, cost:float = np.inf) -> None:
        self.smiles = smiles
        
        if not mol:
            try:
                self.mol = Chem.MolFromSmiles(smiles)
            except:
                self.mol = None
        else:
            self.mol = mol
        
        self.idx = idx
        
        if not pdbqt:
            try:
                self.pdbqt = confgen(smiles, outformat = 'pdbqt')
            except:
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

class Local:
    def __init__(self, mol:Chem.rdchem.Mol, crem_db_path:str, costfunc:object, grow_crem_kwargs:dict = {}, costfunc_kwargs:dict = {}) -> None:
        # Add check to get if the molecules is with Hs in case that some specification of the where to grow is given
        self.mol = mol
        self.crem_db_path = crem_db_path
        self.grow_crem_kwargs = grow_crem_kwargs
        self.grow_crem_kwargs.update({'return_mol':True})
        self.costfunc = costfunc
        self.costfunc_kwargs = costfunc_kwargs
        
        MolNonHs = Chem.RemoveHs(self.mol)
        self.pop = [Individual(Chem.MolToSmiles(MolNonHs), MolNonHs, idx = 0)]
    def __call__(self, njobs:int = 1, pick:int = None):
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
            self.pop.append(Individual(Chem.MolToSmiles(mol), mol, idx = idx0 + i))
        
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
        for Individual in self.pop:
            dictionary = Individual.__dict__.copy()
            del dictionary['mol']
            list_of_dictionaries.append(dictionary)
        return pd.DataFrame(list_of_dictionaries)       

## Problem!!
# Another way to overcome the problem on generating similar molecules is to instead of calculate the similarity respect to the whole molecule give some part of the reference structure, give the pharmacophore in the initialization of GA
# and in the cost function, first find the closer fragments and from there calculate the similarity. because always the generated molecules will be different
# The mutation that I am doing right now is more a crossover but with the CReM library instead with the population itself.
# Code a more conservative mutation (close to the one in the paper of the SMILES) to perform small mutations. (now is user controlled)
# We could also, instead a crossover two mutations with different levels. Because our 'mutate' is indeed some kind of crossover but with the CReM data base. So, what we could do is implement the mutation sof the paper that is at small scale (local optimization): the possibilities are point mutations (atom-by-atom), deletion, add new atoms
# hard_mutate and soft_mutate, our genetic operations
# For the best ligand we could predict the metabolites with sygma (apart of the module)
# Think about how to handle possibles Vina Crash (problems with the type of atoms). i think that is some problems happens with vina the cost function will be just np.inf
# Sometimes the pdbqt structure is not generated (probably problems in the conversion. In this cases the whole simulation crash ans should not be like this, this should rise some warning and continue discarding this structure)
# Apply filter for chemical elements to avoid crash on the vina function, in case that Vina crash we could use the predicted model
# Till now Vina never fails but could happen.
# Add more cpus for the generation of the conformers with RDKit
# The size of the ligands increase with the number of generations (if crossover is used even more)
# How to implement the rationality of where to grow, not just randomness. That could be the "crossover" operation, in fact the grow in a specific direction, based in (for example) the interaction network or the clashing avoid.
# I have to create a filter of atoms in order that vina doesn't fail because B and atoms like that Vina is not able to handle.
class GA:
    
    def __init__(self, seed_smiles:str, costfunc:object, crem_db_path:str, maxiter:int, popsize:int, beta:float = 0.001, pc:float =1, get_similar:bool = False, mutate_crem_kwargs:dict = {}, costfunc_kwargs:dict = {}, save_pop_every_gen:int = 0, pop_file_name:int = 'pop') -> None:
        self.InitIndividual = Individual(seed_smiles, idx=0)
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
        # We need to return the molecule, so we override the possible user definition respect to this keyword
        self.mutate_crem_kwargs['return_mol'] = True
        
        # Saving parameters
        self.save_pop_every_gen = save_pop_every_gen
        self.pop_file_name = pop_file_name
        
        # Tracking parameters
        self.NumCalls = 0
        self.NumGen = 0
        self.SawIndividuals = []
    
    @timeit
    def __call__(self, njobs:int = 1, predictor_model:ScoringPredictor = None):
        # Counting the calls 
        self.NumCalls += 1
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
                pass# I am not sure how to deal with this
            elif len(GenInitStructs) > (self.popsize - 1):
                #Selected random sample from the generation 
                GenInitStructs = random.sample(GenInitStructs, k = self.popsize -1)
            else:
                # Everything is ok!
                pass 

            for i, mol in enumerate(GenInitStructs):
                self.pop.append(Individual(Chem.MolToSmiles(mol), mol, idx = i + 1))# 0 is the InitIndividual
            
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


            # Creating the first model
            # Could be more pretty like creating a method that is make the model, y lo que se hace es que se gurda el modelo
            # Aqui se hace por primera vez pero para no repetir tanto codigo solo se llama a update model
            
            # if predictor_model:
            #     self.Predictor = predictor_model
            #     print('\nUpdating the provided model:\n')
            #     self.Predictor.update(
            #         new_smiles_scoring = dict(((individual.smiles,individual.vina_score) if individual.vina_score != np.inf else (individual.smiles,9999) for individual in self.SawIndividuals)),
            #         receptor = os.path.basename(self.costfunc_kwargs['receptor_path']).split('.')[0],
            #         boxcenter = self.costfunc_kwargs['boxcenter'],
            #         boxsize = self.costfunc_kwargs['boxsize'],
            #         exhaustiveness = self.costfunc_kwargs['exhaustiveness'],
            #     )
            #     self.Predictor(n_estimators=100, n_jobs=njobs*self.costfunc_ncores, random_state=42, oob_score=True)
            #     print('Done!')
            # else:
            #     print('\nCreating the first predicted model...')
            #     self.Predictor = vina.ScoringPredictor(
            #         smiles_scoring = dict(((individual.smiles,individual.vina_score) if individual.vina_score != np.inf else (individual.smiles,9999) for individual in self.SawIndividuals)),
            #         receptor = os.path.basename(self.costfunc_kwargs['receptor_path']).split('.')[0],
            #         boxcenter = self.costfunc_kwargs['boxcenter'],
            #         boxsize = self.costfunc_kwargs['boxsize'],
            #         exhaustiveness = self.costfunc_kwargs['exhaustiveness'],
            #     )
            #     self.Predictor(n_estimators=100, n_jobs=njobs*self.costfunc_ncores, random_state=42, oob_score=True)
            #     print('Done!')

            # print(f'The model presents a oob_score = {self.Predictor.model.oob_score_}\n')
            
            # Best Cost of Iterations
            self.bestcost = []
            self.avg_cost = []
        
        # Saving tracking variables, the first population, outside the if to take into account second calls with different population provided by the user.
        for individual in self.pop:
            if individual not in self.SawIndividuals:
                self.SawIndividuals.append(individual)
        
        # Saving population in disk if it was required
        if self.save_pop_every_gen:
            full_pickle(self.pop_file_name, (self.NumGen,self.pop))
        
        # Main Loop
        number_of_previous_generations = len(self.bestcost) # Another control variable. In case that the __call__ method is used more than ones.
        for iter in range(self.maxiter):
            # Saving Number of Generations
            self.NumGen += 1

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
                # I will save only those offsprings that were not seen 
                if children not in self.SawIndividuals: popc.append(children)

            if popc: # Only if there are new members
                # Calculating cost of each offspring individual (Doing Docking)

                #os.makedirs('.vina_jobs', exist_ok=True)
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
                print(f'\nEvaluating generation {self.NumGen} / {self.maxiter + number_of_previous_generations}:')

                #!!!! Here I have to see if the smiles are in the self.saw_smiles in order to do not perform the docking and just assign the scoring function

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

            # Saving tracking variables and getting new ones for the model update
            # new_smiles_cost = dict()
            for individual in popc:
                if individual not in self.SawIndividuals:
                    # # New variables
                    # if individual.cost == np.inf:
                    #     new_smiles_cost[individual.smiles] = 9999
                    # else:
                    #     new_smiles_cost[individual.smiles] = individual.cost
                    #Tracking variables
                    self.SawIndividuals.append(individual)
            
            # Saving population in disk if it was required
            if self.save_pop_every_gen:
                # Save every save_pop_every_gen and always the last population
                if self.NumGen % self.save_pop_every_gen == 0 or iter + 1 == self.maxiter:
                    full_pickle(self.pop_file_name, (self.NumGen, self.pop))

            # # Update the model
            # print(f'Updating the current model with the information of generation {self.NumGen}...')

            # self.Predictor.update(
            #     new_smiles_scoring = new_smiles_cost.copy(),
            #     receptor = os.path.basename(self.costfunc_kwargs['receptor_path']).split('.')[0],
            #     boxcenter = self.costfunc_kwargs['boxcenter'],
            #     boxsize = self.costfunc_kwargs['boxsize'],
            #     exhaustiveness = self.costfunc_kwargs['exhaustiveness'],
            # )
            # self.Predictor(n_estimators=100, n_jobs=njobs*self.costfunc_ncores, random_state=42, oob_score=True)          
            # print('Done!')
            # print(f'The updated model presents a oob_score = {self.Predictor.model.oob_score_}')
            

            # Show Iteration Information
            print(f"Generation {self.NumGen}: Best Individual: {self.pop[0]}.\n")
            # plt.scatter(self.NumGen, self.pop[0].cost)
        
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

    # def crossover(self, individual1, individual2, ncores = 1, probability = 0, MaxRatioOfIncreaseInWt = 0.25):
    #     # here I have to select some randomness to perform or not the real crossover because I think that we could get far from the solution. It is just a guess.
    #     # How do I control the size of the new offspring? 
    #     # Performing a fragmentation in such a way that the offspring is the same in size
    #     # Here is where more additional information could be used. In order to orient the design of the new offspring. 
    #     # Then, I should control how perform the mutation  in such a way that we could keep or at least evaluate the offspring generated for crossover
    #     if random.random() < probability: # 50% of return the same individuals
    #         fragments1 = fragments(individual1.mol)
    #         fragments2 = fragments(individual2.mol)
    #         all_fragments = list(fragments1) + list(fragments2)
            
    #         # Initialize offspring smiles; cost
    #         offsprings = [
    #                 [None, np.inf],
    #                 [None, np.inf],
    #         ]
    #         for combination in itertools.combinations(all_fragments, 2):
                
    #             # Combine the molecules
    #             try:
    #                 possible_fragments_smiles = list(link_mols(*combination, db_name=self.crem_db_path, radius = 3, min_atoms=1, max_atoms=6, dist = 2, return_mol=False, ncores=ncores))                
    #             except:
    #                 # This is for debugging
    #                 sm1, sm2 = [Chem.MolToSmiles(c) for c in combination]
    #                 raise RuntimeError(f'These are the problematic SMILES: {sm1}, {sm2}')
                
    #             # Perform a filter based on weight. This control the size of the fragments. For now I will test 25 %. Think in the future work with the mols instead of smiles, I have to convert to mols too many times in this section of the code
    #             avg_wt  = 0.5*(Descriptors.ExactMolWt(individual1.mol) + Descriptors.ExactMolWt(individual1.mol))
    #             threshold_wt = (MaxRatioOfIncreaseInWt + 1) * avg_wt
    #             print(f'We had {len(possible_fragments_smiles)} possible fragments')
    #             possible_fragments_smiles = list(filter(lambda x: Descriptors.ExactMolWt(Chem.MolFromSmiles(x)) < threshold_wt, possible_fragments_smiles))
    #             print(f'After the weight filter we have {len(possible_fragments_smiles)} possible fragments')

    #             # In case that it was not possible to link the fragments
    #             if not possible_fragments_smiles:continue

    #             # Bias the searching to similar molecules
    #             if self.get_similar:
    #                 possible_fragments_mols = get_similar_mols(mols = [Chem.MolFromSmiles(smiles) for smiles in possible_fragments_smiles], ref_mol=self.InitIndividual.mol, pick=self.popsize, beta=0.01)
    #                 possible_fragments_smiles = [Chem.MolToSmiles(mol) for mol in possible_fragments_mols]
                
    #             # Here comes the prediction with the model, and get the top two
    #             temp_offsprings = list(zip(possible_fragments_smiles, self.Predictor.predict(possible_fragments_smiles).tolist()))
                
    #             # Merge, Sort and Select
    #             offsprings = sorted(offsprings + temp_offsprings, key = lambda x:x[1])[:2]
    #         # Here I should check that exist offsprings (there not None values as smiles). For now I will assume that we always get at least two. See on the future
    #         return Individual(smiles = offsprings[0][0]), Individual(smiles = offsprings[1][0])
    #     else:
    #         return individual1, individual2    
    
    # Improve

    def mutate(self, individual):
        # See the option max_replacment
        # Or select the mutant based on some criterion
        # try:
            # Here i will pick the molecules based on the model.
        # El problema de seleccionar asi los compuestos es que siempre seleccionamos los mismos. Siempre se esta entrando la misma estructura y terminamos con una pobalcion redundante
        # Esto tengo que pensarlo mejor
        # new_mols = list(mutate_mol(Chem.AddHs(individual.mol), self.crem_db_path, radius=3, min_size=1, max_size=8,min_inc=-3, max_inc=3, return_mol=True, ncores = ncores))
        # new_mols = [Chem.RemoveHs(i[1]) for i in new_mols]
        # best_mol, score = get_top(new_mols + [individual.mol], self.model)
        # smiles = Chem.MolToSmiles(best_mol)
        # mol = best_mol
        # print(score)
        # For now I am generating all the mutants and picking only one at random, this is very inefficient, should be better only generate one, but I am afraid that crem generate always the same or not generate any at all.
        # I think that what first crem does is randomly select on spot and find there all possible mutants. If this spot doesn't generate mutants, then you don't get nothing. But this is a supposition. 
        try:
            mutants = list(mutate_mol(individual.mol, self.crem_db_path, **self.mutate_crem_kwargs))
            # Bias the searching to similar molecules
            if self.get_similar:
                mol = get_similar_mols(mols = [mol for _, mol in mutants], ref_mol=self.InitIndividual.mol, pick=1, beta=0.01)[0]
                smiles = Chem.MolToSmiles(mol)
            else:
                smiles, mol = random.choice(mutants)
        except:
            print('The hard mutation did not work, we returned the same individual')
            smiles, mol = individual.smiles, individual.mol
        return Individual(smiles,mol)
    
    
    def roulette_wheel_selection(self, p):
        c = np.cumsum(p)
        r = sum(p)*np.random.rand()
        ind = np.argwhere(r <= c)
        return ind[0][0]
    
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
        for individual in self.SawIndividuals:
            dictionary = individual.__dict__.copy()
            del dictionary['mol']
            list_of_dictionaries.append(dictionary)
        return pd.DataFrame(list_of_dictionaries)
######################################################################################################################################


if __name__ == '__main__':
    i1 = Individual(smiles = 'CCCO', cost = 10)
    i2 = Individual(smiles = 'CCC', cost = 2)
    i3 = Individual(smiles = 'CCCF', cost = 3)
    i4 = Individual(smiles = 'CCCF', cost = 69)
    i5 = Individual(smiles = 'CCCCCCF', cost = 69)
    # print(i1)
    # print(f"i1 == i2 :{i1 == i2}")
    # print(f"i1 != i2: {i1 != i2}")
    # print(f"i1 <= i2: {i1 <= i2}")
    # print(f"i1 >= i2: {i1 >= i2}")
    # print(f"i1 < i2: {i1 < i2}")
    # print(f"i1 > i2: {i1 > i2}")
    # print(f"i1 is i2: {i1 is i2}")
    # print(f"i1 + i2: {i2 + i1}")
    # print(f"i1 * i2: {i1 * i1}")
    print(f"i1 / i2: { i5//6}")
    print(-i1)
    # print(f"np.mean: {np.mean([i1, i2])}")
    # print(f"bool: {bool(i2)}")
    # print(f"min {min([i1,i2,i3])}")
    # print(f"i4 in [] { 'CCCF' in  [i1,i2,i3]}")
    
    # for i in [i4, i5, Individual('CClCC')]:
    #     if i not in [i1,i2,i3]:
    #         print(i)
    array1=[i1,i2,i3,i4]
    array2=[i1,i2,i3,i4]
    print(sorted(array1), 55555555555)
    # Probabilities Selections
    costs = np.array(array1)
    probs = np.exp((2*costs).astype('float64'))# / np.sum(np.exp(36*costs))
    print(costs[0])
