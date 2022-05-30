#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from rdkit import Chem
from rdkit.Chem import AllChem, rdmolops, DataStructs
from openbabel import openbabel as ob

from copy import deepcopy
import tempfile, subprocess, os, random, time
import numpy as np

#==================================================
# Class to work with lead
#==================================================

class Individual:

    def __init__(self,smiles:str = None, mol:Chem.rdchem.Mol = None, idx:int = 0, pdbqt = None,fragments:list = None, cost:float = np.inf) -> None:
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
        
        if fragments:
            self.fragments = fragments
        else:
            self.fragments = []
        
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

#==================================================
# Define some functions to work with
#==================================================
# Useful as decorator
def timeit(method):
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

def run(command, shell = True, executable = '/bin/bash', Popen = False):
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

def obconvert(inpath, outpath):
    """Convert  molecule ussing openbabel

    Args:
        input (str, path): input molecule.
        output (str, path): must have the extention of the molecule.
    """
    in_ext = os.path.basename(inpath).split('.')[-1]
    out_ext = os.path.basename(outpath).split('.')[-1]

    obConversion = ob.OBConversion()
    obConversion.SetInAndOutFormats(in_ext, out_ext)
    mol = ob.OBMol()
    obConversion.ReadFile(mol, inpath)   # Open Babel will uncompressed automatically
    #mol.AddHydrogens()
    obConversion.WriteFile(mol, outpath)

def confgen(smiles, outformat = "pdbqt"):
    """Create a 3D model from smile to

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

def fragments(mol):
    break_point = int(random.choice(np.where(np.array([b.GetBondType() for b in mol.GetBonds()]) == Chem.rdchem.BondType.SINGLE)[0]))
    # Chem.FragmentOnBonds(mol, break_point) # Could be used to increase randomness give more possible fragments and select two of them
    with Chem.RWMol(mol) as rwmol:
        b = rwmol.GetBondWithIdx(break_point)
        rwmol.RemoveBond(b.GetBeginAtomIdx(), b.GetEndAtomIdx())
    return rdmolops.GetMolFrags(rwmol, asMols = True)

def rdkit_numpy_convert(fp):
    # fp - list of binary fingerprints
    output = []
    for f in fp:
        arr = np.zeros((1,))
        DataStructs.ConvertToNumpyArray(f, arr)
        output.append(arr)
    return np.asarray(output)

def get_top(ms, model):
    # ms - list of molecules
    # model - sklearn model
    fps1 = [AllChem.GetMorganFingerprintAsBitVect(m, 2) for m in ms]
    x1 = rdkit_numpy_convert(fps1)
    pred = model.predict(x1)
    i = np.argmax(pred)
    return ms[i], pred[i]

def get_sim(ms, ref_fps):
    # ms - list of molecules
    # ref_fps - list of fingerprints of reference molecules
    output = []
    fps1 = [AllChem.GetMorganFingerprintAsBitVect(m, 2) for m in ms]
    for fp in fps1:
        v = DataStructs.BulkTanimotoSimilarity(fp, ref_fps)
        i = np.argmax(v)
        output.append([v[i], i])
    return output


if __name__ == '__main__':
    import pickle
    initial_smiles = 'COC(=O)C=1C=CC(=CC1)S(=O)(=O)N'
    i = Individual(initial_smiles)
    
    print(i.pdbqt)
