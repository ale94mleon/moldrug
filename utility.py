#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from rdkit import Chem
from rdkit.Chem import AllChem
from openbabel import openbabel as ob

import copy as _copy
import tempfile, subprocess, os
import numpy as np
from sklearn.multiclass import OutputCodeClassifier

#==================================================
# Class to work with lead
#==================================================
class struct (dict):
    """
    A class to implement the C/C++ or MATLAB-like structures
    Taken from ypstruct package.
    """
    
    def __repr__(self):
        """
        String representation of the struct
        """
        return "struct({})".format(super().__repr__())
    

    def __getattr__(self, field):
        """
        Gets value of a field
        """
        if field not in dir(self):
            if field in self.keys():
                return self[field]
            else:
                return None
        else:
            return None
    
    
    def __setattr__(self, field, value):
        """
        Sets value of a field
        """
        if field not in dir(self):
            self[field] = value
        else:
            return super().__setattr__(field, value)
    
    
    def fields(self):
        """
        Gets the list of defined fields of the struct
        """
        return list(self.keys())

    
    def remove_field(self, field):
        """
        Removes a field from the struct
        """
        if field in self.keys():
            del self[field]
    
    
    def add_field(self, field, value = None):
        """
        Adds a new field to the struct
        """
        if field not in self.keys():
            self[field] = value

    
    def copy(self):
        """
        Creates a shallow copy of the struct
        """
        self_copy = struct()
        for field in self.keys():
            if isinstance(self[field], struct):
                self_copy[field] = self[field].copy()
            else:
                self_copy[field] = _copy.copy(self[field])
        
        return self_copy

    
    def deepcopy(self):
        """
        Creates a deep copy of the struct
        """
        self_copy = struct()
        for field in self.keys():
            if isinstance(self[field], struct):
                self_copy[field] = self[field].deepcopy()
            else:
                self_copy[field] = _copy.deepcopy(self[field])
        
        return self_copy

    
    def repeat(self, n):
        """
        Repeats/replicates the struct to create an array of structs (eg. for initialization)
        """
        return [self.deepcopy() for _ in range(n)]

    
    def __mul__(self, n):
        """
        Overload * operator (multiplication) to repeat/replicate the struct
        """
        if not isinstance(n, int) and not isinstance(n, float):
            raise TypeError("Only integers are allowed.")
        return self.repeat(n)

    
    def __add__(self, other):
        """
        Overload + operator (addition) to merge two struct objects
        """
        if not isinstance(other, dict):
            raise TypeError("Only structure and dict objects are allowed.")
        result = self.deepcopy()
        result.update(other)
        return result

class Individual(struct):
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

#==================================================
# Define some functions to work with
#==================================================
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



if __name__ == '__main__':
    initial_smiles = 'COC1CC[C@@H](C(O)CC(O)C(F)(F)C(F)(F)F)C(O)C1'
    s = confgen(initial_smiles, 'pdb')
    with open('test.pdb', 'w') as mol:
        mol.write(s)