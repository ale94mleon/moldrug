#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from rdkit import Chem
from crem.crem import mutate_mol, grow_mol, link_mols

import numpy as np
import random
import copy as _copy

import subprocess

# Define struct class
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
    def __init__(self,smiles:str = None, mol:Chem.rdchem.Mol = None, fragments:list = None, cost:float = np.inf) -> None:
        self.smiles = smiles
        if not mol:
            try:
                self.mol = Chem.MolFromSmiles(smiles)
            except:
                self.mol = None
        else:
            self.mol = mol
        
        self.cost = cost
        if fragments:
            self.fragments = fragments
        else:
            self.fragments = []

class GA(object):
    
    def __init__(self, smiles, costfunc, maxiter, popsize, crem_db_path) -> None:
        self.InitIndividual = Individual(smiles)
        self.pop = [self.InitIndividual]
        self.costfunc = costfunc
        self.maxiter = maxiter
        self.popsize = popsize
        self.crem_db_path = crem_db_path
    
    def __call__(self):
        # Initialize Population
        GenInitStructs = list(
            mutate_mol(
                Chem.AddHs(self.InitIndividual.mol),
                self.crem_db_path,
                radius=3,
                min_size=0, max_size = 8,
                max_replacements=self.popsize + int(self.popsize/2), # I have to generate more structure
                return_mol= True,
                )
            )

        if len(GenInitStructs) < (self.popsize - 1):
            print('The initial population has repeated elements')
            # temporal solution
            GenInitStructs +=  random.choices(GenInitStructs, k = self.popsize - len(GenInitStructs) -1)
            pass# I am not sure how to deal with this
        elif len(GenInitStructs) > (self.popsize - 1):
            #Selected random sample from the generation 
            GenInitStructs = GenInitStructs = random.sample(GenInitStructs, k = self.popsize -1)
        else:
            # Everything is ok!
            pass 

        self.pop = [self.InitIndividual]
        for smiles, mol in GenInitStructs:
            mol = Chem.RemoveHs(mol)
            # Here I have to keep record of the added fragment
            self.pop.append(Individual(Chem.MolToSmiles(mol), mol))


    
    def crossover(self):
        pass
    
    def mutate(self):
        mols = list(mutate_mol(mol, db_fname, return_mol=True, min_size=1, max_size=1))
        pass

#==================================================
# Define cost functions. for now lets focus on Vina
#==================================================
def VinaCost(smiles):
    return np.random.rand()


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


if __name__ == '__main__':
    ga = GA('O=C(C)Oc1ccccc1C(=O)O', VinaCost, 20, 5, crem_db_path = '/home/ale/GITLAB/bi_crem_database/replacements02_sc2.5.db')
    out = ga()
    # i = Individual('O=C(C)Oc1ccccc1C(=O)O')
    # print(i)
    



