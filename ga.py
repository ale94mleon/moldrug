#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from copy import deepcopy
from rdkit import Chem
import numpy as np

class EmptyIndividual(object):
    def __init__(self,smile:str = None, cost:float = np.inf, fragments:list = None) -> None:
        self.smile = smile
        try:
            self.mol = Chem.MolFromSmile(smile)
        except:
            self.mol = None
        
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

class GA(object):
    def __init__(self, costfunc, maxiter, popsize) -> None:
        self.costfunc = costfunc
        self.maxiter = maxiter
        self.popsize = popsize
    def __call__(self, *args: Any, **kwds: Any) -> Any:
        pass
    def crossover(self):
        pass
    def mutate(self):
        pass

# Define cost functions. for now lets focus on Vina


