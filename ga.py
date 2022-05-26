#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from rdkit import Chem
from crem.crem import mutate_mol, grow_mol, link_mols
import utility, vina, random, os, tqdm, shutil
import multiprocessing as mp
import glob


class GA(object):
    
    def __init__(self, smiles, maxiter, popsize, crem_db_path, costfunc, **costfunc_keywords) -> None:
        self.InitIndividual = utility.Individual(smiles)
        self.pop = [self.InitIndividual]
        self.maxiter = maxiter
        self.popsize = popsize
        self.crem_db_path = crem_db_path
        self.costfunc = costfunc
        self.costfunc_keywords = costfunc_keywords
    
    def __call__(self, njobs = 1):
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
        for i, item in enumerate(GenInitStructs):
            _, mol = item
            mol = Chem.RemoveHs(mol)
            # Here I have to keep record of the added fragment
            self.pop.append(utility.Individual(Chem.MolToSmiles(mol), mol, idx = i + 1))# 0 is the InitIndividual
        
        # Calculating cost of each individual (Doing Docking)
        os.makedirs('.vina', exist_ok=True)
        pool = mp.Pool(njobs)
        # Creating the arguments
        args = []
        for Individual in self.pop:
            args.append(
                (
                Individual.idx,
                Individual.pdbqt,
                self.costfunc_keywords['receptor_path'],
                self.costfunc_keywords['boxcenter'],
                self.costfunc_keywords['boxsize'],
                self.costfunc_keywords['exhaustiveness'],
                self.costfunc_keywords['vina_cpus'],
                self.costfunc_keywords['num_modes'],
                '.vina'
                )
            )

        pdbqts_costs = []
        for pdbqt_cost in tqdm.tqdm(pool.imap(self.costfunc, args), total=len(args)):
            pdbqts_costs.append(pdbqt_cost)
        pool.close()
        # There are some problems with imap and the object that I create, so I have to assign the values postprocessing. 
        for i in range(len(self.pop)):
            self.pop[i].pdbqt, self.pop[i].cost = pdbqts_costs[i] # Have the same order therefore I can use the same index i
        shutil.rmtree('.vina',)

    def crossover(self):
        pass
    
    def mutate(self, mol):
        mols = list(mutate_mol(mol, self.crem_db_path, return_mol=True, min_size=1, max_size=1))
        return mols




if __name__ == '__main__':
    pass
    



