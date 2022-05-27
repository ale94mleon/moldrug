#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import tempfile
from rdkit import Chem
from crem.crem import mutate_mol, grow_mol, link_mols
import utility, random, tqdm, shutil, itertools
import multiprocessing as mp
import numpy as np
import matplotlib.pyplot as plt
import pickle as _pickle

import warnings
warnings.filterwarnings("ignore", message='not removing hydrogen atom with dummy atom neighbors')


## Problem!!
# Apply a filter for redundant molecules.
# The size of the ligands increase with the number of generations (if crossover is used even more)
# How to implement the rationality of where to grow, not just randomness. That could be the "crossover" operation, in fact the grow in a specific direction, based in (for example) the interaction network or the clashing avoid.
# Catch repeated structures. This will help to the total simulation time!!!!!!!
# I have to create a filter of atoms in order that vina doesn't fail because B and atoms like that Vina is not able to handle.

class GA(object):
    
    def __init__(self, smiles, costfunc, crem_db_path, maxiter, popsize, beta = 0.001, pc =1, **costfunc_keywords) -> None:
        self.InitIndividual = utility.Individual(smiles)
        self.costfunc = costfunc
        self.crem_db_path = crem_db_path
        self.pop = [self.InitIndividual]

        self.maxiter = maxiter
        self.popsize = popsize
        self.beta = beta
        self.costfunc_keywords = costfunc_keywords
        self.nc = int(np.round(pc*popsize/2)*2)
        self.exceptions = 0
    
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
        vina_jobs = tempfile.TemporaryDirectory(prefix='vina')
        #os.makedirs('.vina_jobs', exist_ok=True)
        pool = mp.Pool(njobs)
        # Creating the arguments
        args = []
        for Individual in self.pop:
            args.append(
                (
                Individual,
                self.costfunc_keywords['receptor_path'],
                self.costfunc_keywords['boxcenter'],
                self.costfunc_keywords['boxsize'],
                self.costfunc_keywords['exhaustiveness'],
                self.costfunc_keywords['vina_cpus'],
                self.costfunc_keywords['num_modes'],
                vina_jobs.name
                )
            )
        print(f'Creating the first population with {self.popsize} members:')
        self.pop = [Individual for Individual in tqdm.tqdm(pool.imap(self.costfunc, args), total=len(args))]
        print('Done!')
        pool.close()
        shutil.rmtree(vina_jobs.name)

        # Best Cost of Iterations
        self.bestcost = np.empty(self.maxiter)
        # Main Loop
        for iter in range(self.maxiter):
            costs = np.array([Individual.cost for Individual in self.pop])
            avg_cost = np.mean(costs)
            if avg_cost != 0:
                costs = costs/avg_cost
            probs = np.exp(-self.beta*costs)

            popc = []
            for _ in range(self.nc//2):
                # Perform Roulette Wheel Selection
                p1 = self.pop[self.roulette_wheel_selection(probs)]
                p2 = self.pop[self.roulette_wheel_selection(probs)]
                

                # I have to think about this operations, and the parameters to controll them
                # Perform Crossover
                # Implement something that tells me if there are repeated Individuals and change them.
                c1, c2 = self.crossover(p1, p2, ncores=njobs*self.costfunc_keywords['vina_cpus'])

                # Perform Mutation
                c1 = self.mutate(c1, ncores=njobs*self.costfunc_keywords['vina_cpus'])
                c2 = self.mutate(c2, ncores=njobs*self.costfunc_keywords['vina_cpus'])

                # Save offspring population
                popc.append(c1)
                popc.append(c2)

            # Calculating cost of each offspring individual (Doing Docking)
            vina_jobs = tempfile.TemporaryDirectory(prefix='vina')
            #os.makedirs('.vina_jobs', exist_ok=True)
            pool = mp.Pool(njobs)
            # Creating the arguments
            args = []
            for (i, Individual) in enumerate(popc):
                # Add idx label to each Individual
                Individual.idx = i
                args.append(
                    (
                    Individual,
                    self.costfunc_keywords['receptor_path'],
                    self.costfunc_keywords['boxcenter'],
                    self.costfunc_keywords['boxsize'],
                    self.costfunc_keywords['exhaustiveness'],
                    self.costfunc_keywords['vina_cpus'],
                    self.costfunc_keywords['num_modes'],
                    vina_jobs.name
                    )
                )
            print(f'Evaluating generation {iter}:')
            popc = [Individual for Individual in tqdm.tqdm(pool.imap(self.costfunc, args), total=len(args))]  
            pool.close()
            shutil.rmtree(vina_jobs.name)

                
            # Merge, Sort and Select
            self.pop += popc
            self.pop = sorted(self.pop, key=lambda x: x.cost)
            self.pop = self.pop[0:self.popsize]

            # Store Best Cost
            self.bestcost[iter] = self.pop[0].cost

            # Show Iteration Information
            print(f"Generation {iter}: Best Cost = {self.pop[0].cost}; Best individual: {self.pop[0].smiles}")
            plt.scatter(iter, self.pop[0].cost)

   # Improve
    def crossover(self, Individual1, Individual2, ncores = 1):
        # here I have to select some randomness to perform or not the real crossover because I think that we could get far from the solution. It is just a guess.
        if random.choice([False, False]): # Testing without crossover
            fragments1 = utility.fragments(Individual1.mol)
            fragments2 = utility.fragments(Individual2.mol)
            all_fragments = list(fragments1) + list(fragments2)
            combination1, combination2 = list(itertools.combinations(all_fragments, 2))[:2]
            try:
                offspring1_smile, offspring1_mol = random.choice(list(link_mols(*combination1, db_name=self.crem_db_path, radius = 3, max_replacements = 5, return_mol=True, ncores=ncores)))
            except:
                self.exceptions += 1
                offspring1_smile, offspring1_mol = Individual1.smiles, Individual1.mol
            try:
                offspring2_smile, offspring2_mol = random.choice(list(link_mols(*combination2, db_name=self.crem_db_path, radius = 3, max_replacements = 5, return_mol=True, ncores=ncores)))
            except:
                self.exceptions += 1
                offspring2_smile, offspring2_mol = Individual2.smiles, Individual2.mol
            return utility.Individual(offspring1_smile, offspring1_mol), utility.Individual(offspring2_smile, offspring2_mol)
        else:
            return Individual1, Individual2
    
    # Improve
    def mutate(self, Individual, ncores = 1):
        # See the option max_replacment
        # Or select the mutant based on some criterion
        try:
            smiles, mol = random.choice(list(mutate_mol(Individual.mol, self.crem_db_path, radius=3, min_size=1, max_size=8, max_replacements = 5, return_mol=True, ncores = ncores)))
        except:
            self.exceptions += 1
            smiles, mol = Individual.smiles, Individual.mol
        return utility.Individual(smiles,mol)
    
    def roulette_wheel_selection(self, p):
        c = np.cumsum(p)
        r = sum(p)*np.random.rand()
        ind = np.argwhere(r <= c)
        return ind[0][0]
    
    def pickle(self,file):
        cls = self.__class__
        result = cls.__new__(cls)
        result.__dict__.update(self.__dict__)
        with open(file, 'wb') as pkl:
            _pickle.dump(result, pkl)


if __name__ == '__main__':
    pass
    



