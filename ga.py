#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from audioop import avg
from copy import deepcopy
import tempfile, os
from rdkit import Chem
from rdkit.Chem import Descriptors
from crem.crem import mutate_mol, grow_mol, link_mols
from lead import utility, vina
import random, tqdm, shutil, itertools
import multiprocessing as mp
import numpy as np
import pickle as _pickle
import matplotlib.pyplot as plt

from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*') 


## Problem!!
# !!!!!!!!!!!!See the roulet_weel_selection becasue the average part is given me no desired results
# The mutation that I am doing right now is more a crossover but with the CReM library instead with the population itself.
# Code a more conservative mutation (close to the one in the paper of the SMILES) to perform small mutations.
# Implement in the cost fucntino possible restraints to similaroty (I really dont like this option
# Add synthetic accesibility prediciton in the cost fucntion
# Create a better cost fucntion that take into account as much as possible all the necesities of a drug-like molecule,

# We could also, instead a crossover two mutations with different levels. Becasue our 'mutate' is indeed some kind of crossover but with the CReM data base. So, what we could do is implement the mutation sof the paper that is at small scale (local optimazation): the possibilities are point mutations (atom-by-atom), deletion, add new atoms
# I think that I could control this specifying 
# For the best ligand we could predict the metabolites with sygma (apart of the module)
# Think about how to handle possibles Vina Crash
# Sometimes the pdbqt structure is not generated (probably problems in the convertion. In this cases the whole simulation crash ans should not be like this, this should rise some warning and continue discarting this structure)
# Apply filter for chemical elemnts to avoid crash on the vina function, in case that Vina crash we could use the predicted model
# Till now Vina never fails but could happen.
# Add more cpus for the generation of the conformers with RDKit
# The size of the ligands increase with the number of generations (if crossover is used even more)
# How to implement the rationality of where to grow, not just randomness. That could be the "crossover" operation, in fact the grow in a specific direction, based in (for example) the interaction network or the clashing avoid.
# I have to create a filter of atoms in order that vina doesn't fail because B and atoms like that Vina is not able to handle.
# Implement a continuation and automatic saving method for possible crashing and relaunch of the simulation.
class GA(object):
    
    def __init__(self, smiles, costfunc, crem_db_path, maxiter, popsize, beta = 0.001, pc =1, **costfunc_kwargs) -> None:
        self.InitIndividual = utility.Individual(smiles)
        self.costfunc = costfunc
        self.crem_db_path = crem_db_path
        self.pop = [self.InitIndividual]

        self.maxiter = maxiter
        self.popsize = popsize
        self.beta = beta
        self.costfunc_kwargs = costfunc_kwargs
        self.nc = int(np.round(pc*popsize/2)*2)
        self.exceptions = 0

        # Tracking parameters
        self.SawIndividuals = []
    @utility.timeit
    def __call__(self, njobs:int = 1, predictor_model:vina.VinaScoringPredictor = None):
        # Initialize Population
        GenInitStructs = list(
            mutate_mol(
                Chem.AddHs(self.InitIndividual.mol),
                self.crem_db_path,
                radius=3,
                min_size=0, max_size = 8,
                min_inc=-3, max_inc=3,
                #max_replacements=self.popsize + int(self.popsize/2), # I have to generate more structure
                return_mol= True,
                ncores = njobs*self.costfunc_kwargs['vina_cpus']
                )
            )
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

        self.pop = [self.InitIndividual]
        for i, item in enumerate(GenInitStructs):
            _, mol = item
            mol = Chem.RemoveHs(mol)
            # Here I have to keep record of the added fragment
            self.pop.append(utility.Individual(Chem.MolToSmiles(mol), mol, idx = i + 1))# 0 is the InitIndividual
        
        # Calculating cost of each individual (Doing Docking)
        vina_jobs = tempfile.TemporaryDirectory(prefix='vina')
        pool = mp.Pool(njobs)
            # Creating the arguments
        args_list = []
        # Make a copy of the self.costfunc_kwargs
        kwargs_copy = self.costfunc_kwargs.copy()
        kwargs_copy['wd'] = vina_jobs.name
        for Individual in self.pop:
            args_list.append((Individual, kwargs_copy))

        print(f'\n\nCreating the first population with {self.popsize} members:')
        self.pop = [Individual for Individual in tqdm.tqdm(pool.imap(self.__costfunc__, args_list), total=len(args_list))]
        pool.close()
        shutil.rmtree(vina_jobs.name)
        
        # Print some information of the initial population

        BestIndividualOfInitPopulation = min(self.pop, key = lambda x:x.cost)
        print(f"Initial Population: Best individual: {BestIndividualOfInitPopulation.smiles}. Best cost: {BestIndividualOfInitPopulation.cost}")
        # Getting the info of the first individual (Father/Mother) to print at the end how well performed the method
        # Because How the population was initialized and because we are using pool.imap (ordered). The Father/Mother is the first Individual of self.pop
        self.InitIndividual = deepcopy(self.pop[0])

        # Saving tracking variables
        for Individual in self.pop:
            if Individual.smiles not in [si.smiles for si in self.SawIndividuals]:
                self.SawIndividuals.append(Individual)


        # Creating the first model
        # Could be more pretty like creating a method that is make the model, y lo que se hace es que se gurda el modelo
        # Aqui se hace por primera vez pero para no repetir tanto codigo solo se llama a update model
        
        if predictor_model:
            self.Predictor = predictor_model
            print('\nUpdating the provided model:\n')
            self.Predictor.update(
                new_smiles_scoring = dict(((Individual.smiles,Individual.cost) if Individual.cost != np.inf else (Individual.smiles,9999) for Individual in self.SawIndividuals)),
                receptor = os.path.basename(self.costfunc_kwargs['receptor_path']).split('.')[0],
                boxcenter = self.costfunc_kwargs['boxcenter'],
                boxsize = self.costfunc_kwargs['boxsize'],
                exhaustiveness = self.costfunc_kwargs['exhaustiveness'],
            )
            self.Predictor(n_estimators=100, n_jobs=njobs*self.costfunc_kwargs['vina_cpus'], random_state=42, oob_score=True)
            print('Done!')
        else:
            print('\nCreating the first predicted model...')
            self.Predictor = vina.VinaScoringPredictor(
                smiles_scoring = dict(((Individual.smiles,Individual.cost) if Individual.cost != np.inf else (Individual.smiles,9999) for Individual in self.SawIndividuals)),
                receptor = os.path.basename(self.costfunc_kwargs['receptor_path']).split('.')[0],
                boxcenter = self.costfunc_kwargs['boxcenter'],
                boxsize = self.costfunc_kwargs['boxsize'],
                exhaustiveness = self.costfunc_kwargs['exhaustiveness'],
            )
            self.Predictor(n_estimators=100, n_jobs=njobs*self.costfunc_kwargs['vina_cpus'], random_state=42, oob_score=True)
            print('Done!')

        print(f'The model presents a oob_score = {self.Predictor.model.oob_score_}\n')
        
        # Best Cost of Iterations
        self.bestcost = np.empty(self.maxiter)
        self.avg_cost = np.empty(self.maxiter)
        
        # Main Loop
        for iter in range(self.maxiter):
            # Probabilities Selections
            costs = np.array([Individual.cost for Individual in self.pop])
            probs = np.exp(-self.beta*costs) / np.sum(np.exp(-self.beta*costs))

            popc = []
            for _ in range(self.nc//2):
                # Perform Roulette Wheel Selection
                p1 = self.pop[self.roulette_wheel_selection(probs)]
                p2 = self.pop[self.roulette_wheel_selection(probs)]
                

                # I have to think about this operations, and the parameters to controll them
                # Perform Crossover
                # Implement something that tells me if there are repeated Individuals and change them.
                c1, c2 = self.crossover(p1, p2, ncores=njobs*self.costfunc_kwargs['vina_cpus'])

                # Perform Mutation
                c1 = self.mutate(c1, ncores=njobs*self.costfunc_kwargs['vina_cpus'])
                c2 = self.mutate(c2, ncores=njobs*self.costfunc_kwargs['vina_cpus'])
                
                # Save offspring population
                # I will save only those offsprings that were not seen 
                if c1.smiles not in [Individual.smiles for Individual in self.SawIndividuals]: popc.append(c1)
                if c2.smiles not in [Individual.smiles for Individual in self.SawIndividuals]: popc.append(c2)

            if popc: # Only if there are new members
                # Calculating cost of each offspring individual (Doing Docking)
                vina_jobs = tempfile.TemporaryDirectory(prefix='vina')
                #os.makedirs('.vina_jobs', exist_ok=True)
                pool = mp.Pool(njobs)
                # Creating the arguments
                args_list = []
                # Make a copy of the self.costfunc_kwargs
                kwargs_copy = self.costfunc_kwargs.copy()
                kwargs_copy['wd'] = vina_jobs.name
                for (i, Individual) in enumerate(popc):
                    # Add idx label to each Individual
                    Individual.idx = i
                    # The problem here is that we are not being general for other possible Cost functions.
                    args_list.append((Individual,kwargs_copy))
                print(f'\nEvaluating generation {iter + 1}:')

                #!!!! Here I have to see if the smiles are in the self.saw_smiles in order to do not perform the docking and just assign the scoring function

                popc = [Individual for Individual in tqdm.tqdm(pool.imap(self.__costfunc__, args_list), total=len(args_list))]  
                pool.close()
                shutil.rmtree(vina_jobs.name)
                
            # Merge, Sort and Select
            # This could be improved. The problem is that the population could start to get the same individual, 
            # The diversity of the population could be controlled in this steep
            # 
            self.pop += popc
            self.pop = sorted(self.pop, key=lambda x: x.cost)
            self.pop = self.pop[:self.popsize]

            # Store Best Cost
            self.bestcost[iter] = self.pop[0].cost

            # Store Average cost
            self.avg_cost[iter] = np.mean(np.array([Individual.cost for Individual in self.pop]))

            # Saving tracking variables and getting new ones for the model update
            new_smiles_cost = dict()
            for Individual in popc:
                if Individual.smiles not in [sa.smiles for sa in self.SawIndividuals]:
                    # New variables
                    if Individual.cost == np.inf:
                        new_smiles_cost[Individual.smiles] = 9999
                    else:
                        new_smiles_cost[Individual.smiles] = Individual.cost
                    #Tracking variables
                    self.SawIndividuals.append(Individual)
            # Update the model
            print(f'Updating the current model with the information of generation {iter + 1}...')

            self.Predictor.update(
                new_smiles_scoring = new_smiles_cost.copy(),
                receptor = os.path.basename(self.costfunc_kwargs['receptor_path']).split('.')[0],
                boxcenter = self.costfunc_kwargs['boxcenter'],
                boxsize = self.costfunc_kwargs['boxsize'],
                exhaustiveness = self.costfunc_kwargs['exhaustiveness'],
            )
            self.Predictor(n_estimators=100, n_jobs=njobs*self.costfunc_kwargs['vina_cpus'], random_state=42, oob_score=True)          
            print('Done!')
            print(f'The updated model presents a oob_score = {self.Predictor.model.oob_score_}')
            

            # Show Iteration Information
            print(f"Generation {iter + 1}: Best individual: {self.pop[0].smiles}. Best Cost = {self.pop[0].cost}.\n")
            plt.scatter(iter, self.pop[0].cost)
        
        # Printing summary information
        print(f"\n{50*'=+'}\n")
        print(f'The simulation finished successfully after {self.maxiter} generations with a population of {self.popsize} individuals.')
        print(f"Initial Structure: {self.InitIndividual.smiles}. Initial Cost: {self.InitIndividual.cost}")
        print(f"Final Structure: {self.pop[0].smiles}. Final Cost: {self.pop[0].cost}")
        print(f"\n{50*'=+'}\n")

    def __costfunc__(self, args_list):
        Individual, kwargs = args_list
        #This is just to use the progress bar on pool.imap
        return self.costfunc(Individual, **kwargs)

    def crossover(self, Individual1, Individual2, ncores = 1, probability = 0, MaxRatioOfIncreaseInWt = 0.25):
        # here I have to select some randomness to perform or not the real crossover because I think that we could get far from the solution. It is just a guess.
        # How do I control the size of the new offspring? 
        # Performing a fragmentation in such a way that the offspring is the same in size
        # Here is where more additional information could be used. In order to orient the design of the new offspring. 
        # Then, I should control how perform the mutation  in such a way that we could keep or at least evaluate the offspring generated for crossover
        if random.random() < probability: # 50% of return the same individuals
            fragments1 = utility.fragments(Individual1.mol)
            fragments2 = utility.fragments(Individual2.mol)
            all_fragments = list(fragments1) + list(fragments2)
            
            # Initialize offspring smiles; cost
            offsprings = [
                    [None, np.inf],
                    [None, np.inf],
            ]
            for combination in itertools.combinations(all_fragments, 2):
                
                # Combine the molecules
                try:
                    possible_fragments_smiles = list(link_mols(*combination, db_name=self.crem_db_path, radius = 3, min_atoms=1, max_atoms=6, return_mol=False, ncores=ncores))                
                except:
                    # This is for debugging
                    sm1, sm2 = [Chem.MolToSmiles(c) for c in combination]
                    raise RuntimeError(f'These are the problematic SMILES: {sm1}, {sm2}')
                
                # Perform a filter based on weight. This control the size of the fragments. For now I will test 25 %. Think in the future work with the mols instead of smiles, I have to convert to mols too many times in this section of the code
                avg_wt  = 0.5*(Descriptors.ExactMolWt(Individual1.mol) + Descriptors.ExactMolWt(Individual1.mol))
                threshold_wt = (MaxRatioOfIncreaseInWt + 1) * avg_wt
                print(f'We had {len(possible_fragments_smiles)} possible fragments')
                possible_fragments_smiles = list(filter(lambda x: Descriptors.ExactMolWt(Chem.MolFromSmiles(x)) < threshold_wt, possible_fragments_smiles))
                print(f'After the weight filter we have {len(possible_fragments_smiles)} possible fragments')

                # In case that it was not possible to link the fragments
                if not possible_fragments_smiles:continue

                # Here comes the prediction with the model, and get the top two
                temp_offsprings = list(zip(possible_fragments_smiles, self.Predictor.predict(possible_fragments_smiles).tolist()))
                
                # Merge, Sort and Select
                offsprings = sorted(offsprings + temp_offsprings, key = lambda x:x[1])[:2]
            # Here I should check that exist offsprings (there not None values as smiles). For now I will assume that we always get at least two. See on the future
            return utility.Individual(smiles = offsprings[0][0]), utility.Individual(smiles = offsprings[1][0])
        else:
            return Individual1, Individual2    
    
    # Improve
    def mutate(self, Individual, ncores = 1):
        # See the option max_replacment
        # Or select the mutant based on some criterion
        # try:
            # Here i will pick the molecules based on the model.
        # El problema de seleccionar asi los compuestos es que siempre seleccionamos los mismos. Siempre se esta entrando la misma estructura y terminamos con una pobalcion redundante
        # Esto tengo que pensarlo mejor
        # new_mols = list(mutate_mol(Chem.AddHs(Individual.mol), self.crem_db_path, radius=3, min_size=1, max_size=8,min_inc=-3, max_inc=3, return_mol=True, ncores = ncores))
        # new_mols = [Chem.RemoveHs(i[1]) for i in new_mols]
        # best_mol, score = utility.get_top(new_mols + [Individual.mol], self.model)
        # smiles = Chem.MolToSmiles(best_mol)
        # mol = best_mol
        # print(score)
        # For now I am generating all the mutants and picking only one at random, this is very inefficient, should be better only generate one, but I am afraid that crem generate always the same or not generate any at all.
        # I think that what first crem does is randomly select on spot and find there all possible mutants. If this spot doesn't generate mutants, then you don't get nothing. But this is a supposition. 
        try:
            smiles, mol = random.choice(list(mutate_mol(Individual.mol, self.crem_db_path, radius=3, min_size=1, max_size=8,min_inc=-5, max_inc=3, return_mol=True, ncores = ncores)))
        except:
            print('The mutation did not work, we returned the same individual')
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
    



