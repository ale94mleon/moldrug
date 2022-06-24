#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from lead import ga, fitness
import json
from multiprocessing import cpu_count
import os
file_path = os.path.dirname(os.path.realpath(__file__))

receptor = '7e27'#7e27'#'6lu7'#'x0161'#'7e27_periplasm'
maxiter = 5
popsize = 5
njobs = 3
NumbCalls = 1

"""For a local optimization we could use
    min_size=1, max_size=1, min_inc=-1, max_inc=1
    Or use just grow instead
    min_size=0, max_size=0,
    And if we would like to integrate all this option
    min_size=0, max_size=1, min_inc=-1, max_inc=1
    this will add, delate heavy mutate heavy atoms.
    or change Hydrogens for heavy atoms

    with
    min_size=1, max_size=8, min_inc=-5, max_inc=3
    I explore efficiently the chemical space, not good for cost functions with similarity included. 
"""

with open(os.path.join(file_path,'data/box.json'), 'r') as f:
    grid_opt = json.load(f)[receptor]['A']
with open(os.path.join(file_path,'data/smi.json'), 'r') as f:
    init_smiles = json.load(f)[receptor]

out = ga.GA(
    smiles=init_smiles,
    maxiter=maxiter,
    popsize=popsize,
    crem_db_path = '/home/ale/GITLAB/bi_crem_database/replacements02_sc2.5.db',
    pc = 1,
    get_similar = True,
    mutate_crem_kwargs = {
        'radius':3,
        'min_size':0,
        'max_size':1,
        'min_inc':-1,
        'max_inc':1,
        'ncores':cpu_count(),
    },
    costfunc = fitness.Cost,#__CostSimilarity,# __VinaCostLipinski, Cost, __VinaCost, __QedSasVinaCost
    costfunc_kwargs = {
        'receptor_path': f'/home/ale/GITLAB/lead/data/{receptor}.pdbqt',
        'boxcenter' : grid_opt['boxcenter'],
        'boxsize': grid_opt['boxsize'],
        'exhaustiveness': 8,
        'ncores': int(cpu_count() / njobs),
        'num_modes': 1,
        #'ref_smiles': init_smiles,
    },
    save_pop_every_gen = 20,
    pop_file_name = f'/home/ale/GITLAB/lead/pkl/pop',
    )
for i in range(NumbCalls):
    out(njobs = njobs)

for o in out.pop:
    print(o.smiles, o.cost)
out.pickle(f'/home/ale/GITLAB/lead/pkl/desirability_local_{receptor}_NumGen_{out.NumGen}_PopSize_{popsize}', compress=True)
print(out.to_dataframe())