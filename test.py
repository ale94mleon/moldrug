#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from lead import ga, fitness
import json
from multiprocessing import cpu_count
import os
file_path = os.path.dirname(os.path.realpath(__file__))

receptor = '7e27'#7e27'#'6lu7'#'x0161'#'7e27_periplasm'
maxiter = 50
popsize = 50
njobs = 3
NumbCalls = 1

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
        'min_size':2,
        'max_size':5,
        'min_inc':-2,
        'max_inc':2,
        'ncores':cpu_count(),
    },
    costfunc = fitness.__CostSimilarity,# __VinaCostLipinski, Cost, __VinaCost, __QedSasVinaCost
    costfunc_kwargs = {
        'receptor_path': f'/home/ale/GITLAB/lead/data/{receptor}.pdbqt',
        'boxcenter' : grid_opt['boxcenter'],
        'boxsize': grid_opt['boxsize'],
        'exhaustiveness': 8,
        'ncores': int(cpu_count() / njobs),
        'num_modes': 1,
        'ref_smiles': init_smiles,
    },
    save_pop_every_gen = 2,
    pop_file_name = f'pkl/pop',
    )
for i in range(NumbCalls):
    out(njobs = njobs)

for o in out.pop:
    print(o.smiles, o.cost)
out.pickle(f'/home/ale/GITLAB/lead/pkl/desirability_similarity_{receptor}_NumGen_{out.NumGen}_PopSize_{popsize}', compress=True)
print(out.to_dataframe())