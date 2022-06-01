#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from lead import ga, vina
import json

receptor = '7e27'#7e27'#'6lu7'#'x0161'
maxiter = 50
popsize = 50

with open('data/box.json', 'r') as f:
    grid_opt = json.load(f)[receptor]['A']
with open('data/smi.json', 'r') as f:
    init_smiles = json.load(f)[receptor]

out = ga.GA(
    smiles=init_smiles,
    maxiter=maxiter,
    popsize=popsize,
    crem_db_path = '/home/ale/GITLAB/bi_crem_database/replacements02_sc2.5.db',
    pc = 1,
    costfunc = vina.VinaCost,
    receptor_path =f'data/{receptor}.pdbqt',
    boxcenter = grid_opt['boxcenter'],
    boxsize = grid_opt['boxsize'],
    exhaustiveness = 8,
    vina_cpus = 3,
    num_modes = 1,
    apply_lipinski_filter = True,
    )  
out(njobs = 4)
for o in out.pop:
    print(o.smiles, o.cost)
out.pickle(f'lipinski_{receptor}_{maxiter}_{popsize}.pkl')