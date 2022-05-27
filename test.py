#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from logging import raiseExceptions
import ga, vina, json, pickle

with open('data/box.json', 'r') as f:
    grid_opt = json.load(f)['7e27']['A']
initial_smiles = 'COc1ccc(C(=O)/C=C(\\O)C(F)(F)C(F)(F)F)c(O)c1'

out = ga.GA(
    smiles=initial_smiles,
    maxiter=3,
    popsize=4,
    crem_db_path = '/home/ale/GITLAB/bi_crem_database/replacements02_sc2.5.db',
    
    costfunc = vina.VinaCostStar,
    receptor_path ='data/7e27.pdbqt',
    boxcenter = grid_opt['boxcenter'],
    boxsize = grid_opt['boxsize'],
    exhaustiveness = 8,
    vina_cpus = 3,
    num_modes = 1,
    )  
out(njobs = 4)
for o in out.pop:
    print(o.smiles, o.cost)
out.pickle('out.pkl')