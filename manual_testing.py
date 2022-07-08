#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from moldrug import fitness, utils
import json
from multiprocessing import cpu_count
import os

from rdkit import Chem


file_path = os.path.dirname(os.path.realpath(__file__))

TypeOfTest = [
    'single_receptor',
    'multi_receptor',
    'local',
]
TypeOfTest = TypeOfTest[0]

if TypeOfTest == 'single_receptor':

    receptor = '7e27'#7e27'#'6lu7'#'x0161'#'7e27_periplasm'
    maxiter = 2
    popsize = 3
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
        box = json.load(f)
    with open(os.path.join(file_path,'data/smi.json'), 'r') as f:
        init_smiles = json.load(f)[receptor]

    out = utils.GA(
        seed_smiles=init_smiles,
        maxiter=maxiter,
        popsize=popsize,
        crem_db_path = '/home/ale/GITLAB/bi_crem_database/replacements02_sc2.5.db',
        pc = 1,
        get_similar = False,
        mutate_crem_kwargs = {
            'radius':3,
            'min_size':1,
            'max_size':8,
            'min_inc':-5,
            'max_inc':3,
            'ncores':cpu_count(),
        },
        costfunc = fitness.Cost,#__CostSimilarity,# __VinaCostLipinski, Cost, __VinaCost, __QedSasVinaCost, CostMultiReceptors
        costfunc_kwargs = {
            'receptor_path': f'/home/ale/GITLAB/moldrug/data/{receptor}.pdbqt' ,
            'boxcenter' : box[receptor]['A']['boxcenter'] ,
            'boxsize': box[receptor]['A']['boxsize'],
            'exhaustiveness': 8,
            'ncores': int(cpu_count() / njobs),
            'num_modes': 1,
            #'ref_smiles': init_smiles,
        },
        save_pop_every_gen = 20,
        pop_file_name = f'/home/ale/GITLAB/moldrug/pkl/pop',
        )

    for i in range(NumbCalls):
        out(njobs = njobs)

    for o in out.pop:
        print(o.smiles, o.cost)
    out.pickle(f"/home/ale/GITLAB/moldrug/pkl/desirability_{receptor}_NumGen_{out.NumGen}_PopSize_{popsize}", compress=True)
    print(out.to_dataframe())

elif TypeOfTest == 'multi_receptor':
    receptor = ['7e27_periplasm', '7e27']#7e27'#'6lu7'#'x0161'#'7e27_periplasm'
    maxiter = 2
    popsize = 3
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
        box = json.load(f)
    with open(os.path.join(file_path,'data/smi.json'), 'r') as f:
        init_smiles = json.load(f)[receptor[0]]

    out = utils.GA(
        seed_smiles=init_smiles,
        maxiter=maxiter,
        popsize=popsize,
        crem_db_path = '/home/ale/GITLAB/bi_crem_database/replacements02_sc2.5.db',
        pc = 1,
        get_similar = False,
        mutate_crem_kwargs = {
            'radius':3,
            'min_size':1,
            'max_size':8,
            'min_inc':-5,
            'max_inc':3,
            'ncores':cpu_count(),
        },
        costfunc = fitness.CostMultiReceptors,#__CostSimilarity,# __VinaCostLipinski, Cost, __VinaCost, __QedSasVinaCost, CostMultiReceptors
        costfunc_kwargs = {
            'receptor_path': [f'/home/ale/GITLAB/moldrug/data/{r}.pdbqt' for r in receptor],
            'boxcenter' : [box[r]['A']['boxcenter'] for r in receptor],
            'boxsize': [box[r]['A']['boxsize'] for r in receptor],
            'vina_score_types': ['min', 'min'],
            'exhaustiveness': 8,
            'ncores': int(cpu_count() / njobs),
            'num_modes': 1,
            #'ref_smiles': init_smiles,
        },
        save_pop_every_gen = 20,
        pop_file_name = f'/home/ale/GITLAB/moldrug/pkl/pop',
        )

    for i in range(NumbCalls):
        out(njobs = njobs)

    for o in out.pop:
        print(o.smiles, o.cost)
    out.pickle(f"/home/ale/GITLAB/moldrug/pkl/desirability_{'_'.join([r for r in receptor])}_NumGen_{out.NumGen}_PopSize_{popsize}", compress=True)
    print(out.to_dataframe())

elif TypeOfTest == 'local':
    receptor = '7e27'#7e27'#'6lu7'#'x0161'#'7e27_periplasm'
    njobs = 3

    with open(os.path.join(file_path,'data/box.json'), 'r') as f:
        box = json.load(f)
    with open(os.path.join(file_path,'data/smi.json'), 'r') as f:
        init_smiles = json.load(f)[receptor]


    # I have to see how is the position if the refer to the heavy atoms or the hydrogens
    mol = Chem.AddHs(Chem.MolFromSmiles(init_smiles))
    local = utils.Local(
        mol = mol,
        crem_db_path = '/home/ale/GITLAB/bi_crem_database/replacements02_sc2.5.db',
        grow_crem_kwargs = {
                'radius':3,
                'min_atoms':1,
                'max_atoms':3,
                'ncores':cpu_count(),
            },
            costfunc = fitness.Cost,#__CostSimilarity,# __VinaCostLipinski, Cost, __VinaCost, __QedSasVinaCost, CostMultiReceptors
            costfunc_kwargs = {
                'receptor_path': f'/home/ale/GITLAB/moldrug/data/{receptor}.pdbqt' ,
                'boxcenter' : box[receptor]['A']['boxcenter'] ,
                'boxsize': box[receptor]['A']['boxsize'],
                'exhaustiveness': 8,
                'ncores': int(cpu_count() / njobs),
                'num_modes': 1,
            },
    )
    local(njobs = njobs, pick=5)
    print(local.to_dataframe())
