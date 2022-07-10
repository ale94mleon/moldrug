#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from rdkit import Chem
from moldrug import utils, fitness
from moldrug.data import receptors, ligands, boxes
import tempfile, os, gzip, shutil, requests
from multiprocessing import cpu_count

# Creating a temporal directory
tmp_path = tempfile.TemporaryDirectory()
# Creating receptors files
r_x0161_file = os.path.join(tmp_path.name, 'r_x0161.pdbqt')
r_6lu7_file = os.path.join(tmp_path.name, 'r_6lu7.pdbqt')
with open(r_x0161_file, 'w') as r:
    r.write(receptors.r_x0161)
with open(r_6lu7_file, 'w') as r:
    r.write(receptors.r_6lu7)


# Getting the crem data base
url = "http://www.qsar4u.com/files/cremdb/replacements02_sc2.db.gz"
r = requests.get(url, allow_redirects=True)
crem_dbgz_path = os.path.join(tmp_path.name,'crem.db.gz')
crem_db_path = os.path.join(tmp_path.name,'crem.db')
open(crem_dbgz_path, 'wb').write(r.content)
with gzip.open(crem_dbgz_path, 'rb') as f_in:
    with open(crem_db_path, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)




def test_single_receptor():
    maxiter = 1
    popsize = 2
    njobs = 2
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
    
    out = utils.GA(
        seed_smiles=ligands.r_x0161,
        maxiter=maxiter,
        popsize=popsize,
        crem_db_path = crem_db_path,#'/home/ale/GITLAB/bi_crem_database/replacements02_sc2.db',
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
            'vina_executable': 'vina',
            'receptor_path': r_x0161_file,
            'boxcenter' : boxes.r_x0161["A"]['boxcenter'] ,
            'boxsize': boxes.r_x0161["A"]['boxsize'],
            'exhaustiveness': 4,
            'ncores': int(cpu_count() / njobs),
            'num_modes': 1,
        },
        save_pop_every_gen = 20,
        pop_file_name = os.path.join(tmp_path.name, 'pop_test_single_receptor')
        )

    for _ in range(NumbCalls):
        out(njobs = njobs)

    for o in out.pop:
        print(o.smiles, o.cost)
    out.pickle(os.path.join(tmp_path.name, f"result_test_single_receptor_NumGen_{out.NumGen}_PopSize_{popsize}"), compress=True)
    print(out.to_dataframe())

def test_multi_receptor():
    maxiter = 1
    popsize = 2
    njobs = 3
    NumbCalls = 1

    out = utils.GA(
        seed_smiles=ligands.r_x0161,
        maxiter=maxiter,
        popsize=popsize,
        crem_db_path = crem_db_path,#'/home/ale/GITLAB/bi_crem_database/replacements02_sc2.db',
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
            'receptor_path': [r_x0161_file, r_6lu7_file],
            'boxcenter' : [boxes.r_x0161["A"]['boxcenter'], boxes.r_6lu7["A"]['boxcenter']],
            'boxsize': [boxes.r_x0161["A"]['boxsize'], boxes.r_6lu7["A"]['boxsize']],
            'vina_score_types': ['min', 'max'],
            'exhaustiveness': 4,
            'ncores': int(cpu_count() / njobs),
            'num_modes': 1,
        },
        save_pop_every_gen = 20,
        pop_file_name = os.path.join(tmp_path.name, 'pop_test_multi_receptor')
        )

    for _ in range(NumbCalls):
        out(njobs = njobs)

    for o in out.pop:
        print(o.smiles, o.cost)
    out.pickle(os.path.join(tmp_path.name, f"result_test_multi_receptor_NumGen_{out.NumGen}_PopSize_{popsize}"), compress=True)
    print(out.to_dataframe())

def test_local():
    njobs = 2
    local = utils.Local(
        mol = Chem.AddHs(Chem.MolFromSmiles(ligands.r_x0161)),
        crem_db_path = crem_db_path,
        grow_crem_kwargs = {
                'radius':3,
                'min_atoms':1,
                'max_atoms':3,
                'ncores':cpu_count(),
            },
        costfunc = fitness.Cost,
        costfunc_kwargs = {
            'receptor_path': r_x0161_file,
            'boxcenter' : boxes.r_x0161["A"]['boxcenter'],
            'boxsize': boxes.r_x0161["A"]['boxsize'],
            'exhaustiveness': 8,
            'ncores': int(cpu_count() / njobs),
            'num_modes': 1,
        },
    )
    local(njobs = njobs, pick=2)
    print(local.to_dataframe())
