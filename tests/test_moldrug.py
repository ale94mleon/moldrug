#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from rdkit import Chem
from moldrug import utils, fitness, home
from moldrug.data import receptors, ligands, boxes
import tempfile, os, gzip, shutil, requests, yaml, copy, sys
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

def test_single_receptor_command_line():
    Config = {
        "01_grow": {
            "type": "GA",
            "njobs": 3,
            "seed_smiles": ligands.r_x0161,
            "costfunc": "Cost",
            "costfunc_kwargs": {
                "vina_executable": "vina",
                "receptor_path": r_x0161_file,
                "boxcenter": boxes.r_x0161["A"]['boxcenter'],
                "boxsize": boxes.r_x0161["A"]['boxsize'],
                "exhaustiveness": 4,
                "ncores": 4,
                "num_modes": 1
            },
            "crem_db_path": crem_db_path,
            "maxiter": 1,
            "popsize": 2,
            "beta": 0.001,
            "pc": 1,
            "get_similar": False,
            "mutate_crem_kwargs": {
                "radius": 3,
                "min_size": 0,
                "max_size": 0,
                "min_inc": -5,
                "max_inc": 6,
                "ncores": 1
            },
            "save_pop_every_gen": 10,
            "deffnm": "01_grow"
        },
        "02_allow_grow": {
            "mutate_crem_kwargs": {
                "radius": 3,
                "min_size": 1,
                "max_size": 0,
                "min_inc": -5,
                "max_inc": 3,
                "ncores": 12
            },
            "maxiter": 1,
            "deffnm": "02_allow_grow"
        }
    }
    cwd = os.getcwd()
    with open(os.path.join(tmp_path.name, "test_single_receptor.yml"), 'w') as c:
        yaml.dump(Config, c)
    os.chdir(tmp_path.name)
    utils.run('moldrug test_single_receptor.yml', Popen=True)
    os.chdir(cwd)


def test_multi_receptor(maxiter = 1, popsize = 2, njobs = 3, NumbCalls = 1):
    out = utils.GA(
        seed_smiles=ligands.r_x0161,
        maxiter=maxiter,
        popsize=popsize,
        crem_db_path = crem_db_path,
        pc = 1,
        get_similar = True,
        mutate_crem_kwargs = {
            'radius':3,
            'max_size':0,
            'min_inc':-5,
            'max_inc':3,
        },
        costfunc = fitness.CostMultiReceptors,
        costfunc_kwargs = {
            'receptor_paths': [r_x0161_file, r_6lu7_file],
            'boxcenters' : [boxes.r_x0161["A"]['boxcenter'], boxes.r_6lu7["A"]['boxcenter']],
            'boxsizes': [boxes.r_x0161["A"]['boxsize'], boxes.r_6lu7["A"]['boxsize']],
            'vina_score_types': ['min', 'max'],
            'exhaustiveness': 4,
            'ncores': int(cpu_count() / njobs),
            'num_modes': 1,
        },
        save_pop_every_gen = 20,
        deffnm = os.path.join(tmp_path.name, 'test_multi_receptor')
        )

    for _ in range(NumbCalls):
        out(njobs = njobs)

    for o in out.pop:
        print(o.smiles, o.cost)
    out.pickle(os.path.join(tmp_path.name, f"result_test_multi_receptor_NumGens_{out.NumGens}_PopSize_{popsize}"), compress=False)
    print(out.to_dataframe())

    with open(os.path.join(tmp_path.name, 'vina_out.pdbqt'), 'w') as p:
        p.write(out.pop[0].pdbqts[0])

    vina_out = utils.VINA_OUT(os.path.join(tmp_path.name, 'vina_out.pdbqt'))
    vina_out.chunks[0].get_atoms()
    vina_out.chunks[0].write(os.path.join(tmp_path.name, 'chunk.pdbqt'))
    cwd = os.getcwd()
    os.chdir(tmp_path.name)
    vina_out.chunks[0].write()
    os.chdir(cwd)


def test_local_command_line():
    Config = {
        "main": {
            "type": "Local",
            "njobs": 1,
            "pick": 2,
            "mol": Chem.MolToSmiles(Chem.AddHs(Chem.MolFromSmiles(ligands.r_x0161))),
            "costfunc": "CostOnlyVina",
            "costfunc_kwargs": {
                "vina_executable": "vina",
                "receptor_path": r_x0161_file,
                "boxcenter": boxes.r_x0161["A"]['boxcenter'],
                "boxsize": boxes.r_x0161["A"]['boxsize'],
                "exhaustiveness": 4,
                "num_modes": 1
            },
            "crem_db_path": crem_db_path,
            "grow_crem_kwargs": {
                "radius": 3,
                "min_atoms": 1,
                "max_atoms": 3
            }
        }
    }
    cwd = os.getcwd()
    os.chdir(tmp_path.name)
    with open("local_config.yml", 'w') as c:
        yaml.dump(Config, c)

    utils.run(f"moldrug local_config.yml --fitness {os.path.join(home.home(), 'fitness.py')} --outdir results")
    print(os.listdir())
    # This problems with the modules are not so convenient
    os.chdir('results')
    sys.path.append('.')
    result = utils.decompress_pickle('local_result.pbz2')
    result.pickle('local_non_compress', compress=False)
    print(result.to_dataframe())
    os.chdir(cwd)


def test_CostOnlyVina():
    ligand_smiles = ligands.r_x0161
    I = utils.Individual(ligand_smiles)
    receptor_paths = [r_x0161_file,r_6lu7_file]
    boxcenters = [boxes.r_x0161['A']['boxcenter'], boxes.r_6lu7['A']['boxcenter']]
    boxsizes = [boxes.r_x0161['A']['boxsize'], boxes.r_6lu7['A']['boxsize']]
    vina_score_types = ['min', 'max']

    fitness.CostMultiReceptorsOnlyVina(Individual = I,wd = tmp_path.name,receptor_paths = receptor_paths, vina_score_types = vina_score_types, boxcenters = boxcenters,boxsizes = boxsizes,exhaustiveness = 4,ncores = 4)
    fitness.CostMultiReceptorsOnlyVina(Individual = I,wd = tmp_path.name,receptor_paths = receptor_paths, vina_score_types = vina_score_types, boxcenters = boxcenters,boxsizes = boxsizes,exhaustiveness = 4,ncores = 4, wt_cutoff=2)

    fitness.CostOnlyVina(Individual = I,wd = tmp_path.name,receptor_path = r_x0161_file, boxcenter = boxes.r_x0161['A']['boxcenter'], boxsize = boxes.r_x0161['A']['boxsize'],exhaustiveness = 4,ncores = 4)
    fitness.CostOnlyVina(Individual = I,wd = tmp_path.name,receptor_path = r_x0161_file, boxcenter = boxes.r_x0161['A']['boxcenter'], boxsize = boxes.r_x0161['A']['boxsize'],exhaustiveness = 4,ncores = 4, wt_cutoff=2)



def test_home():
    home.home(dataDir='data')


def test_get_sim_utils():
    from rdkit.Chem import AllChem
    mols = []
    ref_fps = []
    for s in ['CC', 'CCO']:
        mol = Chem.MolFromSmiles(s)
        mols.append(mol)
        ref_fps.append(AllChem.GetMorganFingerprintAsBitVect(mol, 2))
    utils.get_sim(mols, ref_fps)


def test_lipinski():
    mol = Chem.MolFromSmiles('CCCO')
    utils.lipinski_filter(Chem.MolFromSmiles('BrCC(COCN)CC(Br)CC(Cl)CCc1ccccc1CCCC(NCCCO)'))
    utils.lipinski_filter(mol)
    utils.lipinski_profile(mol)


def test_Individual():
    I1 = utils.Individual('CC', cost=10)
    I2 = utils.Individual('CCO', cost=2)
    I3 = copy.copy(I1)
    I4 = copy.deepcopy(I2)
    I5 = utils.Individual('CC', cost = 10, pdbqt=I1.pdbqt)
    assert I3 + I4 == 12
    assert I5 - I2 == 8
    assert I1 > I2
    assert I1 >= I2
    assert I2 < I1
    assert I2 <= I1
    assert -I1 == -10
    assert abs(I1) == 10
    assert I1 * I2 == 20
    assert I1 / I2 == 5
    assert I1 // I2 == 5
    assert I1 % I2 == 0
    assert divmod(I1, I2) == (5,0)
    assert I1**I2 == 100



def test_miscellanea():
    obj0 = []
    for i in range(0,50):
        obj0.append(utils.NominalTheBest(Value=i, LowerLimit=10, Target=20, UpperLimit=30))

    utils.full_pickle(os.path.join(tmp_path.name,'test_desirability'), obj0)
    utils.compressed_pickle(os.path.join(tmp_path.name,'test_desirability'), obj0)
    obj1 = utils.loosen(os.path.join(tmp_path.name,'test_desirability.pkl'))
    obj2 = utils.decompress_pickle(os.path.join(tmp_path.name,'test_desirability.pbz2'))

    assert obj0 == obj1
    assert obj0 == obj2


if __name__ == '__main__':
    test_local_command_line()