#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import copy
import gzip
import os
import shutil
import sys
import tempfile
from multiprocessing import cpu_count

import get_vina
import requests
import yaml
from rdkit import Chem

from moldrug import fitness, home, utils
from moldrug.data import get_data

# This is in case vina is not found. You can change the vina executable here
vina_executable = os.path.abspath('vina')

if not os.path.isfile(vina_executable):
    get_vina.download()


# Creating a temporal directory
tmp_path = tempfile.TemporaryDirectory()
wd = tmp_path.name
os.chdir(wd)

# Creating receptors files
TEST_DATA = {
    'x0161': get_data('x0161'),
    '6lu7': get_data('6lu7')
}

# Getting the crem data base
url = "http://www.qsar4u.com/files/cremdb/replacements02_sc2.db.gz"
r = requests.get(url, allow_redirects=True)
crem_dbgz_path = 'crem.db.gz'
crem_db_path = 'crem.db'
open(crem_dbgz_path, 'wb').write(r.content)
with gzip.open(crem_dbgz_path, 'rb') as f_in:
    with open(crem_db_path, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)


def test_single_receptor_command_line():
    Config = {
        "01_grow": {
            "type": "GA",
            "njobs": 3,
            "seed_mol": TEST_DATA['x0161']['smiles'],
            "AddHs": True,
            "costfunc": "Cost",
            "costfunc_kwargs": {
                "vina_executable": vina_executable,
                "receptor_pdbqt_path": TEST_DATA['x0161']['protein']['pdbqt'],
                "boxcenter": TEST_DATA['x0161']['box']['boxcenter'],
                "boxsize": TEST_DATA['x0161']['box']['boxsize'],
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
    with open("test_single_receptor.yml", 'w') as c:
        yaml.dump(Config, c)
    p = utils.run('moldrug test_single_receptor.yml')
    print(p.stdout)
    # Run a second time but with a seed population
    Config['01_grow']['seed_mol'] = ['02_allow_grow_pop.pbz2', '02_allow_grow_pop.pbz2']
    with open("test_single_receptor_init_pop.yml", 'w') as c:
        yaml.dump(Config, c)
    p = utils.run('moldrug test_single_receptor_init_pop.yml')
    print(p.stdout)


def test_multi_receptor(maxiter=1, popsize=2, njobs=3, NumbCalls=1):
    out = utils.GA(
        seed_mol=[Chem.MolFromSmiles(TEST_DATA['x0161']['smiles']), Chem.MolFromSmiles(TEST_DATA['x0161']['smiles'])],
        AddHs=True,
        maxiter=maxiter,
        popsize=popsize,
        crem_db_path=crem_db_path,
        pc=1,
        get_similar=True,
        mutate_crem_kwargs={
            'radius': 3,
            'min_size': 0,
            'min_inc': -5,
            'max_inc': 3,
            'replace_ids': [3, 4, 5, 7],
            'protected_ids': [0],
        },
        costfunc=fitness.CostMultiReceptors,
        costfunc_kwargs={
            'receptor_pdbqt_path': [TEST_DATA['x0161']['protein']['pdbqt'], TEST_DATA['6lu7']['protein']['pdbqt']],
            'boxcenter': [TEST_DATA['x0161']['box']['boxcenter'], TEST_DATA['6lu7']['box']['boxcenter']],
            'boxsize': [TEST_DATA['x0161']['box']['boxsize'], TEST_DATA['6lu7']['box']['boxsize']],
            'vina_score_type': ['min', 'max'],
            'exhaustiveness': 4,
            'vina_executable': vina_executable,
            'ncores': int(cpu_count() / njobs),
            'num_modes': 1,
            'vina_seed': 1234,
        },
        save_pop_every_gen=20,
        deffnm='test_multi_receptor',
        randomseed=123)

    for _ in range(NumbCalls):
        out(njobs=njobs)

    for o in out.pop:
        print(o.smiles, o.cost)
    out.pickle(f"result_test_multi_receptor_NumGens_{out.NumGens}_PopSize_{popsize}", compress=False)
    print(out.to_dataframe())

    with open('vina_out.pdbqt', 'w') as p:
        p.write(out.pop[0].pdbqt[0])

    vina_out = utils.VINA_OUT('vina_out.pdbqt')
    vina_out.chunks[0].get_atoms()
    vina_out.chunks[0].write('chunk.pdbqt')
    vina_out.chunks[0].write()


def test_local_command_line():
    Config = {
        "main": {
            "type": "Local",
            "njobs": 1,
            "pick": 1,
            "AddHs": True,
            "seed_mol": Chem.MolToSmiles(Chem.AddHs(Chem.MolFromSmiles(TEST_DATA['x0161']['smiles']))),
            "costfunc": "CostOnlyVina",
            "costfunc_kwargs": {
                "vina_executable": vina_executable,
                "receptor_pdbqt_path": TEST_DATA['x0161']['protein']['pdbqt'],
                "boxcenter": TEST_DATA['x0161']['box']['boxcenter'],
                "boxsize": TEST_DATA['x0161']['box']['boxsize'],
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
    with open("local_config.yml", 'w') as c:
        yaml.dump(Config, c)

    utils.run(f"moldrug local_config.yml --fitness {os.path.join(home.home(), 'fitness.py')}")
    print(os.listdir())
    # This problems with the modules are not so convenient
    sys.path.append('.')
    result = utils.decompress_pickle('local_result.pbz2')
    result.pickle('local_non_compress', compress=False)
    print(result.to_dataframe())


# @pytest.mark.filterwarnings("ignore:\nVina failed")
def test_fitness_module():
    individual = utils.Individual(Chem.MolFromSmiles(TEST_DATA['x0161']['smiles']))
    individual_corrupted = copy.deepcopy(individual)
    individual_corrupted.pdbqt = 'This is a corrupted pdbqt'
    receptor_pdbqt_path = [TEST_DATA['x0161']['protein']['pdbqt'], TEST_DATA['6lu7']['protein']['pdbqt']]
    boxcenter = [TEST_DATA['x0161']['box']['boxcenter'], TEST_DATA['6lu7']['box']['boxcenter']]
    boxsize = [TEST_DATA['x0161']['box']['boxsize'], TEST_DATA['6lu7']['box']['boxsize']]
    vina_score_type = ['min', 'max']

    fitness.Cost(
        Individual=copy.deepcopy(individual),
        wd=wd,
        vina_executable=vina_executable,
        receptor_pdbqt_path=TEST_DATA['x0161']['protein']['pdbqt'],
        boxcenter=TEST_DATA['x0161']['box']['boxcenter'],
        boxsize=TEST_DATA['x0161']['box']['boxsize'],
        exhaustiveness=4,
        ncores=4)
    fitness.Cost(
        Individual=copy.deepcopy(individual),
        wd=wd,
        vina_executable=vina_executable,
        receptor_pdbqt_path=TEST_DATA['x0161']['protein']['pdbqt'],
        boxcenter=TEST_DATA['x0161']['box']['boxcenter'],
        boxsize=TEST_DATA['x0161']['box']['boxsize'],
        exhaustiveness=4,
        ncores=4,
        constraint=True,
        constraint_type='local_only',
        constraint_ref=Chem.MolFromMolFile(TEST_DATA['x0161']['ligand_3D']),
        constraint_receptor_pdb_path=TEST_DATA['x0161']['protein']['pdb'])

    fitness.Cost(Individual=copy.deepcopy(individual_corrupted), wd=wd,
                 receptor_pdbqt_path=TEST_DATA['x0161']['protein']['pdbqt'],
                 vina_executable=vina_executable,
                 boxcenter=TEST_DATA['x0161']['box']['boxcenter'], boxsize=TEST_DATA['x0161']['box']['boxsize'],
                 exhaustiveness=4, ncores=4, vina_seed=1234)
    fitness.CostMultiReceptors(
        Individual=copy.deepcopy(individual_corrupted), wd=wd, receptor_pdbqt_path=receptor_pdbqt_path,
        vina_executable=vina_executable, vina_score_type=vina_score_type, boxcenter=boxcenter,
        boxsize=boxsize, exhaustiveness=4, ncores=4, vina_seed=1234)
    fitness.CostMultiReceptorsOnlyVina(
        Individual=copy.deepcopy(individual), wd=wd, receptor_pdbqt_path=receptor_pdbqt_path,
        vina_executable=vina_executable, vina_score_type=vina_score_type, boxcenter=boxcenter,
        boxsize=boxsize, exhaustiveness=4, ncores=4, vina_seed=1234)
    fitness.CostMultiReceptorsOnlyVina(
        Individual=copy.deepcopy(individual), wd=wd, receptor_pdbqt_path=receptor_pdbqt_path,
        vina_executable=vina_executable, vina_score_type=vina_score_type, boxcenter=boxcenter,
        boxsize=boxsize, exhaustiveness=4, ncores=4, wt_cutoff=2, vina_seed=1234)
    fitness.CostMultiReceptorsOnlyVina(Individual=copy.deepcopy(
        individual_corrupted), wd=wd, receptor_pdbqt_path=receptor_pdbqt_path,
        vina_executable=vina_executable, vina_score_type=vina_score_type, boxcenter=boxcenter,
        boxsize=boxsize, exhaustiveness=4, ncores=4, vina_seed=1234)

    fitness.CostOnlyVina(
        Individual=copy.deepcopy(individual), wd=wd,
        receptor_pdbqt_path=TEST_DATA['x0161']['protein']['pdbqt'], vina_executable=vina_executable,
        boxcenter=TEST_DATA['x0161']['box']['boxcenter'], boxsize=TEST_DATA['x0161']['box']['boxsize'],
        exhaustiveness=4, ncores=4, vina_seed=1234)
    fitness.CostOnlyVina(
        Individual=copy.deepcopy(individual), wd=wd,
        receptor_pdbqt_path=TEST_DATA['x0161']['protein']['pdbqt'], vina_executable=vina_executable,
        boxcenter=TEST_DATA['x0161']['box']['boxcenter'], boxsize=TEST_DATA['x0161']['box']['boxsize'],
        exhaustiveness=4, ncores=4, wt_cutoff=2, vina_seed=1234)
    fitness.CostOnlyVina(
        Individual=copy.deepcopy(individual_corrupted), wd=wd,
        receptor_pdbqt_path=TEST_DATA['x0161']['protein']['pdbqt'], vina_executable=vina_executable,
        boxcenter=TEST_DATA['x0161']['box']['boxcenter'], boxsize=TEST_DATA['x0161']['box']['boxsize'],
        exhaustiveness=4, ncores=4, vina_seed=1234)

    fitness.__get_mol_cost(
        mol=Chem.MolFromMolFile(TEST_DATA['x0161']['ligand_3D']), wd=wd,
        receptor_pdbqt_path=TEST_DATA['x0161']['protein']['pdbqt'], vina_executable=vina_executable,
        boxcenter=TEST_DATA['x0161']['box']['boxcenter'], boxsize=TEST_DATA['x0161']['box']['boxsize'],)

    # Clean
    utils.tar_errors()
    os.remove('error.tar.gz')


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
    I1 = utils.Individual(Chem.MolFromSmiles('CC'), cost=10)
    I2 = utils.Individual(Chem.MolFromSmiles('CCO'), cost=2)
    I3 = copy.copy(I1)
    I4 = copy.deepcopy(I2)
    I5 = utils.Individual(Chem.MolFromSmiles('CC'), cost=10, pdbqt=I1.pdbqt)
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
    assert divmod(I1, I2) == (5, 0)
    assert I1**I2 == 100


def test_miscellanea():
    obj0 = []
    for i in range(0, 50):
        obj0.append(utils.NominalTheBest(Value=i, LowerLimit=10, Target=20, UpperLimit=30))

    utils.full_pickle('test_desirability', obj0)
    utils.compressed_pickle('test_desirability', obj0)
    obj1 = utils.loosen('test_desirability.pkl')
    obj2 = utils.decompress_pickle('test_desirability.pbz2')

    assert obj0 == obj1
    assert obj0 == obj2


def test_constraintconf():
    from moldrug.constraintconf import constraintconf
    with Chem.SDWriter('fix.sdf') as w:
        mol = Chem.MolFromMolFile(TEST_DATA['x0161']['ligand_3D'])
        w.write(mol)
    with open('mol.smi', 'w') as f:
        f.write(TEST_DATA['x0161']['smiles'])

    constraintconf(
        pdb=TEST_DATA['x0161']['protein']['pdb'],
        smi='mol.smi',
        fix='fix.sdf',
        out='conf.sdf',
        randomseed=1234,
    )
    # Clean
    utils.tar_errors()


def test_generate_conformers():
    from rdkit.Chem import AllChem

    from moldrug.constraintconf import generate_conformers
    ref = Chem.MolFromSmiles('O=S(=O)(Nc1ccc(Cl)cc1)c1ccsc1C(O)O')
    mol = Chem.MolFromSmiles('CN(C)S(=O)(=O)c1cc(NS(=O)(=O)c2ccsc2C(O)O)ccc1Cl')

    AllChem.EmbedMolecule(ref)
    AllChem.MMFFOptimizeMolecule(ref)
    generate_conformers(Chem.RemoveHs(mol), Chem.RemoveHs(ref), 50, randomseed=1234)

    # Clean
    utils.tar_errors()


if __name__ == '__main__':
    pass
