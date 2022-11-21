#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from rdkit import Chem
from moldrug import utils, fitness, home
from moldrug.data import receptor_pdbqt, ligands, boxes, receptor_pdb, constraintref
import tempfile, os, gzip, shutil, requests, yaml, copy, sys
import pytest
from multiprocessing import cpu_count


# Creating a temporal directory
tmp_path = tempfile.TemporaryDirectory()
# Creating receptors files
r_x0161_pdbqt_file = os.path.join(tmp_path.name, 'r_x0161.pdbqt')
r_x0161_pdb_file = os.path.join(tmp_path.name, 'r_x0161.pdb')
r_6lu7_pdbqt_file = os.path.join(tmp_path.name, 'r_6lu7.pdbqt')
with open(r_x0161_pdbqt_file, 'w') as r:
    r.write(receptor_pdbqt.r_x0161)
with open(r_x0161_pdb_file, 'w') as r:
    r.write(receptor_pdb.r_x0161)
with open(r_6lu7_pdbqt_file, 'w') as r:
    r.write(receptor_pdbqt.r_6lu7)

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
            "seed_mol": ligands.r_x0161,
            "AddHs": True,
            "costfunc": "Cost",
            "costfunc_kwargs": {
                "vina_executable": "vina",
                "receptor_pdbqt_path": r_x0161_pdbqt_file,
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
    p = utils.run('moldrug test_single_receptor.yml')
    print(p.stdout)
    os.chdir(cwd)


def test_multi_receptor(maxiter = 1, popsize = 2, njobs = 3, NumbCalls = 1):
    out = utils.GA(
        seed_mol=Chem.MolFromSmiles(ligands.r_x0161),
        AddHs= True,
        maxiter=maxiter,
        popsize=popsize,
        crem_db_path = crem_db_path,
        pc = 1,
        get_similar = True,
        mutate_crem_kwargs = {
            'radius':3,
            'min_size':0,
            'min_inc':-5,
            'max_inc':3,
            'replace_ids': [3, 4, 5, 7],
            'protected_ids': [0],
        },
        costfunc = fitness.CostMultiReceptors,
        costfunc_kwargs = {
            'receptor_pdbqt_path': [r_x0161_pdbqt_file, r_6lu7_pdbqt_file],
            'boxcenter' : [boxes.r_x0161["A"]['boxcenter'], boxes.r_6lu7["A"]['boxcenter']],
            'boxsize': [boxes.r_x0161["A"]['boxsize'], boxes.r_6lu7["A"]['boxsize']],
            'vina_score_type': ['min', 'max'],
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
        p.write(out.pop[0].pdbqt[0])

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
            "pick": 1,
            "AddHs": True,
            "seed_mol": Chem.MolToSmiles(Chem.AddHs(Chem.MolFromSmiles(ligands.r_x0161))),
            "costfunc": "CostOnlyVina",
            "costfunc_kwargs": {
                "vina_executable": "vina",
                "receptor_pdbqt_path": r_x0161_pdbqt_file,
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


@pytest.mark.filterwarnings("ignore:\nVina failed")
def test_fitness_module():
    individual = utils.Individual(Chem.MolFromSmiles(ligands.r_x0161))
    individual_corrupted = copy.deepcopy(individual)
    individual_corrupted.pdbqt = 'This is a corrupted pdbqt'
    receptor_pdbqt_path = [r_x0161_pdbqt_file,r_6lu7_pdbqt_file]
    boxcenter = [boxes.r_x0161['A']['boxcenter'], boxes.r_6lu7['A']['boxcenter']]
    boxsize = [boxes.r_x0161['A']['boxsize'], boxes.r_6lu7['A']['boxsize']]
    vina_score_type = ['min', 'max']

    fitness.Cost(
        Individual = copy.deepcopy(individual),
        wd = tmp_path.name,
        receptor_pdbqt_path = r_x0161_pdbqt_file,
        boxcenter = boxes.r_x0161['A']['boxcenter'],
        boxsize = boxes.r_x0161['A']['boxsize'],
        exhaustiveness = 4,
        ncores = 4)
    fitness.Cost(
        Individual = copy.deepcopy(individual),
        wd = tmp_path.name,
        receptor_pdbqt_path = r_x0161_pdbqt_file,
        boxcenter = boxes.r_x0161['A']['boxcenter'],
        boxsize = boxes.r_x0161['A']['boxsize'],
        exhaustiveness = 4,
        ncores = 4,
        constraint=True,
        constraint_ref=Chem.MolFromMolBlock(constraintref.r_x0161),
        constraint_receptor_pdb_path = r_x0161_pdb_file,
        )

    fitness.Cost(Individual = copy.deepcopy(
        individual_corrupted),wd = tmp_path.name,receptor_pdbqt_path = r_x0161_pdbqt_file,
        boxcenter = boxes.r_x0161['A']['boxcenter'], boxsize = boxes.r_x0161['A']['boxsize'],exhaustiveness = 4,ncores = 4)
    fitness.CostMultiReceptors(
        Individual = copy.deepcopy(individual_corrupted),wd = tmp_path.name,receptor_pdbqt_path = receptor_pdbqt_path,
        vina_score_type = vina_score_type, boxcenter = boxcenter,boxsize = boxsize,exhaustiveness = 4,ncores = 4)
    fitness.CostMultiReceptorsOnlyVina(
        Individual = copy.deepcopy(individual),wd = tmp_path.name,receptor_pdbqt_path = receptor_pdbqt_path,
        vina_score_type = vina_score_type, boxcenter = boxcenter,boxsize = boxsize,exhaustiveness = 4,ncores = 4)
    fitness.CostMultiReceptorsOnlyVina(
        Individual = copy.deepcopy(individual),wd = tmp_path.name,receptor_pdbqt_path = receptor_pdbqt_path,
        vina_score_type = vina_score_type, boxcenter = boxcenter,boxsize = boxsize,exhaustiveness = 4,ncores = 4, wt_cutoff=2)
    fitness.CostMultiReceptorsOnlyVina(Individual = copy.deepcopy(
        individual_corrupted),wd = tmp_path.name,receptor_pdbqt_path = receptor_pdbqt_path,
        vina_score_type = vina_score_type, boxcenter = boxcenter,boxsize = boxsize,exhaustiveness = 4,ncores = 4)


    fitness.CostOnlyVina(
        Individual = copy.deepcopy(individual),wd = tmp_path.name,receptor_pdbqt_path = r_x0161_pdbqt_file,
        boxcenter = boxes.r_x0161['A']['boxcenter'], boxsize = boxes.r_x0161['A']['boxsize'],exhaustiveness = 4,ncores = 4)
    fitness.CostOnlyVina(
        Individual = copy.deepcopy(individual),wd = tmp_path.name,receptor_pdbqt_path = r_x0161_pdbqt_file,
        boxcenter = boxes.r_x0161['A']['boxcenter'], boxsize = boxes.r_x0161['A']['boxsize'],exhaustiveness = 4,ncores = 4, wt_cutoff=2)
    fitness.CostOnlyVina(
        Individual = copy.deepcopy(individual_corrupted),wd = tmp_path.name,receptor_pdbqt_path = r_x0161_pdbqt_file,
        boxcenter = boxes.r_x0161['A']['boxcenter'], boxsize = boxes.r_x0161['A']['boxsize'],exhaustiveness = 4,ncores = 4)

    fitness.__get_mol_cost(
        mol = Chem.MolFromMolBlock(constraintref.r_x0161),wd = tmp_path.name,receptor_pdbqt_path = r_x0161_pdbqt_file,
        boxcenter = boxes.r_x0161['A']['boxcenter'], boxsize = boxes.r_x0161['A']['boxsize'],)

    # Clean
    os.remove('0_error.pbz2')


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
    I5 = utils.Individual(Chem.MolFromSmiles('CC'), cost = 10, pdbqt=I1.pdbqt)
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

# def test_local():
#     os.chdir(tmp_path.name)
#     local = utils.Local(
#         seed_mol = Chem.AddHs(Chem.MolFromSmiles(ligands.r_x0161)),
#         crem_db_path = crem_db_path,
#         costfunc=fitness.CostOnlyVina,
#         costfunc_kwargs = {
#             "vina_executable": "vina",
#             "receptor_pdbqt_path": r_x0161_pdbqt_file,
#             "boxcenter": boxes.r_x0161["A"]['boxcenter'],
#             "boxsize": boxes.r_x0161["A"]['boxsize'],
#             "exhaustiveness": 4,
#             "num_modes": 1
#         },
#         grow_crem_kwargs = {
#             "radius": 3,
#             "min_atoms": 1,
#             "max_atoms": 3
#         },
#         AddHs=True,
#     )
#     local(njobs=1,pick=1)

#     utils.make_sdf(local.pop, sdf_name = os.path.join(tmp_path.name,"local_pop"))
#     print(os.listdir(tmp_path.name))

#     print(local.to_dataframe())
def test_constraintconf():
    from moldrug.constraintconf import constraintconf
    with Chem.SDWriter(os.path.join(tmp_path.name, 'fix.sdf')) as w:
        mol = Chem.MolFromMolBlock(constraintref.r_x0161)
        w.write(mol)
    with open(os.path.join(tmp_path.name, 'mol.smi'), 'w') as f:
        f.write(ligands.r_x0161)

    constraintconf(
        pdb=r_x0161_pdb_file,
        smi = os.path.join(tmp_path.name, 'mol.smi'),
        fix= os.path.join(tmp_path.name, 'fix.sdf'),
        out = os.path.join(tmp_path.name, 'conf.sdf')
    )


if __name__ == '__main__':
    test_multi_receptor()
