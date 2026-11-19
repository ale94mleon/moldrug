#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import copy
import os
import random
import sys
import tempfile
from multiprocessing import cpu_count

import build_crem_test_db
import get_vina
import yaml
from rdkit import Chem
from rdkit.Chem import Descriptors

from moldrug import GA, home
from moldrug.opt import Local
from moldrug.data import get_data
from moldrug.fitness import (Cost, CostMultiReceptors,
                             CostMultiReceptorsOnlyVina, CostOnlyVina,
                             __get_mol_cost)
from moldrug.runner import Runner, RunnerMode
from moldrug.utils import (VINA_OUT, Individual, NominalTheBest,
                           compressed_pickle, decompress_pickle, format_duration,
                           full_pickle,
                           get_sim, lipinski_filter, lipinski_profile, loosen,
                           run, tar_errors)

try:
    from dask_jobqueue import SLURMCluster
    dask_available = True
except ImportError:
    dask_available = False


# This is in case vina is not found. You can change the vina executable here
vina_executable = os.path.abspath('vina')

if not os.path.isfile(vina_executable):
    get_vina.download()


# Creating a temporal directory
tmp_path = tempfile.TemporaryDirectory(dir='.')
wd = tmp_path.name
os.chdir(wd)

# Creating receptors files
TEST_DATA = {
    'x0161': get_data('x0161'),
    '6lu7': get_data('6lu7')
}

# Building a small crem data base.
# See tests/build_crem_test_db.py for why it is built instead of downloaded.
crem_db_path = build_crem_test_db.build(
    os.path.abspath('crem.db'),
    ncpu=cpu_count())


# Execution mode can be changed to "smaug" or "elwe" to calculate cost functions on Smaug or Elwe clusters.
# Both configurations are used only internally by the author, because none of the two clusters is publicly available.
# That being said, an example cluster configuration can be adopted to your needs.
#
# Note on Elwe: Worker nodes of Elwe are unable to reach the controller when moldrug is
# running directly on the frontend. One way to circumvent it is to run tests from an
# interactive job, eg. by calling `elwe slot`, `pytest tests -k test_multi_receptor`.
# Sending the main job to the queue works as well
#   sbatch --partition uds-hub --nodes 1 --mem-per-cpu 1GB --time 00:15:00 \
#       --ntasks-per-node 10 --gres gpu:1 --constraint 'XEON_E5_2630v4' \
#       --wrap 'bash -c "source .venv/bin/activate ; pytest tests -k test_multi_receptor"'
execution_mode = "multiprocessing"


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
                "max_size": 5,
                "min_inc": -5,
                "max_inc": 3,
                "ncores": 12
            },
            "maxiter": 1,
            "deffnm": "02_allow_grow"
        }
    }

    if execution_mode == "multiprocessing":
        pass
    elif execution_mode in ["smaug", "elwe"]:
        assert dask_available

        # This is currently unused
        if execution_mode == "smaug":
            Config['01_grow']['cluster'] = {
                "type": "SLURMCluster",
                "kwargs": {
                    "queue": "short",
                    "cores": 12,
                    "processes": 1,
                    "memory": "8GB",
                    "walltime": "00:15:00",
                    "job_extra_directives": ["--gpus 1"]
                }
            }
        elif execution_mode == "elwe":
            Config['01_grow']['cluster'] = {
                "type": "SLURMCluster",
                "kwargs": {
                    "queue": "uds-hub",
                    "cores": 10,
                    "processes": 1,
                    "memory": "8GB",
                    "walltime": "00:15:00",
                    "job_extra_directives": ["--gpus 1"]
                }
            }
    else:
        raise ValueError(f"Unknown mode {execution_mode}")

    with open("test_single_receptor.yml", 'w') as c:
        yaml.dump(Config, c)
    p = run('moldrug run test_single_receptor.yml')
    print(p.stdout)
    # Run a second time but with a seed population
    Config['01_grow']['seed_mol'] = ['02_allow_grow_pop.pbz2', '02_allow_grow_pop.pbz2']
    with open("test_single_receptor_init_pop.yml", 'w') as c:
        yaml.dump(Config, c)
    p = run('moldrug run test_single_receptor_init_pop.yml')
    print(p.stdout)


def test_dask_local():
    runner_instance = Runner(
        mode=RunnerMode.DASK_LOCAL,
        process_count=6
    )

    result = runner_instance.run(sum, entries=[[1, 2], [2, 3]])
    assert result == [3, 5]
    runner_instance.client.close()


def test_multi_receptor(maxiter=1, popsize=2, njobs=3, NumbCalls=1):
    out = GA(
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
        costfunc=CostMultiReceptors,
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

    if execution_mode == "multiprocessing":
        for _ in range(NumbCalls):
            out(njobs=njobs)
    elif execution_mode in ["smaug", "elwe"]:
        assert dask_available

        if execution_mode == "smaug":
            cluster = SLURMCluster(
                queue='short',
                cores=12,
                processes=1,
                memory='8GB',
                walltime='00:15:00',
                job_extra_directives=['--gpus 1']
            )
        else:
            cluster = SLURMCluster(
                queue='uds-hub',
                cores=10,
                processes=1,
                memory='8GB',
                walltime='00:15:00',
                job_extra_directives=['--gpus 1']
            )

        # Requesting "njobs" workers, each in a separate Slurm slot.
        cluster.scale(njobs)

        # Hint: By increasing the number of processes one can split Slurm allocations
        # into multiple workers. It might be useful when requesting entire nodes with
        # `--exclusive` switch and then running as many processes as GPUs.

        runner = Runner(RunnerMode.DASK_JOB_QUEUE_CLUSTER, dask_cluster=cluster)

        for _ in range(NumbCalls):
            out(runner=runner)
    else:
        assert False, f"Unknown mode: {execution_mode}"

    for o in out.pop:
        print(o.smiles, o.cost)
    out.pickle(f"result_test_multi_receptor_NumGens_{out.NumGens}_PopSize_{popsize}", compress=False)
    print(out.to_dataframe())

    with open('vina_out.pdbqt', 'w') as p:
        p.write(out.pop[0].pdbqt[0])

    vina_out = VINA_OUT('vina_out.pdbqt')
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

    run(f"moldrug run local_config.yml --fitness {os.path.join(home.home(), 'fitness.py')}")
    print(os.listdir())
    # This problems with the modules are not so convenient
    sys.path.append('.')
    result = decompress_pickle('local_result.pbz2')
    result.pickle('local_non_compress', compress=False)
    print(result.to_dataframe())


# @pytest.mark.filterwarnings("ignore:\nVina failed")
def test_fitness_module():
    individual = Individual(Chem.MolFromSmiles(TEST_DATA['x0161']['smiles']))
    individual_corrupted = copy.deepcopy(individual)
    individual_corrupted.pdbqt = 'This is a corrupted pdbqt'
    receptor_pdbqt_path = [TEST_DATA['x0161']['protein']['pdbqt'], TEST_DATA['6lu7']['protein']['pdbqt']]
    boxcenter = [TEST_DATA['x0161']['box']['boxcenter'], TEST_DATA['6lu7']['box']['boxcenter']]
    boxsize = [TEST_DATA['x0161']['box']['boxsize'], TEST_DATA['6lu7']['box']['boxsize']]
    vina_score_type = ['min', 'max']

    Cost(
        Individual=copy.deepcopy(individual),
        wd=wd,
        vina_executable=vina_executable,
        receptor_pdbqt_path=TEST_DATA['x0161']['protein']['pdbqt'],
        boxcenter=TEST_DATA['x0161']['box']['boxcenter'],
        boxsize=TEST_DATA['x0161']['box']['boxsize'],
        exhaustiveness=4,
        ncores=4)
    Cost(
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

    Cost(Individual=copy.deepcopy(individual_corrupted), wd=wd,
         receptor_pdbqt_path=TEST_DATA['x0161']['protein']['pdbqt'],
         vina_executable=vina_executable,
         boxcenter=TEST_DATA['x0161']['box']['boxcenter'], boxsize=TEST_DATA['x0161']['box']['boxsize'],
         exhaustiveness=4, ncores=4, vina_seed=1234)
    CostMultiReceptors(
        Individual=copy.deepcopy(individual_corrupted), wd=wd, receptor_pdbqt_path=receptor_pdbqt_path,
        vina_executable=vina_executable, vina_score_type=vina_score_type, boxcenter=boxcenter,
        boxsize=boxsize, exhaustiveness=4, ncores=4, vina_seed=1234)
    CostMultiReceptorsOnlyVina(
        Individual=copy.deepcopy(individual), wd=wd, receptor_pdbqt_path=receptor_pdbqt_path,
        vina_executable=vina_executable, vina_score_type=vina_score_type, boxcenter=boxcenter,
        boxsize=boxsize, exhaustiveness=4, ncores=4, vina_seed=1234)
    CostMultiReceptorsOnlyVina(
        Individual=copy.deepcopy(individual), wd=wd, receptor_pdbqt_path=receptor_pdbqt_path,
        vina_executable=vina_executable, vina_score_type=vina_score_type, boxcenter=boxcenter,
        boxsize=boxsize, exhaustiveness=4, ncores=4, wt_cutoff=2, vina_seed=1234)
    CostMultiReceptorsOnlyVina(Individual=copy.deepcopy(
        individual_corrupted), wd=wd, receptor_pdbqt_path=receptor_pdbqt_path,
        vina_executable=vina_executable, vina_score_type=vina_score_type, boxcenter=boxcenter,
        boxsize=boxsize, exhaustiveness=4, ncores=4, vina_seed=1234)

    CostOnlyVina(
        Individual=copy.deepcopy(individual), wd=wd,
        receptor_pdbqt_path=TEST_DATA['x0161']['protein']['pdbqt'], vina_executable=vina_executable,
        boxcenter=TEST_DATA['x0161']['box']['boxcenter'], boxsize=TEST_DATA['x0161']['box']['boxsize'],
        exhaustiveness=4, ncores=4, vina_seed=1234)
    CostOnlyVina(
        Individual=copy.deepcopy(individual), wd=wd,
        receptor_pdbqt_path=TEST_DATA['x0161']['protein']['pdbqt'], vina_executable=vina_executable,
        boxcenter=TEST_DATA['x0161']['box']['boxcenter'], boxsize=TEST_DATA['x0161']['box']['boxsize'],
        exhaustiveness=4, ncores=4, wt_cutoff=2, vina_seed=1234)
    CostOnlyVina(
        Individual=copy.deepcopy(individual_corrupted), wd=wd,
        receptor_pdbqt_path=TEST_DATA['x0161']['protein']['pdbqt'], vina_executable=vina_executable,
        boxcenter=TEST_DATA['x0161']['box']['boxcenter'], boxsize=TEST_DATA['x0161']['box']['boxsize'],
        exhaustiveness=4, ncores=4, vina_seed=1234)

    __get_mol_cost(
        mol=Chem.MolFromMolFile(TEST_DATA['x0161']['ligand_3D']), wd=wd,
        receptor_pdbqt_path=TEST_DATA['x0161']['protein']['pdbqt'], vina_executable=vina_executable,
        boxcenter=TEST_DATA['x0161']['box']['boxcenter'], boxsize=TEST_DATA['x0161']['box']['boxsize'],)

    # Clean
    tar_errors()
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
    get_sim(mols, ref_fps)


def test_lipinski():
    mol = Chem.MolFromSmiles('CCCO')
    lipinski_filter(Chem.MolFromSmiles('BrCC(COCN)CC(Br)CC(Cl)CCc1ccccc1CCCC(NCCCO)'))
    lipinski_filter(mol)
    lipinski_profile(mol)


def test_Individual():
    I1 = Individual(Chem.MolFromSmiles('CC'), cost=10)
    I2 = Individual(Chem.MolFromSmiles('CCO'), cost=2)
    I3 = copy.copy(I1)
    I4 = copy.deepcopy(I2)
    I5 = Individual(Chem.MolFromSmiles('CC'), cost=10, pdbqt=I1.pdbqt)
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
        obj0.append(NominalTheBest(Value=i, LowerLimit=10, Target=20, UpperLimit=30))

    full_pickle('test_desirability', obj0)
    compressed_pickle('test_desirability', obj0)
    obj1 = loosen('test_desirability.pkl')
    obj2 = decompress_pickle('test_desirability.pbz2')

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
    tar_errors()


def test_generate_conformers():
    from rdkit.Chem import AllChem

    from moldrug.constraintconf import generate_conformers
    ref = Chem.MolFromSmiles('O=S(=O)(Nc1ccc(Cl)cc1)c1ccsc1C(O)O')
    mol = Chem.MolFromSmiles('CN(C)S(=O)(=O)c1cc(NS(=O)(=O)c2ccsc2C(O)O)ccc1Cl')

    AllChem.EmbedMolecule(ref)
    AllChem.MMFFOptimizeMolecule(ref)
    generate_conformers(Chem.RemoveHs(mol), Chem.RemoveHs(ref), 50, randomseed=1234)

    # Clean
    tar_errors()


def molecular_weight_cost(Individual, **kwargs):
    """A cost function without docking, so that the reproducibility tests stay fast."""
    Individual.cost = Descriptors.MolWt(Individual.mol)
    return Individual


def reproducible_ga():
    """A short GA run, cheap enough to repeat, with every random decision seeded."""
    return GA(
        seed_mol=Chem.MolFromSmiles(TEST_DATA['x0161']['smiles']),
        costfunc=molecular_weight_cost,
        costfunc_kwargs={},
        crem_db_path=crem_db_path,
        maxiter=1,
        popsize=4,
        randomseed=42)


def test_randomseed_is_reproducible():
    """The same randomseed gives the same population, whatever the global random state is."""
    random.seed(1)
    first = reproducible_ga()
    first()

    random.seed(999)
    second = reproducible_ga()
    second()

    assert [individual.smiles for individual in first.pop] == \
        [individual.smiles for individual in second.pop]


def test_format_duration():
    """The duration is reported in the largest unit that keeps it above one."""
    assert format_duration(0.4) == '0.40 s'
    assert format_duration(59.9) == '59.90 s'
    # Each boundary belongs to the larger unit.
    assert format_duration(60) == '1.00 min'
    assert format_duration(3600) == '1.00 h'
    assert format_duration(86400) == '1.00 days'
    # Days are the largest unit, so a run of a week stays in days.
    assert format_duration(7 * 86400) == '7.00 days'


def test_timing_is_tracked():
    """Every generation records how long it spent generating and evaluating."""
    ga = reproducible_ga()
    ga()

    # Generation 0 is the initialization of the population, 1 is the only iteration.
    assert sorted(ga.timing) == [0, 1]
    for generation in ga.timing.values():
        assert generation['evaluation'] > 0
    # Nothing is generated during the initialization of a population built from the seed,
    # but the single iteration does mutate.
    assert ga.timing[1]['generation'] > 0


def test_local_timing_is_tracked():
    """Local times the same two phases, and accumulates them over calls."""
    local = Local(
        seed_mol=Chem.MolFromSmiles(TEST_DATA['x0161']['smiles']),
        costfunc=molecular_weight_cost,
        costfunc_kwargs={},
        crem_db_path=crem_db_path,
        grow_crem_kwargs={'radius': 3, 'min_atoms': 1, 'max_atoms': 3},
        randomseed=5)
    local(pick=2)

    assert sorted(local.timing) == ['evaluation', 'generation']
    assert all(seconds > 0 for seconds in local.timing.values())

    # A second call adds to the same counters rather than replacing them.
    first_call = dict(local.timing)
    local(pick=2)
    assert all(local.timing[phase] > first_call[phase] for phase in first_call)


def test_checkpoint_resume_is_reproducible():
    """A run restored from a checkpoint carries on from where its generator left off.

    The generator lives on the class, so it is pickled with it. Were it the global one,
    the two continuations below would diverge, because they are seeded differently.
    """
    ga = reproducible_ga()
    ga()
    ga.pickle('test_checkpoint_resume', compress=True)

    random.seed(1)
    first = decompress_pickle('test_checkpoint_resume.pbz2')
    first()

    random.seed(999)
    second = decompress_pickle('test_checkpoint_resume.pbz2')
    second()

    assert [individual.smiles for individual in first.pop] == \
        [individual.smiles for individual in second.pop]
    # And the resumed run is not simply repeating the generation it loaded.
    assert first.NumGens == ga.NumGens + 1


if __name__ == '__main__':
    pass
