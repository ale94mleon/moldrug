#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from moldrug import utils
from moldrug.fitness import _vinadock
from rdkit import Chem
from typing import Dict, List
from rdkit import RDLogger
from molskill.scorer import MolSkillScorer
RDLogger.DisableLog('rdApp.*')


def Cost(
    Individual:utils.Individual,
    wd:str = '.vina_jobs',
    vina_executable:str = 'vina',
    receptor_pdbqt_path:str = None,
    boxcenter:List[float] = None,
    boxsize:List[float] = None,
    exhaustiveness:int = 8,
    ad4map:str = None,
    ncores:int = 1,
    num_modes:int = 1,
    constraint:bool = False,
    constraint_type:str = 'score_only', # score_only, local_only
    constraint_ref:Chem.rdchem.Mol = None,
    constraint_receptor_pdb_path:str = None,
    constraint_num_conf:int = 100,
    constraint_minimum_conf_rms:int = 0.01,
    desirability:Dict = None,
    ):
    """
    This is the main Cost function of the module. It use the concept of desirability functions. The response variables are:

    #. `Vina score. <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3041641/>`_
    #. `MolSkill. <https://github.com/microsoft/molskill>`_
    #. `Synthetic accessibility score.  <https://jcheminf.biomedcentral.com/articles/10.1186/1758-2946-1-8)>`_

    If ad4map is set, the last version of vina (`releases <https://github.com/ccsb-scripps/AutoDock-Vina/releases>`_)
    must be installed. To see how to use AutoDock4 force fields in the new version of vina, follow
    `this tutorial <https://autodock-vina.readthedocs.io/en/latest/docking_zinc.html>_`

    Parameters
    ----------
    Individual : utils.Individual
        A Individual with the pdbqt attribute
    wd : str, optional
        The working directory to execute the docking jobs, by default '.vina_jobs'
    vina_executable : str, optional
        This is the name of the vina executable, could be a path to the binary object,
        which must have  execution permits (chmod a+x <your binary file>),  by default 'vina'
    receptor_pdbqt_path : str, optional
        Where the receptor pdbqt file is located, by default None
    boxcenter : list[float], optional
        A list of three floats with the definition of the center of the box in angstrom for docking (x, y, z), by default None
    boxsize : list[float], optional
        A list of three floats with the definition of the box size in angstrom of the docking box (x, y, z), by default None
    exhaustiveness : int, optional
        Parameter of vina that controls the accuracy of the docking searching, by default 8
    ad4map : str, optional
        Affinity maps for the autodock4.2 (ad4) or vina scoring function, by default None
    ncores : int, optional
        Number of cpus to use in Vina, by default 1
    num_modes : int, optional
        How many modes should Vina export, by default 1
    constraint : bool, optional
        Controls if constraint docking will be perform, by default False
    constraint_type : str, optional
        This is the type of constraint docking.
        Could be local_only (vina will perform local optimization and score the resulted pose)
        or score_only (in this case the provided pose by the internal conformer generator will only be scored),
        by default 'score_only'
    constraint_ref : Chem.rdchem.Mol, optional
        The part of the molecule that we would like to constraint, by default None
    constraint_receptor_pdb_path : str, optional
        The same as constraint_receptor_pdbqt_path but in pdb format, by default None
    constraint_num_conf : int, optional
        Maximum number of conformer to be generated internally by MolDrug , by default 100
    constraint_minimum_conf_rms : int, optional
        RMS to filter duplicate conformers, by default 0.01
    desirability : dict, optional
        The definition of the desirability to use for each used variable = [qed, sa_score, vina_score].
        Each variable only will accept the keys [w, and the name of the desirability function of :meth:`moldrug.utils.DerringerSuichDesirability`],
        by default None which means that it will be used:
        desirability = {
        'qed': {'w': 1,'LargerTheBest': {'LowerLimit': 0.1,'Target': 0.75, 'r': 1}
        },
        'sa_score': {'w': 1,'SmallerTheBest': {'Target': 3,'UpperLimit': 7,'r': 1}
        },
        'vina_score': {'w': 1,'SmallerTheBest': {'Target': -12,'UpperLimit': -6,'r': 1}
        }
        }
    Returns
    -------
    utils.Individual
        A new instance of the original Individual with the the new attributes:
        pdbqt, qed, vina_score, sa_score and cost.
        cost attribute will be a number between 0 and 1, been 0 the optimal value.
    Example
    -------
    .. ipython:: python

        from moldrug import utils, fitness
        from rdkit import Chem
        import tempfile, os
        from moldrug.data import ligands, boxes, receptor_pdbqt
        tmp_path = tempfile.TemporaryDirectory()
        ligand_mol = Chem.MolFromSmiles(ligands.r_x0161)
        I = utils.Individual(ligand_mol)
        receptor_path = os.path.join(tmp_path.name,'receptor.pdbqt')
        with open(receptor_path, 'w') as r: r.write(receptor_pdbqt.r_x0161)
        box = boxes.r_x0161['A']
        # Using the default desirability
        NewI = fitness.Cost(
            Individual = I,wd = tmp_path.name,
            receptor_pdbqt_path = receptor_path,boxcenter = box['boxcenter'],
            boxsize = box['boxsize'],exhaustiveness = 4,ncores = 4)
        print(NewI.cost, NewI.vina_score, NewI.qed, NewI.sa_score)
    """
    if not desirability:
        desirability = {
            'molskill_score': {
                'w': 1,
                'SmallerTheBest': {
                    'Target': -15,
                    'UpperLimit': 0,
                    'r': 1
                }
            },
            'sa_score': {
                'w': 1,
                'SmallerTheBest': {
                    'Target': 3,
                    'UpperLimit': 7,
                    'r': 1
                }
            },
            'vina_score': {
                'w': 1,
                'SmallerTheBest': {
                    'Target': -12,
                    'UpperLimit': -6,
                    'r': 1
                }
            }
        }


    # Multicriteria optimization,Optimization of Several Response Variables
    # Estimation of SkillScorer
    molskill_scorer = MolSkillScorer()
    Individual.molskill_score = molskill_scorer.score([Individual.smiles])

    # Getting synthetic accessibility score
    sascorer = utils.import_sascorer()
    Individual.sa_score = sascorer.calculateScore(Individual.mol)

    # Getting vina_score and update pdbqt
    Individual.vina_score, Individual.pdbqt = _vinadock(
        Individual = Individual,
        wd = wd,
        vina_executable = vina_executable,
        receptor_pdbqt_path =  receptor_pdbqt_path,
        boxcenter = boxcenter,
        boxsize = boxsize,
        exhaustiveness = exhaustiveness,
        ad4map = ad4map,
        ncores = ncores,
        num_modes = num_modes,
        constraint = constraint,
        constraint_type = constraint_type,
        constraint_ref = constraint_ref,
        constraint_receptor_pdb_path = constraint_receptor_pdb_path,
        constraint_num_conf = constraint_num_conf,
        constraint_minimum_conf_rms = constraint_minimum_conf_rms,
    )
    # Adding the cost using all the information of qed, sas and vina_cost
    # Construct the desirability
    # Quantitative estimation of drug-likeness (ranges from 0 to 1). We could use just the value perse, but using LargerTheBest we are more permissible.

    base = 1
    exponent = 0
    for variable in desirability:
        for key in desirability[variable]:
            if key == 'w':
                w = desirability[variable][key]
            elif key in utils.DerringerSuichDesirability():
                d = utils.DerringerSuichDesirability()[key](getattr(Individual, variable), **desirability[variable][key])
            else:
                raise RuntimeError(f"Inside the desirability dictionary you provided for the variable = {variable} "\
                f"a non implemented key = {key}. Only are possible: 'w' (standing for weight) and any "\
                f"possible Derringer-Suich desirability function: {utils.DerringerSuichDesirability().keys()}")
        base *= d**w
        exponent += w

    # We are using a geometric mean. And because we are minimizing we have to return
    Individual.cost = 1 - base**(1/exponent)
    return Individual