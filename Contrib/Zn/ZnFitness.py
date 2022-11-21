#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from copy import deepcopy
from moldrug import utils, constraintconf
from rdkit import Chem
from rdkit.Chem import QED
import os
import numpy as np
from typing import Dict, List
import warnings
from meeko import MoleculePreparation


def __vinadock(
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
    constraint_type = 'score_only', # score_only, local_only
    constraint_ref:Chem.rdchem.Mol = None,
    constraint_receptor_pdb_path:str = None,
    constraint_num_conf:int = 100,
    constraint_minimum_conf_rms:int = 0.01):
    """
    This function is intend to be used to perform docking for all the cost functions implemented here.
    If ad4map is set, the last version of vina (`releases <https://github.com/ccsb-scripps/AutoDock-Vina/releases>`_)
    must be installed. To see how to use AutoDock4 force fields in the new version of vina, follow
    `this tutorial <https://autodock-vina.readthedocs.io/en/latest/docking_zinc.html>_`

    Parameters
    ----------
    Individual : utils.Individual
        An Individual with the pdbqt attribute (only needed for non constraint docking).
    wd : str, optional
        The working directory to execute the docking jobs, by default '.vina_jobs'
    vina_executable : str, optional
        This is the name of the vina executable, could be a path to the binary object, by default 'vina'
    receptor_pdbqt_path : str, optional
        Where the receptor pdbqt file is located, by default None
    boxcenter : List[float], optional
        A list of three floats with the definition of the center of the box in angstrom for docking (x, y, z), by default None
    boxsize : List[float], optional
        A list of three floats with the definition of the box size in angstrom of the docking box (x, y, z), by default None
    exhaustiveness : int, optional
         Parameter of vina that controls the accuracy of the docking searching, by default 8
    ad4map : str, optional
        The path where the AD4 maps are, by default None
    ncores : int, optional
        Number of cpus to use in Vina, by default 1
    num_modes : int, optional
        How many modes should Vina export, by default 1
    constraint : bool, optional
        Controls if constraint docking will be perform, by default False
    constraint_type : str, optional
        This is the type of constraint docking. Could be local_only
        (vina will perform local optimization and score the resulted pose)
        or score_only (in this case the provided pose by the internal conformer
        generator will only be scored), by default 'score_only'
    constraint_ref : Chem.rdchem.Mol, optional
        The part of the molecule that we would like to constraint, by default None
    constraint_receptor_pdb_path : str, optional
        The same as constraint_receptor_pdbqt_path but in pdb format, by default None
    constraint_num_conf : int, optional
        Maximum number of conformer to be generated internally by MolDrug , by default 100
    constraint_minimum_conf_rms : int, optional
        RMS to filter duplicate conformers, by default 0.01

    Returns
    -------
    tuple
        A tuple with two elements:
        (vina score, pdbqt string)

    Raises
    ------
    Exception
        Inappropriate constraint_type. must be local_only or score_only. Only will be checked if constraint is set to True.
    """

    # Creating the working directory if needed
    if not os.path.exists(wd):
        os.makedirs(wd)
    # Creating the command line for vina
    cmd_vina_str = f"{vina_executable}"\
        f" --cpu {ncores} --exhaustiveness {exhaustiveness} --num_modes {num_modes}"

    if ad4map:
        cmd_vina_str += f" --scoring ad4 --maps {os.path.abspath(ad4map)}"
    else:
        cmd_vina_str += f" --receptor {receptor_pdbqt_path}"\
            f" --center_x {boxcenter[0]} --center_y {boxcenter[1]} --center_z {boxcenter[2]}"\
            f" --size_x {boxsize[0]} --size_y {boxsize[1]} --size_z {boxsize[2]}"\

    if constraint:
        # Check for the correct type of docking
        if constraint_type in ['score_only', 'local_only']:
            cmd_vina_str += f" --{constraint_type}"
        else:
            raise Exception("constraint_type only admit two possible values: score_only, local_only.")

        # Generate constrained conformer
        try:
            out_mol = constraintconf.generate_conformers(
                mol = Chem.RemoveHs(Individual.mol),
                ref_mol = Chem.RemoveHs(constraint_ref),
                num_conf = constraint_num_conf,
                #ref_smi=Chem.MolToSmiles(constraint_ref),
                minimum_conf_rms=constraint_minimum_conf_rms)
            out_mol = Chem.AddHs(out_mol)
        except Exception:
            vina_score_pdbqt = (np.inf, "NonValidConformer")
            return vina_score_pdbqt

        # Remove conformers that clash with the protein
        clash_filter = constraintconf.ProteinLigandClashFilter(protein_pdbpath = constraint_receptor_pdb_path, distance=1.5)
        clashIds = [conf.GetId() for conf in out_mol.GetConformers() if clash_filter(conf)]
        _ = [out_mol.RemoveConformer(clashId) for clashId in clashIds]

        # Check first if some valid conformer exist
        if len(out_mol.GetConformers()):
            vina_score_pdbqt = (np.inf, None)
            for conf in out_mol.GetConformers():
                temp_mol = deepcopy(out_mol)
                temp_mol.RemoveAllConformers()
                temp_mol.AddConformer(out_mol.GetConformer(conf.GetId()), assignId=True)

                preparator = MoleculePreparation()
                preparator.prepare(temp_mol)
                preparator.write_pdbqt_file(os.path.join(wd, f'{Individual.idx}_conf_{conf.GetId()}.pdbqt'))

                # Make a copy to the vina command string and add the out (is needed) and ligand options
                cmd_vina_str_tmp = cmd_vina_str[:]
                cmd_vina_str_tmp += f" --ligand {os.path.join(wd, f'{Individual.idx}_conf_{conf.GetId()}.pdbqt')}"
                if constraint_type == 'local_only':
                    cmd_vina_str_tmp += f" --out {os.path.join(wd, f'{Individual.idx}_conf_{conf.GetId()}_out.pdbqt')}"

                try:
                    cmd_vina_result = utils.run(cmd_vina_str_tmp)
                except Exception as e:
                    if os.path.isfile(receptor_pdbqt_path):
                        with open(receptor_pdbqt_path, 'r') as f:
                            receptor_str = f.read()
                    else:
                        receptor_str = None

                    error = {
                        'Exception': e,
                        'Individual': Individual,
                        f'used_mol_conf_{conf.GetId()}': temp_mol,
                        f'used_ligand_pdbqt_conf_{conf.GetId()}': preparator.write_pdbqt_string(),
                        'receptor_str': receptor_str,
                        'boxcenter': boxcenter,
                        'boxsize': boxsize,
                    }
                    utils.compressed_pickle(f'{Individual.idx}_conf_{conf.GetId()}_error', error)
                    warnings.warn(f"\nVina failed! Check: {Individual.idx}_conf_{conf.GetId()}_error.pbz2 file.\n")
                    vina_score_pdbqt = (np.inf, preparator.write_pdbqt_string())
                    return vina_score_pdbqt

                vina_score = np.inf
                for line in cmd_vina_result.stdout.split('\n'):
                    # Check over different vina versions
                    if line.startswith('Affinity'):
                        vina_score = float(line.split()[1])
                        break
                    elif 'Estimated Free Energy of Binding' in line:
                        vina_score = float(line.split(':')[1].split()[0])
                        break
                if vina_score < vina_score_pdbqt[0]:
                    if constraint_type == 'local_only':
                        if os.path.isfile(os.path.join(wd, f'{Individual.idx}_conf_{conf.GetId()}_out.pdbqt')):
                            with open(os.path.join(wd, f'{Individual.idx}_conf_{conf.GetId()}_out.pdbqt'), 'r') as f:
                                pdbqt = f.read()
                        else:
                            pdbqt = "NonExistedFileToRead"
                        vina_score_pdbqt = (vina_score, pdbqt)
                    else:
                        vina_score_pdbqt = (vina_score, preparator.write_pdbqt_string())
        else:
            vina_score_pdbqt = (np.inf, "NonValidConformer")
    # "Normal" docking
    else:
        cmd_vina_str += f" --ligand {os.path.join(wd, f'{Individual.idx}.pdbqt')} --out {os.path.join(wd, f'{Individual.idx}_out.pdbqt')}"
        with open(os.path.join(wd, f'{Individual.idx}.pdbqt'), 'w') as lig_pdbqt:
            lig_pdbqt.write(Individual.pdbqt)
        try:
            utils.run(cmd_vina_str)
        except Exception as e:
            if os.path.isfile(receptor_pdbqt_path):
                with open(receptor_pdbqt_path, 'r') as f:
                    receptor_str = f.read()
            else:
                receptor_str = None

            error = {
                'Exception': e,
                'Individual': Individual,
                'receptor_str': receptor_str,
                'boxcenter': boxcenter,
                'boxsize': boxsize,
            }
            utils.compressed_pickle(f'{Individual.idx}_error', error)
            warnings.warn(f"\nVina failed! Check: {Individual.idx}_error.pbz2 file.\n")

            vina_score_pdbqt = (np.inf, 'VinaFailed')
            return vina_score_pdbqt

        # Getting the information
        best_energy = utils.VINA_OUT(os.path.join(wd, f'{Individual.idx}_out.pdbqt')).BestEnergy()
        vina_score_pdbqt = (best_energy.freeEnergy, ''.join(best_energy.chunk))
    return vina_score_pdbqt

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
    constraint_type = 'score_only', # score_only, local_only
    constraint_ref:Chem.rdchem.Mol = None,
    constraint_receptor_pdb_path:str = None,
    constraint_num_conf:int = 100,
    constraint_minimum_conf_rms:int = 0.01,
    desirability:Dict = None,
    ):
    """
    This is the main Cost function of the module. It use the concept of desirability functions. The response variables are:

    #. `Vina score. <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3041641/>`_
    #. `Quantitative Estimation of Drug-likeness (QED). <https://www.nature.com/articles/nchem.1243>`_
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
        This is the name of the vina executable, could be a path to the binary object,  by default 'vina'
    receptor_pdbqt_path : str, optional
        Where the receptor pdbqt file is located, by default None
    boxcenter : list[float], optional
        A list of three floats with the definition of the center of the box in angstrom for docking (x, y, z), by default None
    boxsize : list[float], optional
        A list of three floats with the definition of the box size in angstrom of the docking box (x, y, z), by default None
    exhaustiveness : int, optional
        Parameter of vina that controls the accuracy of the docking searching, by default 8
    ad4map : str, optional
        The path where the AD4 maps are, by default None
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
            'qed': {
                'w': 1,
                'LargerTheBest': {
                    'LowerLimit': 0.1,
                    'Target': 0.75,
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

    sascorer = utils.import_sascorer()
    # multicriteria optimization,Optimization of Several Response Variables
    # Getting estimate of drug-likness
    Individual.qed = QED.weights_mean(Individual.mol)

    # Getting synthetic accessibility score
    Individual.sa_score = sascorer.calculateScore(Individual.mol)

    # Getting vina_score and update pdbqt
    Individual.vina_score, Individual.pdbqt = __vinadock(
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