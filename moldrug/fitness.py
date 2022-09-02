#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from copy import deepcopy
from moldrug import utils, constraintconf
from rdkit import Chem
from rdkit.Chem import QED, Descriptors
import os
import numpy as np
from typing import Dict, List
import warnings
from meeko import MoleculePreparation

def get_mol_cost(
    mol,
    wd:str = '.vina_jobs',
    vina_executable:str = 'vina',
    receptor_pdbqt_path:str = None,
    boxcenter:List[float] = None,
    boxsize:List[float] = None,
    desirability:Dict = {
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
    ):
    if not os.path.exists(wd):
        os.makedirs(wd)
    # Initializing result dict
    results = dict()
    # Importing sascorer
    sascorer = utils.import_sascorer()
    # multicriteria optimization,Optimization of Several Response Variables

    # Getting estimate of drug-likness
    results['qed'] = QED.weights_mean(mol)

    # Getting synthetic accessibility score
    results['sa_score'] = sascorer.calculateScore(mol)

    # Getting vina_score and update pdbqt
    # Making the ligand pdbqt
    preparator = MoleculePreparation()
    preparator.prepare(mol)
    preparator.write_pdbqt_file(os.path.join(wd, 'ligand.pdbqt'))

    # Creating the command line for vina
    cmd_vina_str = f"{vina_executable} --receptor {receptor_pdbqt_path}"\
        f" --center_x {boxcenter[0]} --center_y {boxcenter[1]} --center_z {boxcenter[2]}"\
        f" --size_x {boxsize[0]} --size_y {boxsize[1]} --size_z {boxsize[2]}"\
        f" --score_only --ligand {os.path.join(wd, 'ligand.pdbqt')}"

    cmd_vina_result = utils.run(cmd_vina_str)
    for line in cmd_vina_result.stdout.split('\n'):
        # Check over different vina versions
        if line.startswith('Affinity'):
            results['vina_score'] = float(line.split()[1])
            break
        elif 'Estimated Free Energy of Binding' in line:
            results['vina_score'] = float(line.split(':')[1].split()[0])
            break    
    # Getting the desirability
    base = 1
    exponent = 0
    for variable in desirability:
        for key in desirability[variable]:
            if key == 'w':
                w = desirability[variable][key]
            elif key in utils.DerringerSuichDesirability():
                d = utils.DerringerSuichDesirability()[key](results[variable], **desirability[variable][key])
            else:
                raise RuntimeError(f"Inside the desirability dictionary you provided for the variable = {variable} a non implemented key = {key}. Only are possible: 'w' (standing for weight) and any possible Derringer-Suich desirability function: {utils.DerringerSuichDesirability().keys()}")
        base *= d**w
        exponent += w
    # Average
    #D = (w_qed*d_qed + w_sa_score*d_sa_score + w_vina_score*d_vina_score) / (w_qed + w_sa_score + w_vina_score)
    # Geometric mean
    # D = (d_qed**w_qed * d_sa_score**w_sa_score * d_vina_score**w_vina_score)**(1/(w_qed + w_sa_score + w_vina_score))
    # We are using Geometric mean
    D = base**(1/exponent)
    #  And because we are minimizing we have to return
    results['desirability'] = 1 - D
    return results

def vinadock(
    Individual:utils.Individual,
    wd:str = '.vina_jobs',
    vina_executable:str = 'vina',
    receptor_pdbqt_path:str = None,
    boxcenter:List[float] = None,
    boxsize:List[float] = None,
    exhaustiveness:int = 8,
    ncores:int = 1,
    num_modes:int = 1,
    constraint:bool = False,
    constraint_type = 'score_only', # score_only, local_only
    constraint_ref:Chem.rdchem.Mol = None,
    constraint_receptor_pdb_path:str = None,
    constraint_num_conf:int = 100,
    constraint_minimum_conf_rms:int = 0.01):
    """This function is intend to be used to perform docking fro all the cost functions implemented on :mod:`moldrug.fitness`

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
    ncores : int, optional
        Number of cpus to use in Vina, by default 1
    num_modes : int, optional
        How many modes should Vina export, by default 1
    constraint : bool, optional
        Controls if constraint docking will be perform, by default False
    constraint_type : str, optional
        This is the type of constraint docking. Could be local_only (vina will perform local optimization and score the resulted pose) or score_only (in this case the provided pose by the internal conformer generator will only be scored), by default 'score_only'
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
    cmd_vina_str = f"{vina_executable} --receptor {receptor_pdbqt_path}"\
        f" --center_x {boxcenter[0]} --center_y {boxcenter[1]} --center_z {boxcenter[2]}"\
        f" --size_x {boxsize[0]} --size_y {boxsize[1]} --size_z {boxsize[2]}"\
        f" --cpu {ncores} --exhaustiveness {exhaustiveness} --num_modes {num_modes}"

    if constraint:
        # Check for the correct type of docking
        if constraint_type in ['score_only', 'local_only']:
            cmd_vina_str += f" --{constraint_type}"
        else:
            raise Exception(f"constraint_type only admit two possible values: score_only, local_only.")

        # Generate constrained conformer
        try:
            out_mol = constraintconf.generate_conformers(
                mol = Chem.RemoveHs(Individual.mol),
                ref_mol = Chem.RemoveHs(constraint_ref),
                num_conf = constraint_num_conf,
                #ref_smi=Chem.MolToSmiles(constraint_ref),
                minimum_conf_rms=constraint_minimum_conf_rms)
        except Exception as e:
            vina_score_pdbqt = (np.inf, "NonValidConformer")
            return vina_score_pdbqt

        # Remove conformers that clash with the protein
        clash_filter = constraintconf.ProteinLigandClashFilter(protein_pdbpath = constraint_receptor_pdb_path, distance=1.5)
        clashIds = [conf.GetId() for conf in out_mol.GetConformers() if clash_filter(conf)]
        [out_mol.RemoveConformer(clashId) for clashId in clashIds]

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
        with open(os.path.join(wd, f'{Individual.idx}.pdbqt'), 'w') as l:
            l.write(Individual.pdbqt)
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

            vina_score_pdbqt = (np.inf, Individual.pdbqt)
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
    boxsize:List[float] =None,
    exhaustiveness:int = 8,
    ncores:int = 1,
    num_modes:int = 1,
    constraint:bool = False,
    constraint_type = 'score_only', # score_only, local_only
    constraint_ref:Chem.rdchem.Mol = None,
    constraint_receptor_pdb_path:str = None,
    constraint_num_conf:int = 100,
    constraint_minimum_conf_rms:int = 0.01,
    desirability:Dict = {
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
    ):
    """This is the main Cost function of the module. It use the concept of desirability functions. The response variables are:

    #. `Vina score. <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3041641/>`_
    #. `Quantitative Estimation of Drug-likeness (QED). <https://www.nature.com/articles/nchem.1243>`_
    #. `Synthetic accessibility score.  <https://jcheminf.biomedcentral.com/articles/10.1186/1758-2946-1-8)>`_

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
    ncores : int, optional
        Number of cpus to use in Vina, by default 1
    num_modes : int, optional
        How many modes should Vina export, by default 1
    constraint : bool, optional
        Controls if constraint docking will be perform, by default False
    constraint_type : str, optional
        This is the type of constraint docking. Could be local_only (vina will perform local optimization and score the resulted pose) or score_only (in this case the provided pose by the internal conformer generator will only be scored), by default 'score_only'
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
        Each variable only will accept the keys [w, and the name of the desirability function of :meth:`moldrug.utils.DerringerSuichDesirability`]
        ,by default { 'qed': { 'w': 1, 'LargerTheBest': { 'LowerLimit': 0.1, 'Target': 0.75, 'r': 1 } }, 'sa_score': { 'w': 1, 'SmallerTheBest': { 'Target': 3, 'UpperLimit': 7, 'r': 1 } }, 'vina_score': { 'w': 1, 'SmallerTheBest': { 'Target': -12, 'UpperLimit': -6, 'r': 1 } } }

    Returns
    -------
    utils.Individual
        A new instance of the original Individual with the the new attributes: pdbqt, qed, vina_score, sa_score and cost. cost attribute will be a number between 0 and 1, been 0 the optimal value.
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
        NewI = fitness.Cost(Individual = I,wd = tmp_path.name,receptor_pdbqt_path = receptor_path,boxcenter = box['boxcenter'],boxsize = box['boxsize'],exhaustiveness = 4,ncores = 4)
        print(NewI.cost, NewI.vina_score, NewI.qed, NewI.sa_score)
    """
    sascorer = utils.import_sascorer()
    # multicriteria optimization,Optimization of Several Response Variables
    # Getting estimate of drug-likness
    Individual.qed = QED.weights_mean(Individual.mol)

    # Getting synthetic accessibility score
    Individual.sa_score = sascorer.calculateScore(Individual.mol)

    # Getting vina_score and update pdbqt
    Individual.vina_score, Individual.pdbqt = vinadock(
        Individual = Individual,
        wd = wd,
        vina_executable = vina_executable,
        receptor_pdbqt_path =  receptor_pdbqt_path,
        boxcenter = boxcenter,
        boxsize = boxsize,
        exhaustiveness = exhaustiveness,
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
    # Quantitative estimation of drug-likness (ranges from 0 to 1). We could use just the value perse, but using LargerTheBest we are more permissible.

    base = 1
    exponent = 0
    for variable in desirability:
        for key in desirability[variable]:
            if key == 'w':
                w = desirability[variable][key]
            elif key in utils.DerringerSuichDesirability():
                d = utils.DerringerSuichDesirability()[key](getattr(Individual, variable), **desirability[variable][key])
            else:
                raise RuntimeError(f"Inside the desirability dictionary you provided for the variable = {variable} a non implemented key = {key}. Only are possible: 'w' (standing for weight) and any possible Derringer-Suich desirability function: {utils.DerringerSuichDesirability().keys()}")
        base *= d**w
        exponent += w
    # Average
    #D = (w_qed*d_qed + w_sa_score*d_sa_score + w_vina_score*d_vina_score) / (w_qed + w_sa_score + w_vina_score)
    # Geometric mean
    # D = (d_qed**w_qed * d_sa_score**w_sa_score * d_vina_score**w_vina_score)**(1/(w_qed + w_sa_score + w_vina_score))
    # We are using Geometric mean
    D = base**(1/exponent)
    #  And because we are minimizing we have to return
    Individual.cost = 1 - D
    return Individual

def CostOnlyVina(
    Individual:utils.Individual,
    wd:str = '.vina_jobs',
    vina_executable:str = 'vina',
    receptor_pdbqt_path:str = None,
    boxcenter:List[float] = None,
    boxsize:List[float] =None,
    exhaustiveness:int = 8,
    ncores:int = 1,
    num_modes:int = 1,
    constraint:bool = False,
    constraint_type = 'score_only', # score_only, local_only
    constraint_ref:Chem.rdchem.Mol = None,
    constraint_receptor_pdb_path:str = None,
    constraint_num_conf:int = 100,
    constraint_minimum_conf_rms:int = 0.01,
    wt_cutoff:float = None,
    ):
    """This Cost function performs Docking and return the vina_score as cost.

    Parameters
    ----------
    Individual : utils.Individual
        A Individual with the pdbqt attribute
    wd : str, optional
        The working directory to execute the docking jobs, by default '.vina_jobs'
    vina_executable : str, optional
        This is the name of the vina executable, could be a path to the binary object (x, y, z),  by default 'vina'
    receptor_path : str, optional
        Where the receptor pdbqt file is located, by default None
    boxcenter : list[float], optional
        A list of three floats with the definition of the center of the box in angstrom for docking (x, y, z), by default None
    boxsize : list[float], optional
        A list of three floats with the definition of the box size in angstrom of the docking box (x, y, z), by default None
    exhaustiveness : int, optional
        Parameter of vina that controls the accuracy of the docking searching, by default 8
    ncores : int, optional
        Number of cpus to use in Vina, by default 1
    num_modes : int, optional
        How many modes should Vina export, by default 1
    constraint : bool, optional
        Controls if constraint docking will be perform, by default False
    constraint_type : str, optional
        This is the type of constraint docking. Could be local_only (vina will perform local optimization and score the resulted pose) or score_only (in this case the provided pose by the internal conformer generator will only be scored), by default 'score_only'
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
        Each variable only will accept the keys [w, and the name of the desirability function of :meth:`moldrug.utils.DerringerSuichDesirability`]
        ,by default { 'qed': { 'w': 1, 'LargerTheBest': { 'LowerLimit': 0.1, 'Target': 0.75, 'r': 1 } }, 'sa_score': { 'w': 1, 'SmallerTheBest': { 'Target': 3, 'UpperLimit': 7, 'r': 1 } }, 'vina_score': { 'w': 1, 'SmallerTheBest': { 'Target': -12, 'UpperLimit': -6, 'r': 1 } } }
    wt_cutoff : float, optional
        If some number is provided the molecules with a molecular weight higher than wt_cutoff will get as vina_score = cost = np.inf. Vina will not be invoked, by default None
    Returns
    -------
    utils.Individual
        A new instance of the original Individual with the the new attributes: pdbqt, vina_score and cost. In this case cost = vina_score, the lowest the values the best individual.
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
        NewI = fitness.CostOnlyVina(Individual = I,wd = tmp_path.name,receptor_pdbqt_path = receptor_path,boxcenter = box['boxcenter'],boxsize = box['boxsize'],exhaustiveness = 4,ncores = 4)
        print(NewI.cost, NewI.vina_score)
    """
    # If the molecule is heavy, don't perform docking and assign infinite to the cost attribute. Add the pdbqt to pdbqts and np.inf to vina_scores
    if wt_cutoff:
        if Descriptors.MolWt(Individual.mol) > wt_cutoff:
            Individual.vina_score = np.inf
            Individual.cost = np.inf
            return Individual

    # Getting vina_score and update pdbqt
    Individual.vina_score, Individual.pdbqt = vinadock(
        Individual = Individual,
        wd = wd,
        vina_executable = vina_executable,
        receptor_pdbqt_path =  receptor_pdbqt_path,
        boxcenter = boxcenter,
        boxsize = boxsize,
        exhaustiveness = exhaustiveness,
        ncores = ncores,
        num_modes = num_modes,
        constraint = constraint,
        constraint_type = constraint_type,
        constraint_ref = constraint_ref,
        constraint_receptor_pdb_path = constraint_receptor_pdb_path,
        constraint_num_conf = constraint_num_conf,
        constraint_minimum_conf_rms = constraint_minimum_conf_rms,
    )
    Individual.cost = Individual.vina_score
    return Individual



def CostMultiReceptors(
    Individual:utils.Individual,
    wd:str = '.vina_jobs',
    vina_executable:str = 'vina',
    receptor_pdbqt_path:List[str] = None,
    vina_score_type:List[str] = None,
    boxcenter:List[float] = None,
    boxsize:List[float] =None,
    exhaustiveness:int = 8,
    ncores:int = 1,
    num_modes:int = 1,
    constraint:bool = False,
    constraint_type = 'score_only', # score_only, local_only
    constraint_ref:Chem.rdchem.Mol = None,
    constraint_receptor_pdb_path:List[str] = None,
    constraint_num_conf:int = 100,
    constraint_minimum_conf_rms:int = 0.01,
    desirability:Dict = {
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
        'vina_scores': {
            'min':{
                'w': 1,
                'SmallerTheBest': {
                    'Target': -12,
                    'UpperLimit': -6,
                    'r': 1
                }
            },
            'max':{
                'w': 1,
                'LargerTheBest': {
                    'LowerLimit': -4,
                    'Target': 0,
                    'r': 1
                }
            }
        }
    }
    ):
    """This function is similar to :meth:`moldrug.fitness.Cost` but it will add the possibility to work with more than one receptor. It also use the concept of desirability and the response variables are:

    #. `Vina scores. <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3041641/>`_
    #. `Quantitative Estimation of Drug-likeness (QED). <https://www.nature.com/articles/nchem.1243>`_
    #. `Synthetic accessibility score.  <https://jcheminf.biomedcentral.com/articles/10.1186/1758-2946-1-8)>`_

    In this case every vina score (for all the provided receptors) will be used for the construction of the desirability.

    Parameters
    ----------
    Individual : utils.Individual
        A Individual with the pdbqt attribute
    wd : str, optional
        The working directory to execute the docking jobs, by default '.vina_jobs'
    vina_executable : str, optional
        This is the name of the vina executable, could be a path to the binary object (x, y, z), by default 'vina'
    receptor_pdbqt_path : list[str], optional
        A list of location of the receptors pdbqt files, by default None
    vina_score_type : list[str], optional
        This is a list with the keywords 'min' and/or 'max'. E.g. If two receptor were provided and for the first one we would like to find a minimum in the vina scoring function and for the other one a maximum (selectivity for the first receptor); we must provided the list: ['min', 'max'], by default None
    boxcenter : list[float], optional
        A list of three floats with the definition of the center of the box in angstrom for docking (x, y, z), by default None
    boxsize : list[float], optional
        A list of three floats with the definition of the box size in angstrom of the docking box (x, y, z), by default None
    exhaustiveness : int, optional
        Parameter of vina that controls the accuracy of the docking searching, by default 8
    ncores : int, optional
         Number of cpus to use in Vina, by default 1
    num_modes : int, optional
        How many modes should Vina export, by default 1
    constraint : bool, optional
        Controls if constraint docking will be perform, by default False
    constraint_type : str, optional
        This is the type of constraint docking. Could be local_only (vina will perform local optimization and score the resulted pose) or score_only (in this case the provided pose by the internal conformer generator will only be scored), by default 'score_only'
    constraint_ref : Chem.rdchem.Mol, optional
        The part of the molecule that we would like to constraint, by default None
    constraint_receptor_pdb_path : list[str], optional
        The same as constraint_receptor_pdbqt_path but in pdb format, by default None
    constraint_num_conf : int, optional
        Maximum number of conformer to be generated internally by MolDrug , by default 100
    constraint_minimum_conf_rms : int, optional
        RMS to filter duplicate conformers, by default 0.01
    desirability : dict, optional
        The definition of the desirability to use for each used variable = [qed, sa_score, vina_scores].
        Each variable only will accept the keys [w, and the name of the desirability function of :meth:`moldrug.utils.DerringerSuichDesirability`].
        In the case of vina_scores there is another layer for the vina_score_type= [min, max],
        by default { 'qed': { 'w': 1, 'LargerTheBest': { 'LowerLimit': 0.1, 'Target': 0.75, 'r': 1 } }, 'sa_score': { 'w': 1, 'SmallerTheBest': { 'Target': 3, 'UpperLimit': 7, 'r': 1 } }, 'vina_score': { 'min':{ 'w': 1, 'SmallerTheBest': { 'Target': -12, 'UpperLimit': -6, 'r': 1 } }, 'max':{ 'w': 1, 'LargerTheBest': { 'LowerLimit': -4, 'Target': 0, 'r': 1 } } } }

    Returns
    -------
    utils.Individual
        A new instance of the original Individual with the the new attributes: pdbqts [a list of pdbqt], qed, vina_scores [a list of vina_score], sa_score and cost. cost attribute will be a number between 0 and 1, been 0 the optimal value.

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
        receptor_paths = [os.path.join(tmp_path.name,'receptor1.pdbqt'),os.path.join(tmp_path.name,'receptor2.pdbqt')]
        with open(receptor_paths[0], 'w') as r: r.write(receptor_pdbqt.r_x0161)
        with open(receptor_paths[1], 'w') as r: r.write(receptor_pdbqt.r_6lu7)
        boxcenter = [boxes.r_x0161['A']['boxcenter'], boxes.r_6lu7['A']['boxcenter']]
        boxsize = [boxes.r_x0161['A']['boxsize'], boxes.r_6lu7['A']['boxsize']]
        vina_score_type = ['min', 'max']
        # Using the default desirability
        NewI = fitness.CostMultiReceptors(Individual = I,wd = tmp_path.name,receptor_pdbqt_path = receptor_paths, vina_score_type = vina_score_type, boxcenter = boxcenter,boxsize = boxsize,exhaustiveness = 4,ncores = 4)
        print(NewI.cost, NewI.vina_score, NewI.qed, NewI.sa_score)
    """
    sascorer = utils.import_sascorer()
    Individual.qed = QED.weights_mean(Individual.mol)

    # Getting synthetic accessibility score
    Individual.sa_score = sascorer.calculateScore(Individual.mol)

    # Getting Vina score
    pdbqt_list = []
    Individual.vina_score = []
    for i in range(len(receptor_pdbqt_path)):
    # Getting vina_score and update pdbqt
        if constraint:
            vina_score, pdbqt = vinadock(
                Individual = Individual,
                wd = wd,
                vina_executable = vina_executable,
                receptor_pdbqt_path =  receptor_pdbqt_path[i],
                boxcenter = boxcenter[i],
                boxsize = boxsize[i],
                exhaustiveness = exhaustiveness,
                ncores = ncores,
                num_modes = num_modes,
                constraint = constraint,
                constraint_type = constraint_type,
                constraint_ref = constraint_ref,
                constraint_receptor_pdb_path = constraint_receptor_pdb_path[i],
                constraint_num_conf = constraint_num_conf,
                constraint_minimum_conf_rms = constraint_minimum_conf_rms,
            )
        else:
            vina_score, pdbqt = vinadock(
                Individual = Individual,
                wd = wd,
                vina_executable = vina_executable,
                receptor_pdbqt_path =  receptor_pdbqt_path[i],
                boxcenter = boxcenter[i],
                boxsize = boxsize[i],
                exhaustiveness = exhaustiveness,
                ncores = ncores,
                num_modes = num_modes,
            )
        Individual.vina_score.append(vina_score)
        pdbqt_list.append(pdbqt)
    # Update the pdbqt attribute
    Individual.pdbqt = pdbqt_list

    # make a copy of the default values of desirability
    # pops the region of vina_scores
    desirability_to_work_with = desirability.copy()
    vina_desirability_section = desirability_to_work_with.pop('vina_scores')
    # Initialize base and exponent
    base = 1
    exponent = 0
    # Runs for all properties different to vina_scores
    for variable in desirability_to_work_with:
        for key in desirability_to_work_with[variable]:
            if key == 'w':
                w = desirability_to_work_with[variable][key]
            elif key in utils.DerringerSuichDesirability():
                d = utils.DerringerSuichDesirability()[key](getattr(Individual, variable), **desirability_to_work_with[variable][key])
            else:
                raise RuntimeError(f"Inside the desirability dictionary you provided for the variable = {variable} a non implemented key = {key}. Only are possible: 'w' (standing for weight) and any possible Derringer-Suich desirability function: {utils.DerringerSuichDesirability().keys()}. Only in the case of vina_scores [min and max] keys")
        base *= d**w
        exponent += w

    # Run only for vina_scores
    for vs, vst in zip(Individual.vina_score, vina_score_type):
        for key in vina_desirability_section[vst]:
            if key == 'w':
                w = vina_desirability_section[vst][key]
            elif key in utils.DerringerSuichDesirability():
                d = utils.DerringerSuichDesirability()[key](vs, **vina_desirability_section[vst][key])
            else:
                raise RuntimeError(f"Inside the desirability dictionary you provided for the variable = vina_scores[{vst}] a non implemented key = {key}. Only are possible: 'w' (standing for weight) and any possible Derringer-Suich desirability function: {utils.DerringerSuichDesirability().keys()}.")
        base *= d**w
        exponent += w

    # Average
    #D = (w_qed*d_qed + w_sa_score*d_sa_score + sum([w_vina_score*d_vina_score for w_vina_score, d_vina_score in zip(w_vina_scores, d_vina_scores)])) / (w_qed + w_sa_score + sum(w_vina_scores))
    # Geometric mean
    # D = (d_qed**w_qed * d_sa_score**w_sa_score * np.prod([d_vina_score**w_vina_score for d_vina_score, w_vina_score in zip(d_vina_scores, w_vina_scores)]))** (1/(w_qed + w_sa_score + sum(w_vina_scores)))
    D = base**(1/exponent)
    # And because we are minimizing we have to return
    Individual.cost = 1 - D
    return Individual

def CostMultiReceptorsOnlyVina(
    Individual:utils.Individual,
    wd:str = '.vina_jobs',
    vina_executable:str = 'vina',
    receptor_pdbqt_path:List[str] = None,
    vina_score_type:List[str] = None,
    boxcenter:List[float] = None,
    boxsize:List[float] =None,
    exhaustiveness:int = 8,
    ncores:int = 1,
    num_modes:int = 1,
    constraint:bool = False,
    constraint_type = 'score_only', # score_only, local_only
    constraint_ref:Chem.rdchem.Mol = None,
    constraint_receptor_pdb_path:List[str] = None,
    constraint_num_conf:int = 100,
    constraint_minimum_conf_rms:int = 0.01,
    desirability:Dict = {
        'min':{
            'w': 1,
            'SmallerTheBest': {
                'Target': -12,
                'UpperLimit': -6,
                'r': 1
            }
        },
        'max':{
            'w': 1,
            'LargerTheBest': {
                'LowerLimit': -4,
                'Target': 0,
                'r': 1
            }
        }
    },
    wt_cutoff = None,
    ):
    """This function is similar to :meth:`moldrug.fitness.CostOnlyVina` but it will add the possibility to work with more than one receptor. It also use the concept of desirability.
    The response variables are the vina scores on each receptor.

    Parameters
    ----------
    Individual : utils.Individual
        A Individual with the pdbqt attribute
    wd : str, optional
        The working directory to execute the docking jobs, by default '.vina_jobs'
    vina_executable : str, optional
        This is the name of the vina executable, could be a path to the binary object (x, y, z), by default 'vina'
    receptor_pdbqt_path : list[str], optional
        A list of location of the receptors pdbqt files, by default None
    vina_score_type : list[str], optional
        This is a list with the keywords 'min' and/or 'max'. E.g. If two receptor were provided and for the first one we would like to find a minimum in the vina scoring function and for the other one a maximum (selectivity for the first receptor); we must provided the list: ['min', 'max'], by default None
    boxcenter : list[float], optional
        A list of three floats with the definition of the center of the box in angstrom for docking (x, y, z), by default None
    boxsize : list[float], optional
        A list of three floats with the definition of the box size in angstrom of the docking box (x, y, z), by default None
    exhaustiveness : int, optional
        Parameter of vina that controls the accuracy of the docking searching, by default 8
    ncores : int, optional
         Number of cpus to use in Vina, by default 1
    num_modes : int, optional
        How many modes should Vina export, by default 1
    constraint : bool, optional
        Controls if constraint docking will be perform, by default False
    constraint_type : str, optional
        This is the type of constraint docking. Could be local_only (vina will perform local optimization and score the resulted pose) or score_only (in this case the provided pose by the internal conformer generator will only be scored), by default 'score_only'
    constraint_ref : Chem.rdchem.Mol, optional
        The part of the molecule that we would like to constraint, by default None
    constraint_receptor_pdb_path : list[str], optional
        The same as constraint_receptor_pdbqt_path but in pdb format, by default None
    constraint_num_conf : int, optional
        Maximum number of conformer to be generated internally by MolDrug , by default 100
    constraint_minimum_conf_rms : int, optional
        RMS to filter duplicate conformers, by default 0.01
    desirability : dict, optional
        The definition of the desirability when min or max is used.
        Each variable only will accept the keys [w, and the name of the desirability function of :meth:`moldrug.utils.DerringerSuichDesirability`].
        by default { 'min':{ 'w': 1, 'SmallerTheBest': { 'Target': -12, 'UpperLimit': -6, 'r': 1 } }, 'max':{ 'w': 1, 'LargerTheBest': { 'LowerLimit': -4, 'Target': 0, 'r': 1 } } }
    wt_cutoff : float, optional
        If some number is provided the molecules with a molecular weight higher than wt_cutoff will get as vina_score = cost = np.inf. Vina will not be invoked, by default None

    Returns
    -------
    utils.Individual
        A new instance of the original Individual with the the new attributes: pdbqts [a list of pdbqt], vina_scores [a list of vina_score], and cost. cost attribute will be a number between 0 and 1, been 0 the optimal value.

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
        receptor_paths = [os.path.join(tmp_path.name,'receptor1.pdbqt'),os.path.join(tmp_path.name,'receptor2.pdbqt')]
        with open(receptor_paths[0], 'w') as r: r.write(receptor_pdbqt.r_x0161)
        with open(receptor_paths[1], 'w') as r: r.write(receptor_pdbqt.r_6lu7)
        boxcenter = [boxes.r_x0161['A']['boxcenter'], boxes.r_6lu7['A']['boxcenter']]
        boxsize = [boxes.r_x0161['A']['boxsize'], boxes.r_6lu7['A']['boxsize']]
        vina_score_type = ['min', 'max']
        # Using the default desirability
        NewI = fitness.CostMultiReceptorsOnlyVina(Individual = I,wd = tmp_path.name,receptor_pdbqt_path = receptor_paths, vina_score_type = vina_score_type, boxcenter = boxcenter,boxsize = boxsize,exhaustiveness = 4,ncores = 4)
        print(NewI.cost, NewI.vina_score)
    """
    pdbqt_list = []
    Individual.vina_score = []

    # If the molecule is heavy, don't perform docking and assign infinite to the cost attribute. Add the pdbqt to pdbqts and np.inf to vina_scores
    if wt_cutoff:
        if Descriptors.MolWt(Individual.mol) > wt_cutoff:
            for _ in range(len(receptor_pdbqt_path)):
                pdbqt_list.append(Individual.pdbqt)
                Individual.vina_score.append(np.inf)
            Individual.cost = np.inf
            # Update the pdbqt attribute
            Individual.pdbqt = pdbqt_list
            return Individual

    # Getting Vina score
    pdbqt_list = []
    Individual.vina_score = []
    for i in range(len(receptor_pdbqt_path)):
    # Getting vina_score and update pdbqt
        if constraint:
            vina_score, pdbqt = vinadock(
                Individual = Individual,
                wd = wd,
                vina_executable = vina_executable,
                receptor_pdbqt_path =  receptor_pdbqt_path[i],
                boxcenter = boxcenter[i],
                boxsize = boxsize[i],
                exhaustiveness = exhaustiveness,
                ncores = ncores,
                num_modes = num_modes,
                constraint = constraint,
                constraint_type = constraint_type,
                constraint_ref = constraint_ref,
                constraint_receptor_pdb_path = constraint_receptor_pdb_path[i],
                constraint_num_conf = constraint_num_conf,
                constraint_minimum_conf_rms = constraint_minimum_conf_rms,
            )
        else:
            vina_score, pdbqt = vinadock(
                Individual = Individual,
                wd = wd,
                vina_executable = vina_executable,
                receptor_pdbqt_path =  receptor_pdbqt_path[i],
                boxcenter = boxcenter[i],
                boxsize = boxsize[i],
                exhaustiveness = exhaustiveness,
                ncores = ncores,
                num_modes = num_modes,
            )
        Individual.vina_score.append(vina_score)
        pdbqt_list.append(pdbqt)
    # Update the pdbqt attribute
    Individual.pdbqt = pdbqt_list

    # Initialize base and exponent
    base = 1
    exponent = 0
    # Run for vina_scores
    for vs, vst in zip(Individual.vina_score, vina_score_type):
        for key in desirability[vst]:
            if key == 'w':
                w = desirability[vst][key]
            elif key in utils.DerringerSuichDesirability():
                d = utils.DerringerSuichDesirability()[key](vs, **desirability[vst][key])
            else:
                raise RuntimeError(f"Inside the desirability dictionary you provided for the variable = vina_scores[{vst}] a non implemented key = {key}. Only are possible: 'w' (standing for weight) and any possible Derringer-Suich desirability function: {utils.DerringerSuichDesirability().keys()}.")
        base *= d**w
        exponent += w

    # Average
    #D = (w_qed*d_qed + w_sa_score*d_sa_score + sum([w_vina_score*d_vina_score for w_vina_score, d_vina_score in zip(w_vina_scores, d_vina_scores)])) / (w_qed + w_sa_score + sum(w_vina_scores))
    # Geometric mean
    # D = (d_qed**w_qed * d_sa_score**w_sa_score * np.prod([d_vina_score**w_vina_score for d_vina_score, w_vina_score in zip(d_vina_scores, w_vina_scores)]))** (1/(w_qed + w_sa_score + sum(w_vina_scores)))
    D = base**(1/exponent)
    # And because we are minimizing we have to return
    Individual.cost = 1 - D
    return Individual


if __name__ == '__main__':
    pass