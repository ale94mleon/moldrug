#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from random import expovariate
from moldrug import utils
from rdkit.Chem import QED
import os, numpy as np
from typing import Dict, List



def Cost(
    Individual:utils.Individual,
    wd:str = '.vina_jobs',
    vina_executable:str = 'vina',
    receptor_path:str = None,
    boxcenter:List[float] = None,
    boxsize:List[float] =None,
    exhaustiveness:int = 8,
    ncores:int = 1, 
    num_modes:int = 1,
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
    """This is the main Cost function of the module

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
    desirability : dict, optional
        The definition of the desirability to use for each used variable = [qed, sa_score, vina_score].
        Each variable only will accept the keys [w, and the name of the desirability function of :meth:`moldrug.utils.DerringerSuichDesirability`]
        ,by default { 'qed': { 'w': 1, 'LargerTheBest': { 'LowerLimit': 0.1, 'Target': 0.75, 'r': 1 } }, 'sa_score': { 'w': 1, 'SmallerTheBest': { 'Target': 3, 'UpperLimit': 7, 'r': 1 } }, 'vina_score': { 'w': 1, 'SmallerTheBest': { 'Target': -12, 'UpperLimit': -6, 'r': 1 } } }

    Returns
    -------
    utils.Individual
        A new instance of the original Individual with the the new attributes: pdbqt, qed, vina_score, sa_score and cost
    Example
    -------
    .. ipython:: python
        from moldrug import utils, fitness
        import tempfile, os
        from moldrug.data import ligands, boxes, receptors
        tmp_path = tempfile.TemporaryDirectory()
        ligand_smiles = ligands.r_x0161
        I = utils.Individual(ligand_smiles)
        receptor_path = os.path.join(tmp_path.name,'receptor.pdbqt')
        with open(receptor_path, 'w') as r: r.write(receptors.r_x0161)
        box = boxes.r_x0161['A']
        # Using the default desirability
        NewI = fitness.Cost(Individual = I,wd = tmp_path.name,receptor_path = receptor_path,boxcenter = box['boxcenter'],boxsize = box['boxsize'],exhaustiveness = 4,ncores = 4)
        print(NewI.cost, NewI.vina_score, NewI.qed, NewI.sa_score)
    """
    sascorer = utils.import_sascorer()
    # multicriteria optimization,Optimization of Several Response Variables
    # Getting estimate of drug-likness
    Individual.qed = QED.weights_mean(Individual.mol)

    # Getting synthetic accessibility score
    Individual.sa_score = sascorer.calculateScore(Individual.mol)

    # Getting Vina score
    cmd = f"{vina_executable} --receptor {receptor_path} --ligand {os.path.join(wd, f'{Individual.idx}.pdbqt')} "\
        f"--center_x {boxcenter[0]} --center_y {boxcenter[1]} --center_z {boxcenter[2]} "\
        f"--size_x {boxsize[0]} --size_y {boxsize[1]} --size_z {boxsize[2]} "\
        f"--out {os.path.join(wd, f'{Individual.idx}_out.pdbqt')} --cpu {ncores} --exhaustiveness {exhaustiveness} --num_modes {num_modes}"
    #print(cmd)
    # Creating the ligand pdbqt
    with open(os.path.join(wd, f'{Individual.idx}.pdbqt'), 'w') as l:
        l.write(Individual.pdbqt)
    # try:
    #     # This should be always the case. Consider to use some atom filter in order to eliminate incompatible vina molecules
    utils.run(cmd)

    # Getting the information
    best_energy = utils.VINA_OUT(os.path.join(wd, f'{Individual.idx}_out.pdbqt')).BestEnergy()
    # Changing the 3D conformation by the conformation of the binding pose
    Individual.pdbqt = ''.join(best_energy.chunk)

    # Getting the Scoring function of Vina
    Individual.vina_score = best_energy.freeEnergy
# except:
    #     # If something happens with vina, the pdbqt is wrong or any other error. This could be dangerous because hide some possibles user errors. Think about it.
    #     pass

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

def CostMultiReceptors(
    Individual:utils.Individual,
    wd:str = '.vina_jobs',
    vina_executable:str = 'vina',
    receptor_paths:List[str] = None,
    vina_score_types:List[str] = None,
    boxcenters:List[float] = None,
    boxsizes:List[float] =None,
    exhaustiveness:int = 8,
    ncores:int = 1,
    num_modes:int = 1,
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
    """This function is similar to :meth:`moldrug.fitness.Cost` but it will add the possibility to work with more than one receptor.

    Parameters
    ----------
    Individual : utils.Individual
        A Individual with the pdbqt attribute
    wd : str, optional
        The working directory to execute the docking jobs, by default '.vina_jobs'
    vina_executable : str, optional
        This is the name of the vina executable, could be a path to the binary object (x, y, z), by default 'vina'
    receptor_paths : list[str], optional
        A list of location of the receptors pdbqt files, by default None
    vina_score_types : list[str], optional
        This is a list with the keywords 'min' and/or 'max'. E.g. If two receptor were provided and for the first one we would like to find a minimum in the vina scoring function and for the other one a maximum (selectivity for the first receptor); we must provided the list: ['min', 'max'], by default None
    boxcenters : list[float], optional
        A list of three floats with the definition of the center of the box in angstrom for docking (x, y, z), by default None
    boxsizes : list[float], optional
        A list of three floats with the definition of the box size in angstrom of the docking box (x, y, z), by default None
    exhaustiveness : int, optional
        Parameter of vina that controls the accuracy of the docking searching, by default 8
    ncores : int, optional
         Number of cpus to use in Vina, by default 1
    num_modes : int, optional
        How many modes should Vina export, by default 1
    desirability : dict, optional
        The definition of the desirability to use for each used variable = [qed, sa_score, vina_scores].
        Each variable only will accept the keys [w, and the name of the desirability function of :meth:`moldrug.utils.DerringerSuichDesirability`].
        In the case of vina_scores there is another layer for the vina_score_type= [min, max],
        by default { 'qed': { 'w': 1, 'LargerTheBest': { 'LowerLimit': 0.1, 'Target': 0.75, 'r': 1 } }, 'sa_score': { 'w': 1, 'SmallerTheBest': { 'Target': 3, 'UpperLimit': 7, 'r': 1 } }, 'vina_score': { 'min':{ 'w': 1, 'SmallerTheBest': { 'Target': -12, 'UpperLimit': -6, 'r': 1 } }, 'max':{ 'w': 1, 'LargerTheBest': { 'LowerLimit': -4, 'Target': 0, 'r': 1 } } } }

    Returns
    -------
    utils.Individual
        A new instance of the original Individual with the the new attributes: pdbqts [a list of pdbqt], qed, vina_scores [a list of vina_score], sa_score and cost

    Example
    -------
    .. ipython:: python

        from moldrug import utils, fitness
        import tempfile, os
        from moldrug.data import ligands, boxes, receptors
        tmp_path = tempfile.TemporaryDirectory()
        ligand_smiles = ligands.r_x0161
        I = utils.Individual(ligand_smiles)
        receptor_paths = [os.path.join(tmp_path.name,'receptor1.pdbqt'),os.path.join(tmp_path.name,'receptor2.pdbqt')]
        with open(receptor_paths[0], 'w') as r: r.write(receptors.r_x0161)
        with open(receptor_paths[1], 'w') as r: r.write(receptors.r_6lu7)
        boxcenters = [boxes.r_x0161['A']['boxcenter'], boxes.r_6lu7['A']['boxcenter']]
        boxsizes = [boxes.r_x0161['A']['boxsize'], boxes.r_6lu7['A']['boxsize']]
        vina_score_types = ['min', 'max']
        # Using the default desirability
        NewI = fitness.CostMultiReceptors(Individual = I,wd = tmp_path.name,receptor_paths = receptor_paths, vina_score_types = vina_score_types, boxcenters = boxcenters,boxsizes = boxsizes,exhaustiveness = 4,ncores = 4)
        print(NewI.cost, NewI.vina_scores, NewI.qed, NewI.sa_score)
    """
    sascorer = utils.import_sascorer()
    Individual.qed = QED.weights_mean(Individual.mol)

    # Getting synthetic accessibility score
    Individual.sa_score = sascorer.calculateScore(Individual.mol)

    # Getting Vina score
    Individual.pdbqts = []
    Individual.vina_scores = []
    for i in range(len(receptor_paths)):
        cmd = f"{vina_executable} --receptor {receptor_paths[i]} --ligand {os.path.join(wd, f'{Individual.idx}_{i}.pdbqt')} "\
            f"--center_x {boxcenters[i][0]} --center_y {boxcenters[i][1]} --center_z {boxcenters[i][2]} "\
            f"--size_x {boxsizes[i][0]} --size_y {boxsizes[i][1]} --size_z {boxsizes[i][2]} "\
            f"--out {os.path.join(wd, f'{Individual.idx}_{i}_out.pdbqt')} --cpu {ncores} --exhaustiveness {exhaustiveness} --num_modes {num_modes}"
        #print(cmd)
        # Creating the ligand pdbqt
        with open(os.path.join(wd, f'{Individual.idx}_{i}.pdbqt'), 'w') as l:
            l.write(Individual.pdbqt)
        utils.run(cmd)

        # Getting the information
        best_energy = utils.VINA_OUT(os.path.join(wd, f'{Individual.idx}_{i}_out.pdbqt')).BestEnergy()
        # Changing the 3D conformation by the conformation of the binding pose
        Individual.pdbqts.append(''.join(best_energy.chunk))

        # Getting the Scoring function of Vina
        Individual.vina_scores.append(best_energy.freeEnergy)

    # pops the region of vina_scores
    vina_desirability_section = desirability.pop('vina_scores')
    # Initialize base and exponent
    base = 1
    exponent = 0
    # Runs for all properties different to vina_scores
    for variable in desirability:
        for key in desirability[variable]:
            if key == 'w':
                w = desirability[variable][key]
            elif key in utils.DerringerSuichDesirability():
                d = utils.DerringerSuichDesirability()[key](getattr(Individual, variable), **desirability[variable][key])
            else:
                raise RuntimeError(f"Inside the desirability dictionary you provided for the variable = {variable} a non implemented key = {key}. Only are possible: 'w' (standing for weight) and any possible Derringer-Suich desirability function: {utils.DerringerSuichDesirability().keys()}. Only in the case of vina_scores [min and max] keys")
        base *= d**w
        exponent += w

    # Run only for vina_scores
    for vina_score, vina_score_type in zip(Individual.vina_scores, vina_score_types):
        for key in vina_desirability_section[vina_score_type]:
            if key == 'w':
                w = vina_desirability_section[vina_score_type][key]
            elif key in utils.DerringerSuichDesirability():
                d = utils.DerringerSuichDesirability()[key](vina_score, **vina_desirability_section[vina_score_type][key])
            else:
                RuntimeError(f"Inside the desirability dictionary you provided for the variable = vina_scores[{vina_score_type}] a non implemented key = {key}. Only are possible: 'w' (standing for weight) and any possible Derringer-Suich desirability function: {utils.DerringerSuichDesirability().keys()}.")
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