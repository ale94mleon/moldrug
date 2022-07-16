#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from moldrug import utils
from rdkit.Chem import QED
import os, numpy as np



def Cost(Individual:utils.Individual, wd:str = '.vina_jobs', vina_executable:str = 'vina', receptor_path:str = None, boxcenter:list[float] = None, boxsize:list[float] =None, exhaustiveness:int = 8, ncores:int = 1,  num_modes:int = 1):
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
        Where the receptor.pdbqt is located, by default None
    boxcenter : list[float], optional
        A list of three floats with the definition of the center of the box for docking, by default None
    boxsize : list[float], optional
        A list of three floats with the definition of the box size in angstrom of the docking box (x, y, z), by default None
    exhaustiveness : int, optional
        Parameter of vina that controls the accuracy of the searching on Vina, by default 8
    ncores : int, optional
        Number of cpus to use in Vina, by default 1
    num_modes : int, optional
        How many modes should Vina export, by default 1

    Returns
    -------
    utils.Individual
        A new instance of the original Individual with the the new attributes: pdbqt, qed, vina_score, sa_score and cost
    
    Example
    -------
    .. ipython:: python
        :okwarning:
        :okexcept:

        from moldrug import utils, fitness
        import tempfile, os
        from moldrug.data import ligands, boxes, receptors
        tmp_path = tempfile.TemporaryDirectory()
        ligand_smiles = ligands.r_x0161
        I = utils.Individual(ligand_smiles)
        print(type(I))
        print(ligand_smiles)
        print(I.mol)
        print(I.pdbqt)
        receptor_path = os.path.join(tmp_path.name,'receptor.pdbqt')
        with open(receptor_path, 'w') as r: r.write(receptors.r_x0161)
        box = boxes.r_x0161['A']
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
    utils.run(cmd)

    # Getting the information
    best_energy = utils.VINA_OUT(os.path.join(wd, f'{Individual.idx}_out.pdbqt')).BestEnergy()
    # Changing the xyz conformation by the conformation of the binding pose
    Individual.pdbqt = ''.join(best_energy.chunk)

    # Getting the Scoring function of Vina
    Individual.vina_score = best_energy.freeEnergy

    # Adding the cost using all the information of qed, sas and vina_cost


    # Construct the desirability
    # Quantitative estimation of drug-likness (ranges from 0 to 1). We could use just the value perse, but using LargerTheBest we are more permissible. 
    w_qed = 1
    d_qed = utils.LargerTheBest(Individual.qed, LowerLimit=0.1, Target=0.75, r = 1)
    # Synthetic accessibility score (ranges from 1 to 10). We would like to have 3 or lower and beyond 7 we assume that is not good.
    w_sa_score = 1
    d_sa_score = utils.SmallerTheBest(Individual.sa_score, Target = 3, UpperLimit=7, r = 1)
    # Vina. In this case is difficult to create a desirability because the range, we will use as target -12 upperlimit -6.
    w_vina_score = 1
    d_vina_score = utils.SmallerTheBest(Individual.vina_score, Target = -12, UpperLimit=-6, r = 1)
    # Average 
    #D = (w_qed*d_qed + w_sa_score*d_sa_score + w_vina_score*d_vina_score) / (w_qed + w_sa_score + w_vina_score)
    # Geometric mean
    D = (d_qed**w_qed * d_sa_score**w_sa_score * d_vina_score**w_vina_score)**(1/(w_qed + w_sa_score + w_vina_score))
    # And because we are minimizing we have to return 
    Individual.cost = 1 - D
    return Individual

def CostMultiReceptors(Individual:utils.Individual, wd:str = '.vina_jobs', vina_executable = 'vina', receptor_path:list = None, boxcenter:list = None, vina_score_types:list = None, boxsize:list =None, exhaustiveness:int = 8, ncores:int = 1,  num_modes:int = 1):
    """A fitness function that add multiple receptor on the responses variables.

    Args:
        Individual (utils.Individual): An Individual Instance
        wd (str, optional): The path were all the vina jobs will be executed. Defaults to '.vina_jobs'.
        receptor_path (list, optional): A list of paths of the corresponded pdbqt files. Defaults to None.
        boxcenter (list, optional): A list of list with the values x, y, z of the center. Defaults to None.
        vina_score_type (list, optional): A list of strings: 'min' or 'max'. For example if we are using two receptors: the first one we are looking for a minimum of the vina scoring function and in the other one a maximum (we are looking for specificity) then you should provided ['min', 'max']. Defaults to None.
        boxsize (list, optional): A list of list with the values x, y, z of the for the size. Defaults to None.. Defaults to None.
        exhaustiveness (int, optional): The exhaustiveness of Vina. Defaults to 8.
        ncores (int, optional): Number of cores to use on Vina. Defaults to 1.
        num_modes (int, optional): Number of modes to return from the Vina simulation. Defaults to 1.

    Returns:
        Individual: A modified Individual instance with the added attributes: qed, sa_score, pdbqt, vina_score and cost
    """
    sascorer = utils.import_sascorer()
    Individual.qed = QED.weights_mean(Individual.mol)
    
    # Getting synthetic accessibility score
    Individual.sa_score = sascorer.calculateScore(Individual.mol)
    
    # Getting Vina score
    Individual.vina_pdbqts = []
    Individual.vina_scores = []
    for i in range(len(receptor_path)):
        cmd = f"{vina_executable} --receptor {receptor_path[i]} --ligand {os.path.join(wd, f'{Individual.idx}_{i}.pdbqt')} "\
            f"--center_x {boxcenter[i][0]} --center_y {boxcenter[i][1]} --center_z {boxcenter[i][2]} "\
            f"--size_x {boxsize[i][0]} --size_y {boxsize[i][1]} --size_z {boxsize[i][2]} "\
            f"--out {os.path.join(wd, f'{Individual.idx}_{i}_out.pdbqt')} --cpu {ncores} --exhaustiveness {exhaustiveness} --num_modes {num_modes}"
        #print(cmd)
        # Creating the ligand pdbqt
        with open(os.path.join(wd, f'{Individual.idx}_{i}.pdbqt'), 'w') as l:
            l.write(Individual.pdbqt)
        utils.run(cmd)

        # Getting the information
        best_energy = utils.VINA_OUT(os.path.join(wd, f'{Individual.idx}_{i}_out.pdbqt')).BestEnergy()
        # Changing the xyz conformation by the conformation of the binding pose
        Individual.vina_pdbqts.append(''.join(best_energy.chunk))

        # Getting the Scoring function of Vina
        Individual.vina_scores.append(best_energy.freeEnergy)

    # Adding the cost using all the information of qed, sas and vina_cost


    # Construct the desirability
    # Quantitative estimation of drug-likness (ranges from 0 to 1). We could use just the value perse, but using LargerTheBest we are more permissible. 
    w_qed = 1
    d_qed = utils.LargerTheBest(Individual.qed, LowerLimit=0.1, Target=0.75, r = 1)
    # Synthetic accessibility score (ranges from 1 to 10). We would like to have 3 or lower and beyond 7 we assume that is not good.
    w_sa_score = 1
    d_sa_score = utils.SmallerTheBest(Individual.sa_score, Target = 3, UpperLimit=7, r = 1)
    # Vina. In this case is difficult to create a desirability because the range, we will use as target -12 upperlimit -6.
    w_vina_scores = len(vina_score_types)*[1]
    d_vina_scores = []
    for vina_score, vina_score_type in zip(Individual.vina_scores, vina_score_types):
        if vina_score_type == 'min':
            d_vina_scores.append(utils.SmallerTheBest(vina_score, Target = -12, UpperLimit=-6, r = 1))
        elif vina_score_type == 'max':
            d_vina_scores.append(utils.LargerTheBest(vina_score, LowerLimit = -4, Target = 0, r = 1))
        else:
            # I have to check first if the list vina_score_types only contains min and max strings
            pass
    # Average 
    #D = (w_qed*d_qed + w_sa_score*d_sa_score + sum([w_vina_score*d_vina_score for w_vina_score, d_vina_score in zip(w_vina_scores, d_vina_scores)])) / (w_qed + w_sa_score + sum(w_vina_scores))
    # Geometric mean
    D = (d_qed**w_qed * d_sa_score**w_sa_score * np.prod([d_vina_score**w_vina_score for d_vina_score, w_vina_score in zip(d_vina_scores, w_vina_scores)]))** (1/(w_qed + w_sa_score + sum(w_vina_scores)))
    # And because we are minimizing we have to return 
    Individual.cost = 1 - D
    return Individual
if __name__ == '__main__':
    pass