#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from lead import vina, utils
from rdkit import Chem
import numpy as np
from rdkit.Chem import QED
import sys, os
from rdkit.Chem import RDConfig
sys.path.append(os.path.join(RDConfig.RDContribDir, 'SA_Score'))
# now you can import sascore!
import sascorer


def __VinaCost(Individual, wd = '.vina_jobs', receptor_path = None, boxcenter = None, boxsize =None, exhaustiveness = 8, vina_cpus = 1,  num_modes = 1):
    cmd = f"{vina.vina_executable} --receptor {receptor_path} --ligand {os.path.join(wd, f'{Individual.idx}.pdbqt')} "\
        f"--center_x {boxcenter[0]} --center_y {boxcenter[1]} --center_z {boxcenter[2]} "\
        f"--size_x {boxsize[0]} --size_y {boxsize[1]} --size_z {boxsize[2]} "\
        f"--out {os.path.join(wd, f'{Individual.idx}_out.pdbqt')} --cpu {vina_cpus} --exhaustiveness {exhaustiveness} --num_modes {num_modes}"
    #print(cmd)
    # Creating the ligand pdbqt
    with open(os.path.join(wd, f'{Individual.idx}.pdbqt'), 'w') as l:
        l.write(Individual.pdbqt)
    utils.run(cmd)

    # Getting the information
    best_energy = vina.VINA_OUT(os.path.join(wd, f'{Individual.idx}_out.pdbqt')).BestEnergy()
    # Changing the xyz conformation by the conformation of the binding pose
    Individual.pdbqt = best_energy.chunk

    # Getting the Scoring function of Vina
    Individual.vina_score = best_energy.freeEnergy
    Individual.cost = best_energy.freeEnergy

    # Adding the cost using all the information of qed, sas and vina_cost

    # Construct the desirability.
    # Merge the desirabilities
    # qed ranges between 0 and 1
    # sas 0 to 10
    # Vina is more difficult to capture But i already have an idea.
    # save cost
    # Think in a more clean way to this cost function
    # We could have vina but also a cost function that integrate all of them is cleaner 
    return Individual
def __VinaCostLipinski(Individual, wd = '.vina_jobs', receptor_path = None, boxcenter = None, boxsize =None, exhaustiveness = 8, vina_cpus = 1,  num_modes = 1):
    if utils.lipinski_filter(Individual.mol):
        cmd = f"{vina.vina_executable} --receptor {receptor_path} --ligand {os.path.join(wd, f'{Individual.idx}.pdbqt')} "\
            f"--center_x {boxcenter[0]} --center_y {boxcenter[1]} --center_z {boxcenter[2]} "\
            f"--size_x {boxsize[0]} --size_y {boxsize[1]} --size_z {boxsize[2]} "\
            f"--out {os.path.join(wd, f'{Individual.idx}_out.pdbqt')} --cpu {vina_cpus} --exhaustiveness {exhaustiveness} --num_modes {num_modes}"
        #print(cmd)
        # Creating the ligand pdbqt
        with open(os.path.join(wd, f'{Individual.idx}.pdbqt'), 'w') as l:
            l.write(Individual.pdbqt)
        utils.run(cmd)

        # Getting the information
        best_energy = vina.VINA_OUT(os.path.join(wd, f'{Individual.idx}_out.pdbqt')).BestEnergy()
        # Changing the xyz conformation by the conformation of the binding pose
        Individual.pdbqt = best_energy.chunk

        # Getting the Scoring function of Vina
        Individual.vina_score = best_energy.freeEnergy
        Individual.cost = best_energy.freeEnergy
    else:
        Individual.vina_score = np.inf
        Individual.cost = np.inf


    # Adding the cost using all the information of qed, sas and vina_cost

    # Construct the desirability.
    # Merge the desirabilities
    # qed ranges between 0 and 1
    # sas 0 to 10
    # Vina is more difficult to capture But i already have an idea.
    # save cost
    # Think in a more clean way to this cost function
    # We could have vina but also a cost function that integrate all of them is cleaner 
    return Individual

def QED_SAS_vina_cost(Individual, wd = '.vina_jobs', receptor_path = None, boxcenter = None, boxsize =None, exhaustiveness = 8, vina_cpus = 1,  num_modes = 1):
    # multicriteria optimization,Optimization of Several Response Variables
    # Getting estimate of drug-likness
    Individual.qed = QED.weights_mean(Individual.mol)
    
    # Getting synthetic accessibility score
    Individual.sa_score = sascorer.calculateScore(Individual.mol)
    
    if Individual.qed >=0.75 and Individual.sa_score <= 4:

        # Getting Vina score
        cmd = f"{vina.vina_executable} --receptor {receptor_path} --ligand {os.path.join(wd, f'{Individual.idx}.pdbqt')} "\
            f"--center_x {boxcenter[0]} --center_y {boxcenter[1]} --center_z {boxcenter[2]} "\
            f"--size_x {boxsize[0]} --size_y {boxsize[1]} --size_z {boxsize[2]} "\
            f"--out {os.path.join(wd, f'{Individual.idx}_out.pdbqt')} --cpu {vina_cpus} --exhaustiveness {exhaustiveness} --num_modes {num_modes}"
        #print(cmd)
        # Creating the ligand pdbqt
        with open(os.path.join(wd, f'{Individual.idx}.pdbqt'), 'w') as l:
            l.write(Individual.pdbqt)
        utils.run(cmd)

        # Getting the information
        best_energy = vina.VINA_OUT(os.path.join(wd, f'{Individual.idx}_out.pdbqt')).BestEnergy()
        # Changing the xyz conformation by the conformation of the binding pose
        Individual.pdbqt = best_energy.chunk

        # Getting the Scoring function of Vina
        Individual.vina_score = best_energy.freeEnergy
        Individual.cost = best_energy.freeEnergy
    else:
        # Getting the Scoring function of Vina
        Individual.vina_score =np.inf
        Individual.cost = np.inf
    return Individual


def Cost(Individual, wd = '.vina_jobs', receptor_path = None, boxcenter = None, boxsize =None, exhaustiveness = 8, vina_cpus = 1,  num_modes = 1):
    # multicriteria optimization,Optimization of Several Response Variables
    # Getting estimate of drug-likness
    Individual.qed = QED.weights_mean(Individual.mol)
    
    # Getting synthetic accessibility score
    Individual.sa_score = sascorer.calculateScore(Individual.mol)
    
    # Getting Vina score
    cmd = f"{vina.vina_executable} --receptor {receptor_path} --ligand {os.path.join(wd, f'{Individual.idx}.pdbqt')} "\
        f"--center_x {boxcenter[0]} --center_y {boxcenter[1]} --center_z {boxcenter[2]} "\
        f"--size_x {boxsize[0]} --size_y {boxsize[1]} --size_z {boxsize[2]} "\
        f"--out {os.path.join(wd, f'{Individual.idx}_out.pdbqt')} --cpu {vina_cpus} --exhaustiveness {exhaustiveness} --num_modes {num_modes}"
    #print(cmd)
    # Creating the ligand pdbqt
    with open(os.path.join(wd, f'{Individual.idx}.pdbqt'), 'w') as l:
        l.write(Individual.pdbqt)
    utils.run(cmd)

    # Getting the information
    best_energy = vina.VINA_OUT(os.path.join(wd, f'{Individual.idx}_out.pdbqt')).BestEnergy()
    # Changing the xyz conformation by the conformation of the binding pose
    Individual.pdbqt = best_energy.chunk

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
    D = (d_qed**w_qed * d_sa_score**w_sa_score * d_vina_score**w_vina_score)** 1/(w_qed + w_sa_score + w_vina_score)
    # And because we are minimizing we have to return 
    Individual.cost = 1 - D
    return Individual
