#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Clean the path and all subdirectories and create the distribution plot.
"""
import os
import seaborn as sns
import matplotlib.pyplot as plt
from typing import List
from moldrug.utils import Individual, decompress_pickle, tar_errors
import numpy as np
from moldrug import utils
import pandas as pd
import seaborn as sns
import argparse

# This function is possible from version 2.2.0

def plot_dist(individuals:List[Individual], properties:List[str], every_gen:int = 1):
    """Create the violin plot for the MolDrug run

    Parameters
    ----------
    individuals : List[Individual]
        A list of individuals
    properties : List[str]
        A list of the properties to be graph (must be attributes of the provided individuals)
    every_gen : int, optional
        Frequency to plot the distrubution: every how many generations, by default 1

    Returns
    -------
    tuple
        fig, axes
    """


    # Set up the matplotlib figure
    sns.set_theme(style="whitegrid")
    fig, axes = plt.subplots(nrows = len(properties), figsize=(25, 25))

    SawIndividuals = utils.to_dataframe(individuals).drop(['pdbqt'], axis = 1).replace([np.inf, -np.inf], np.nan).dropna()
    SawIndividuals = SawIndividuals[SawIndividuals['kept_gens'].map(len) != 0].reset_index(drop=True)
    gen_idxs = sorted(SawIndividuals.genID.unique())
    NumGens = max(gen_idxs)

    # Set pop to the initial population and pops out the first gen
    pop = SawIndividuals[SawIndividuals.genID == gen_idxs.pop(0)].sort_values(by=["cost"])
    pops = pop.copy()
    for gen_idx in gen_idxs:
        idx = [i for i in range(SawIndividuals.shape[0]) if gen_idx in SawIndividuals.loc[i,'kept_gens']]
        pop = SawIndividuals.copy().iloc[idx,:].assign(genID=gen_idx)
        pops = pd.concat([pops, pop.copy()])
    # Draw a violinplot with a narrow bandwidth than the default
    pops = pops.loc[pops['genID'].isin([gen for gen in range(0, NumGens+every_gen, every_gen)])]
    for i, prop in enumerate(properties):
        sns.violinplot(x = 'genID', y = prop, data=pops, palette="Set3", bw=.2, cut=0, linewidth=1, ax=axes[i])

    return fig, axes

def plot_dist_path(path:str, properties:List[str], every_gen:int = 1, pbz_file_name = '04_local_result.pbz2'):
    """Is a wraper arounf plot_dist()
    It will look in the root folder and all the subdirectories for pbz_file_name.
    If it is found plot.svg will be created on that subdirectory.

    Parameters
    ----------
    path : str
        The root path form where the pbz_file_name will be looking for
    properties : List[str]
         A list of the properties to be graph (must be attributes of the provided individuals)
    every_gen : int, optional
        Frequency to plot the distrubution: every how many generations, by default 1
    pbz_file_name : str, optional
        Name of the output binary result of MolDrug, by default '04_local_result.pbz2'
    """
    cwd = os.getcwd()
    for (root, _, filenames) in list(os.walk(path)): #dirpath, dirnames, filenames
        os.chdir(root)

        if pbz_file_name in filenames:
            print(root)
            try:
                SawIndividuals = decompress_pickle(pbz_file_name).SawIndividuals
            except Exception as e:
                print('No distribution.')
                continue
            # print(utils.to_dataframe(SawIndividuals).columns)
            # raise RuntimeError
            try:
                fig, _ = plot_dist(SawIndividuals, properties, every_gen=every_gen)
            except Exception:
                try:
                    fig, _ = plot_dist(SawIndividuals, ['vina_score', 'cost'], every_gen=every_gen)
                except Exception:
                    print('No distribution.')
            try:
                fig.savefig('dist.png')
                plt.close()
            except Exception:
                pass
    os.chdir(cwd)

# def cleanup_errors(path:str, error:str = '*error*'):
#     """It looks in the root path provided and look for the
#     files with the pattern of error in path and all it subdirectories.
#     If there are some of them; a directory error will be created and all
#     the error files will be modved there.

#     Parameters
#     ----------
#     path : str
#         Root path to look for error
#     error : str, optional
#         The pattern to look for the errors, by default '*error*'
#     """
#     cwd = os.getcwd()
#     for (root, _, _) in list(os.walk(path)):
#         os.chdir(root)
#         possible_error_files = [f for f in glob(error) if os.path.isfile(f)]
#         if possible_error_files and os.path.basename(root) != 'error':
#             os.makedirs('error', exist_ok=True)
#             print(root)
#             for possible_error_file in tqdm(possible_error_files):
#                 copy2(possible_error_file, 'error')
#                 os.remove(possible_error_file)
#     os.chdir(cwd)

def cleanup_errors(path:str, error_path:str = 'error'):
    """It looks in the root path provided and look for the
    files with the pattern of error in path and all it subdirectories.
    If there are some of them; a directory error will be created and all
    the error files will be modved there.

    Parameters
    ----------
    path : str
        Root path to look for error
    error_path : str, optional
        Where the errors are storged, by default 'error'
    """
    cwd = os.getcwd()
    for (root, _, _) in list(os.walk(path)):
        if os.path.basename(root) == error_path: continue
        os.chdir(root)
        tar_errors(error_path=error_path)

    os.chdir(cwd)

# def taka(path):
#     cwd = os.getcwd()
#     for (root, _, _) in list(os.walk(path)):
#         os.chdir(root)

#         if os.path.basename(root) == 'error':
#             print(root)
#             utils.run('cd ..; tar -czf error.tar.gz error/ && rm -r error/', executable='/bin/bash')
#     os.chdir(cwd)

# def taka2(path):
#     cwd = os.getcwd()
#     for (root, _, _) in list(os.walk(path)):
#         os.chdir(root)
#         try:
#             utils.run('rm dist.svg', executable='/bin/bash')
#         except:...
#     os.chdir(cwd)

def postpro_cmd():
    """CLI for plot_dist_path and cleanup_errors
    """
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        help='The root path form where the pbz_file_name will be looking for and the graph will be created',
        dest='path',
        type=str)
    parser.add_argument(
        '-p', '--properties',
        help='A list of the properties to be graph (must be attributes of the provided individuals), by default %(default)s',
        nargs='+',
        dest='properties',
        default=['qed', 'sa_score', 'vina_score', 'cost'],
        type=str)
    parser.add_argument(
        '--every_gen',
        help='Frequency to plot the distrubution: every how many generations, by default %(default)s',
        dest='every_gen',
        default=10,
        type=int)
    parser.add_argument(
        '--pbz_file_name',
        help='Name of the output binary result of MolDrug, by default %(default)s',
        dest='pbz_file_name',
        default='04_local_result.pbz2',
        type=str)
    parser.add_argument(
        '--error',
        help='Directory where the errors are storged, by default %(default)s',
        dest='error',
        default='*error*',
        type=str)
    args = parser.parse_args()

    # To clean the directories and subdirectories
    print('Cleaning directories...')
    cleanup_errors(
        path=args.path,
        error_path=args.error
    )
    # To create the distribution plot
    print('Creating distribution plots ...')
    plot_dist_path(
        path=args.path,
        properties=args.properties,
        every_gen=args.every_gen,
        pbz_file_name=args.pbz_file_name,
    )
    print('Done!')


if __name__ == '__main__':
    #/home/ale/mnt/snowden2/MolDrug/BI
    postpro_cmd()