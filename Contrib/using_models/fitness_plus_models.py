#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from copy import deepcopy
from typing import Dict, List

import joblib
import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Descriptors

from moldrug import utils
from moldrug.fitness import _vinadock

RDLogger.DisableLog('rdApp.*')


class Featurizer:
    def __init__(self, fpb=2048, calc_descriptors=True, scale = False):
        """
        :param fpb: number of Morgan bits
        :param descriptors: boolean - calculate descriptors.
        """
        self.descriptor_list = ['MaxEStateIndex', 'MinEStateIndex', 'MaxAbsEStateIndex', 'MinAbsEStateIndex', 'qed',
                                'MolWt', 'HeavyAtomMolWt', 'ExactMolWt', 'NumValenceElectrons', 'NumRadicalElectrons',
                                'MaxPartialCharge', 'MinPartialCharge', 'MaxAbsPartialCharge', 'MinAbsPartialCharge',
                                'FpDensityMorgan1', 'FpDensityMorgan2', 'FpDensityMorgan3', 'BalabanJ', 'BertzCT',
                                'Chi0', 'Chi0n', 'Chi0v', 'Chi1', 'Chi1n', 'Chi1v', 'Chi2n', 'Chi2v', 'Chi3n', 'Chi3v',
                                'Chi4n', 'Chi4v', 'HallKierAlpha', 'Ipc', 'Kappa1', 'Kappa2', 'Kappa3', 'LabuteASA',
                                'PEOE_VSA1', 'PEOE_VSA10', 'PEOE_VSA11', 'PEOE_VSA12', 'PEOE_VSA13', 'PEOE_VSA14',
                                'PEOE_VSA2', 'PEOE_VSA3', 'PEOE_VSA4', 'PEOE_VSA5', 'PEOE_VSA6', 'PEOE_VSA7',
                                'PEOE_VSA8', 'PEOE_VSA9', 'SMR_VSA1', 'SMR_VSA10', 'SMR_VSA2', 'SMR_VSA3', 'SMR_VSA4',
                                'SMR_VSA5', 'SMR_VSA6', 'SMR_VSA7', 'SMR_VSA8', 'SMR_VSA9', 'SlogP_VSA1', 'SlogP_VSA10',
                                'SlogP_VSA11', 'SlogP_VSA12', 'SlogP_VSA2', 'SlogP_VSA3', 'SlogP_VSA4', 'SlogP_VSA5',
                                'SlogP_VSA6', 'SlogP_VSA7', 'SlogP_VSA8', 'SlogP_VSA9', 'TPSA', 'EState_VSA1',
                                'EState_VSA10', 'EState_VSA11', 'EState_VSA2', 'EState_VSA3', 'EState_VSA4',
                                'EState_VSA5', 'EState_VSA6', 'EState_VSA7', 'EState_VSA8', 'EState_VSA9',
                                'VSA_EState1', 'VSA_EState10', 'VSA_EState2', 'VSA_EState3', 'VSA_EState4',
                                'VSA_EState5', 'VSA_EState6', 'VSA_EState7', 'VSA_EState8', 'VSA_EState9',
                                'FractionCSP3', 'HeavyAtomCount', 'NHOHCount', 'NOCount', 'NumAliphaticCarbocycles',
                                'NumAliphaticHeterocycles', 'NumAliphaticRings', 'NumAromaticCarbocycles',
                                'NumAromaticHeterocycles', 'NumAromaticRings', 'NumHAcceptors', 'NumHDonors',
                                'NumHeteroatoms', 'NumRotatableBonds', 'NumSaturatedCarbocycles',
                                'NumSaturatedHeterocycles', 'NumSaturatedRings', 'RingCount', 'MolLogP', 'MolMR',
                                'fr_Al_COO', 'fr_Al_OH', 'fr_Al_OH_noTert', 'fr_ArN', 'fr_Ar_COO', 'fr_Ar_N',
                                'fr_Ar_NH', 'fr_Ar_OH', 'fr_COO', 'fr_COO2', 'fr_C_O', 'fr_C_O_noCOO', 'fr_C_S',
                                'fr_HOCCN', 'fr_Imine', 'fr_NH0', 'fr_NH1', 'fr_NH2', 'fr_N_O', 'fr_Ndealkylation1',
                                'fr_Ndealkylation2', 'fr_Nhpyrrole', 'fr_SH', 'fr_aldehyde', 'fr_alkyl_carbamate',
                                'fr_alkyl_halide', 'fr_allylic_oxid', 'fr_amide', 'fr_amidine', 'fr_aniline',
                                'fr_aryl_methyl', 'fr_azide', 'fr_azo', 'fr_barbitur', 'fr_benzene',
                                'fr_benzodiazepine', 'fr_bicyclic', 'fr_diazo', 'fr_dihydropyridine', 'fr_epoxide',
                                'fr_ester', 'fr_ether', 'fr_furan', 'fr_guanido', 'fr_halogen', 'fr_hdrzine',
                                'fr_hdrzone', 'fr_imidazole', 'fr_imide', 'fr_isocyan', 'fr_isothiocyan', 'fr_ketone',
                                'fr_ketone_Topliss', 'fr_lactam', 'fr_lactone', 'fr_methoxy', 'fr_morpholine',
                                'fr_nitrile', 'fr_nitro', 'fr_nitro_arom', 'fr_nitro_arom_nonortho', 'fr_nitroso',
                                'fr_oxazole', 'fr_oxime', 'fr_para_hydroxylation', 'fr_phenol',
                                'fr_phenol_noOrthoHbond', 'fr_phos_acid', 'fr_phos_ester', 'fr_piperdine',
                                'fr_piperzine', 'fr_priamide', 'fr_prisulfonamd', 'fr_pyridine', 'fr_quatN',
                                'fr_sulfide', 'fr_sulfonamd', 'fr_sulfone', 'fr_term_acetylene', 'fr_tetrazole',
                                'fr_thiazole', 'fr_thiocyan', 'fr_thiophene', 'fr_unbrch_alkane', 'fr_urea']
        self.fpb = fpb
        self.calc_descriptors = calc_descriptors

    def featurize(self, mol):
        """
        Returns features of your molecule.
        :param mol: rdkit molecule object
        :return: np.array
        """
        if self.fpb:
            fp_arr = np.array(AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=self.fpb))
            if not self.calc_descriptors:
                return fp_arr

        desc_arr = np.array([getattr(Descriptors, desc)(mol) for desc in self.descriptor_list])
        if self.fpb:
            return np.append(fp_arr, desc_arr)
        else:
            return desc_arr

    def __call__(self, mol):
        """
        Returns features of your molecule.
        :param mol: rdkit molecule object
        :return: np.array
        """
        return self.featurize(mol)


class Predictors():
    def __init__(self, models, featurizer):
        """

        :param models: list[str], locations of trained predictive models
        :param featurizer: Featurizer object
        """
        self.models = [joblib.load(model) for model in models]
        self.featurizer = featurizer

    def predict(self, molecule):
        """
        Based on molecule calculates prediction for each model
        :param molecule: rdkit.molecule
        :return:
        """
        X = self.featurizer.featurize(molecule)
        to_return = []
        for item in self.models:
            # To avoid inplace modifications:
            model = item['model']
            scaler = item['scaler']
            if scaler:
                X_to_use = scaler.transform(np.copy(X).reshape(1, -1))
                X_to_use = X_to_use.flatten()
                prediction = model.predict(X_to_use[np.newaxis,:])[0]
            else:
                prediction = model.predict(X[np.newaxis,:])[0]
            to_return.append(prediction)

        return to_return

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
    models:Dict = None,
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
        Maximum number of conformer to be generated internally by moldrug , by default 100
    constraint_minimum_conf_rms : int, optional
        RMS to filter duplicate conformers, by default 0.01
    models : dict = None,
        The definition of the models,
        by default None which means that it will be used:
        models = {
        'hppb': 'hppb.jlib',
        'clearance': 'clearance.jlib',
        }
        You could define as many models as you want, but the desirability must be defined accordantly. The key values
        must be the same as in desirability. In case that provided more models than hppb and clearance
        the same must be done in desirability.
    desirability : dict, optional
        The definition of the desirability to use for each used variable = [qed, sa_score, vina_score].
        Each variable only will accept the keys [w, and the name of the desirability function of :meth:`moldrug.utils.DerringerSuichDesirability`],
        by default None which means that it will be used:
        desirability = {
        'hppb': {'w': 1,'LargerTheBest': {'LowerLimit': 25,'Target': 75, 'r': 1}
        },
        'clearance': {'w': 1,'SmallerTheBest': {'Target': 20,'UpperLimit': 125,'r': 1}
        },
        'vina_score': {'w': 1,'SmallerTheBest': {'Target': -12,'UpperLimit': -6,'r': 1}
        }
        }
        If vina_score is not provided. The default values will be set anyway.
        The date set used have the following limits hppb = [11.18; 99.95] (higher is better);
        clearance = [3.0; 150.0] (lower is better). But your models could change.
    Returns
    -------
    utils.Individual
        A new instance of the original Individual with the the new attributes:
        pdbqt, qed, vina_score, sa_score and cost.
        cost attribute will be a number between 0 and 1, been 0 the optimal value.
    """
    if not models:
        models = {
            'hppb':  'hppb.jlib',
            'clearance': 'clearance.jlib',
        }
    if not desirability:
        desirability = {
            'hppb': {
                'w': 1,
                'LargerTheBest': {
                    'LowerLimit': 25,
                    'Target':75,
                    'r': 1
                }
            },
            'clearance': {
                'w': 1,
                'SmallerTheBest': {
                    'Target': 20,
                    'UpperLimit': 125,
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
    # Check that everything is ok with naming in models and desirability
    diff = list(set(desirability) - set(models))
    if len(diff) == 0:
        # vina_score was not defined. the default values will be used.
        desirability['vina_score'] = {
            'w': 1,
            'SmallerTheBest': {
                'Target': -12,
                'UpperLimit': -6,
                'r': 1
            }
        }
    else:
        if len(diff) != 1:
            raise Exception(f"You provided models = {models.keys()} and desirability = {desirability.keys()}. "\
                "However, desirability must have the same keywords as models and optionally the keyword vina_score")
        elif diff[0] != 'vina_score':
            raise Exception(f"desirability has the keyword: '{diff[0]}' which is not defined in models =  {models.keys()} and is not vina_score")

    # Getting and setting properties on the Individual
    predictor = Predictors(featurizer=Featurizer(), models=models.values())

    # MUST be a copy of Individual.mol becasue if not creazy stuffs will happen!!
    predictions = predictor.predict(deepcopy(Individual.mol))
    _ = [setattr(Individual, name, value) for name, value in zip(models, predictions)]


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