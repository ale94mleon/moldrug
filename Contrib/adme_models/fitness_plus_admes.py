#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from copy import deepcopy
from moldrug import utils, constraintconf
from rdkit import Chem
import os
import numpy as np
from typing import Dict, List
import warnings
from meeko import MoleculePreparation

import joblib
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')


class Featurizer:
    def __init__(self, fpb=2048, calc_descriptors=True):
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
        return [model.predict(X[np.newaxis,:])[0] for model in self.models]




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
    adme_models:Dict = None,
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
    adme_models : dict = None,
        The definition of the models,
        by default None which means that it will be used:    
        adme_models = {
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
        The date set used have the following limits hppb = [11.18; 99.95] (higher is better); clearance = [3.0; 150.0] (lower is better). But your models could change.
    Returns
    -------
    utils.Individual
        A new instance of the original Individual with the the new attributes:
        pdbqt, qed, vina_score, sa_score and cost.
        cost attribute will be a number between 0 and 1, been 0 the optimal value.
    """
    if not adme_models:
        adme_models = {
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
    # Check that everything is ok with naming in adme_models and desirability
    diff = list(set(desirability) - set(adme_models))
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
            raise Exception(f"You provided adme_models = {adme_models.keys()} and desirability = {desirability.keys()}. "\
                "However, desirability must have the same keywords as adme_models and optionally the keyword vina_score")
        elif diff[0] != 'vina_score':
            raise Exception(f"desirability has the keyword: '{diff[0]}' which is not defined in adme_models =  {adme_models.keys()} and is not vina_score")

    # Getting and setting properties on the Individual
    predictor = Predictors(featurizer=Featurizer(), models=adme_models.values())
    predictions = predictor.predict(Individual.mol)
    _ = [setattr(Individual, name, value) for name, value in zip(adme_models, predictions)]


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