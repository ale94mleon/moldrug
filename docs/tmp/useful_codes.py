from rdkit import Chem
import random
import numpy as np
from rdkit.Chem import AllChem, rdmolops, DataStructs, Lipinski, Descriptors
from datetime import datetime
from sklearn.ensemble import RandomForestRegressor
import warnings
import os
from moldrug import utils
# utils

def fragments(mol:Chem.rdchem.Mol):
    """Create the fragments of the molecule based on the single bonds

    Args:
        mol (Chem.rdchem.Mol): An RDKit molecule.

    Returns:
        Chem.rdmolops.GetMolFrags: The fragments. you could cast the result using list and get all the fragments.
    """
    break_point = int(random.choice(np.where(np.array([b.GetBondType() for b in mol.GetBonds()]) == Chem.rdchem.BondType.SINGLE)[0]))
    # Chem.FragmentOnBonds(mol, break_point) # Could be used to increase randomness give more possible fragments and select two of them
    with Chem.RWMol(mol) as rwmol:
        b = rwmol.GetBondWithIdx(break_point)
        rwmol.RemoveBond(b.GetBeginAtomIdx(), b.GetEndAtomIdx())
    return rdmolops.GetMolFrags(rwmol, asMols = True)

def rdkit_numpy_convert(fp):
    """Convert a list of binary fingerprint to numpy arrays.

    Args:
        fp (_type_): list of binary fingerprints

    Returns:
        numpy.array: An array with dimension (number of elements in the list, number of bits on the fingerprint).
    """
    # fp - list of binary fingerprints
    output = []
    for f in fp:
        arr = np.zeros((1,))
        DataStructs.ConvertToNumpyArray(f, arr)
        output.append(arr)
    return np.asarray(output)

def get_top(ms:list, model):
    """Get the molecule with higher value predicted by the model.

    Args:
        ms (list): list of molecules
        model (sklearn model): A model for with the training set was a set of numpy array based on the AllChem.GetMorganFingerprintAsBitVect 

    Returns:
        _type_: The molecule with the higher predicted value.
    """
    fps1 = [AllChem.GetMorganFingerprintAsBitVect(m, 2) for m in ms]
    x1 = rdkit_numpy_convert(fps1)
    pred = model.predict(x1)
    i = np.argmax(pred)
    return ms[i], pred[i]

class ScoringPredictor:
    """This class is used to predict a scoring based on the information of the SMILES.
    """
    def __init__(self, smiles_scoring:dict, receptor:str, boxcenter:list, boxsize:list, exhaustiveness:int) -> None:
        # Esto puede servir para simplemente buscar en la base de datos por la moelcula, y si esta dar el valor exacto de Vina sin tenr que hacer el calculo
        # Y por supuesto para predecir
        # La prediccion en la parte del "crossover" y el si esta la moelcuela para el caclulo real
        self.lastupdate = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        self.smiles_scoring = smiles_scoring
        # This are control variables to save for future use
        self.receptor = receptor
        self.boxcenter = boxcenter
        self.boxsize = boxsize
        self.exhaustiveness = exhaustiveness
        self.bfp = rdkit_numpy_convert([AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(s), 2) for s in self.smiles_scoring])
        self.model = None

    
    def __call__(self, **keywords):
        # n_estimators=100, n_jobs=njobs*self.costfunc_keywords['vina_cpus'], random_state=42, oob_score=True
        self.model = RandomForestRegressor(**keywords)
        self.model.fit(self.bfp, list(self.smiles_scoring.values()))
    
    def update(self, new_smiles_scoring:dict, receptor:str, boxcenter:list, boxsize:list, exhaustiveness:int) -> None:
        # Control that the new introduced data belongs to the same model
        assert self.receptor == receptor, f"The original object was constructed with the receptor {self.receptor} and you introduced {receptor}. Consider to generate a new object."
        assert self.boxcenter == boxcenter, f"The original object was constructed with the boxcenter {self.boxcenter} and you introduced {boxcenter}. Consider to generate a new object."
        assert self.boxsize == boxsize, f"The original object was constructed with the boxsize {self.boxsize} and you introduced {boxsize}. Consider to generate a new object."
        assert self.exhaustiveness == exhaustiveness, f"The original object was constructed with the exhaustiveness {self.exhaustiveness} and you introduced {exhaustiveness}. Consider to generate a new object."
        
        # Update the date 
        self.lastupdate = datetime.now().strftime
        
        # Look for new structures
        smiles_scoring2use = dict()
        for sm_sc in new_smiles_scoring:
            if sm_sc not in self.smiles_scoring:
                smiles_scoring2use[sm_sc] = new_smiles_scoring[sm_sc]
        if smiles_scoring2use:
            self.smiles_scoring.update(smiles_scoring2use)

            self.bfp = np.vstack((self.bfp, rdkit_numpy_convert([AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(s), 2) for s in smiles_scoring2use])))
            self.model = None
            print(f"{len(smiles_scoring2use)} new structures will be incorporate to the model.")
        else:
            print("The introduced smiles are already in the data base. No need to update.")

    def predict(self, smiles:list):
        if self.model and smiles:
            if type(smiles) == str: smiles = [smiles]
            fps = [AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(s), 2) for s in smiles]
            return self.model.predict(rdkit_numpy_convert(fps))
        else:
            if not self.model:
                warnings.warn('There are not model on this object, please call it (__call__) to create it in order to get a prediction.  If you update the class, you must call the class again in order to create a new model based on the updated information, right now was set to None" Right now you just got None!')
                return None
            if not smiles:
                warnings.warn('You did not provide a smiles. You will just get None as prediction')
                return None
    
    def pickle(self,title, compress = False):
        cls = self.__class__
        result = cls.__new__(cls)
        result.__dict__.update(self.__dict__)
        if compress:
            compressed_pickle(title, result)
        else:
            full_pickle(title, result)


# Fitness
def __VinaCost(Individual, wd = '.vina_jobs', vina_executable = 'vina', receptor_path = None, boxcenter = None, boxsize =None, exhaustiveness = 8, ncores = 1,  num_modes = 1):
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
    Individual.cost = best_energy.freeEnergy
    return Individual

def __VinaCostLipinski(Individual, wd = '.vina_jobs', vina_executable = 'vina', receptor_path = None, boxcenter = None, boxsize =None, exhaustiveness = 8, ncores = 1,  num_modes = 1):
    if utils.lipinski_filter(Individual.mol):
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
        Individual.cost = best_energy.freeEnergy
    else:
        Individual.vina_score = np.inf
        Individual.cost = np.inf
    return Individual

def __QedSasVinaCost(Individual, wd = '.vina_jobs', vina_executable = 'vina', receptor_path = None, boxcenter = None, boxsize =None, exhaustiveness = 8, ncores = 1,  num_modes = 1):
    Individual.qed = QED.weights_mean(Individual.mol)
    
    # Getting synthetic accessibility score
    Individual.sa_score = sascorer.calculateScore(Individual.mol)
    
    if Individual.qed >=0.75 and Individual.sa_score <= 4:

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
        Individual.cost = best_energy.freeEnergy
    else:
        # Getting the Scoring function of Vina
        Individual.vina_score =np.inf
        Individual.cost = np.inf
    return Individual


def __CostSimilarity(Individual, ref_smiles, wd = '.vina_jobs', vina_executable = 'vina',  receptor_path = None, boxcenter = None, boxsize =None, exhaustiveness = 8, ncores = 1,  num_modes = 1):
    # multicriteria optimization,Optimization of Several Response Variables
    # Getting estimate of drug-likness
    Individual.qed = QED.weights_mean(Individual.mol)
    
    # Getting synthetic accessibility score
    Individual.sa_score = sascorer.calculateScore(Individual.mol)
    
    # Getting similarity

    Individual.similarity_to_ref = DataStructs.FingerprintSimilarity(
        AllChem.GetMorganFingerprintAsBitVect(Individual.mol, 2),
        AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(ref_smiles), 2)
    )


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
    # Similarity range from 0 to 1
    w_similarity = 1
    d_similarity = utils.LargerTheBest(Individual.qed, LowerLimit=0.5, Target=0.9, r = 1)
    # Vina. In this case is difficult to create a desirability because the range, we will use as target -12 upperlimit -6.
    w_vina_score = 1
    d_vina_score = utils.SmallerTheBest(Individual.vina_score, Target = -12, UpperLimit=-6, r = 1)
    # Average 
    #D = (w_qed*d_qed + w_sa_score*d_sa_score + w_vina_score*d_vina_score) / (w_qed + w_sa_score + w_vina_score)
    # Geometric mean
    D = (d_qed**w_qed * d_sa_score**w_sa_score * d_similarity**w_similarity * d_vina_score**w_vina_score)**(1/(w_qed + w_sa_score + w_similarity + w_vina_score))
    # And because we are minimizing we have to return 
    Individual.cost = 1 - D
    return Individual