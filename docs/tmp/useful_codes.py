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
import inspect
print(inspect.getargspec(Cost).args)

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


### Para el GA
## Problem!!
# Another way to overcome the problem on generating similar molecules is to instead of calculate the similarity respect to the whole molecule give some part of the reference structure, give the pharmacophore in the initialization of GA
# and in the cost function, first find the closer fragments and from there calculate the similarity. because always the generated molecules will be different
# The mutation that I am doing right now is more a crossover but with the CReM library instead with the population itself.
# Code a more conservative mutation (close to the one in the paper of the SMILES) to perform small mutations. (now is user controlled)
# We could also, instead a crossover two mutations with different levels. Because our 'mutate' is indeed some kind of crossover but with the CReM data base. So, what we could do is implement the mutation sof the paper that is at small scale (local optimization): the possibilities are point mutations (atom-by-atom), deletion, add new atoms
# hard_mutate and soft_mutate, our genetic operations
# For the best ligand we could predict the metabolites with sygma (apart of the module)
# Think about how to handle possibles Vina Crash (problems with the type of atoms). i think that is some problems happens with vina the cost function will be just np.inf
# Sometimes the pdbqt structure is not generated (probably problems in the conversion. In this cases the whole simulation crash ans should not be like this, this should rise some warning and continue discarding this structure)
# Apply filter for chemical elements to avoid crash on the vina function, in case that Vina crash we could use the predicted model
# Till now Vina never fails but could happen.
# Add more cpus for the generation of the conformers with RDKit
# The size of the ligands increase with the number of generations (if crossover is used even more)
# How to implement the rationality of where to grow, not just randomness. That could be the "crossover" operation, in fact the grow in a specific direction, based in (for example) the interaction network or the clashing avoid.
# I have to create a filter of atoms in order that vina doesn't fail because B and atoms like that Vina is not able to handle.





            # Creating the first model
            # Could be more pretty like creating a method that is make the model, y lo que se hace es que se gurda el modelo
            # Aqui se hace por primera vez pero para no repetir tanto codigo solo se llama a update model
            
            # if predictor_model:
            #     self.Predictor = predictor_model
            #     print('\nUpdating the provided model:\n')
            #     self.Predictor.update(
            #         new_smiles_scoring = dict(((individual.smiles,individual.vina_score) if individual.vina_score != np.inf else (individual.smiles,9999) for individual in self.SawIndividuals)),
            #         receptor = os.path.basename(self.costfunc_kwargs['receptor_path']).split('.')[0],
            #         boxcenter = self.costfunc_kwargs['boxcenter'],
            #         boxsize = self.costfunc_kwargs['boxsize'],
            #         exhaustiveness = self.costfunc_kwargs['exhaustiveness'],
            #     )
            #     self.Predictor(n_estimators=100, n_jobs=njobs*self.costfunc_ncores, random_state=42, oob_score=True)
            #     print('Done!')
            # else:
            #     print('\nCreating the first predicted model...')
            #     self.Predictor = vina.ScoringPredictor(
            #         smiles_scoring = dict(((individual.smiles,individual.vina_score) if individual.vina_score != np.inf else (individual.smiles,9999) for individual in self.SawIndividuals)),
            #         receptor = os.path.basename(self.costfunc_kwargs['receptor_path']).split('.')[0],
            #         boxcenter = self.costfunc_kwargs['boxcenter'],
            #         boxsize = self.costfunc_kwargs['boxsize'],
            #         exhaustiveness = self.costfunc_kwargs['exhaustiveness'],
            #     )
            #     self.Predictor(n_estimators=100, n_jobs=njobs*self.costfunc_ncores, random_state=42, oob_score=True)
            #     print('Done!')

            # print(f'The model presents a oob_score = {self.Predictor.model.oob_score_}\n')

            # new_smiles_cost = dict()
                    # # New variables
                    # if individual.cost == np.inf:
                    #     new_smiles_cost[individual.smiles] = 9999
                    # else:
                    #     new_smiles_cost[individual.smiles] = individual.cost

            # # Update the model
            # print(f'Updating the current model with the information of generation {self.NumGen}...')

            # self.Predictor.update(
            #     new_smiles_scoring = new_smiles_cost.copy(),
            #     receptor = os.path.basename(self.costfunc_kwargs['receptor_path']).split('.')[0],
            #     boxcenter = self.costfunc_kwargs['boxcenter'],
            #     boxsize = self.costfunc_kwargs['boxsize'],
            #     exhaustiveness = self.costfunc_kwargs['exhaustiveness'],
            # )
            # self.Predictor(n_estimators=100, n_jobs=njobs*self.costfunc_ncores, random_state=42, oob_score=True)          
            # print('Done!')
            # print(f'The updated model presents a oob_score = {self.Predictor.model.oob_score_}')

    # def crossover(self, individual1, individual2, ncores = 1, probability = 0, MaxRatioOfIncreaseInWt = 0.25):
    #     # here I have to select some randomness to perform or not the real crossover because I think that we could get far from the solution. It is just a guess.
    #     # How do I control the size of the new offspring? 
    #     # Performing a fragmentation in such a way that the offspring is the same in size
    #     # Here is where more additional information could be used. In order to orient the design of the new offspring. 
    #     # Then, I should control how perform the mutation  in such a way that we could keep or at least evaluate the offspring generated for crossover
    #     if random.random() < probability: # 50% of return the same individuals
    #         fragments1 = fragments(individual1.mol)
    #         fragments2 = fragments(individual2.mol)
    #         all_fragments = list(fragments1) + list(fragments2)
            
    #         # Initialize offspring smiles; cost
    #         offsprings = [
    #                 [None, np.inf],
    #                 [None, np.inf],
    #         ]
    #         for combination in itertools.combinations(all_fragments, 2):
                
    #             # Combine the molecules
    #             try:
    #                 possible_fragments_smiles = list(link_mols(*combination, db_name=self.crem_db_path, radius = 3, min_atoms=1, max_atoms=6, dist = 2, return_mol=False, ncores=ncores))                
    #             except:
    #                 # This is for debugging
    #                 sm1, sm2 = [Chem.MolToSmiles(c) for c in combination]
    #                 raise RuntimeError(f'These are the problematic SMILES: {sm1}, {sm2}')
                
    #             # Perform a filter based on weight. This control the size of the fragments. For now I will test 25 %. Think in the future work with the mols instead of smiles, I have to convert to mols too many times in this section of the code
    #             avg_wt  = 0.5*(Descriptors.ExactMolWt(individual1.mol) + Descriptors.ExactMolWt(individual1.mol))
    #             threshold_wt = (MaxRatioOfIncreaseInWt + 1) * avg_wt
    #             print(f'We had {len(possible_fragments_smiles)} possible fragments')
    #             possible_fragments_smiles = list(filter(lambda x: Descriptors.ExactMolWt(Chem.MolFromSmiles(x)) < threshold_wt, possible_fragments_smiles))
    #             print(f'After the weight filter we have {len(possible_fragments_smiles)} possible fragments')

    #             # In case that it was not possible to link the fragments
    #             if not possible_fragments_smiles:continue

    #             # Bias the searching to similar molecules
    #             if self.get_similar:
    #                 possible_fragments_mols = get_similar_mols(mols = [Chem.MolFromSmiles(smiles) for smiles in possible_fragments_smiles], ref_mol=self.InitIndividual.mol, pick=self.popsize, beta=0.01)
    #                 possible_fragments_smiles = [Chem.MolToSmiles(mol) for mol in possible_fragments_mols]
                
    #             # Here comes the prediction with the model, and get the top two
    #             temp_offsprings = list(zip(possible_fragments_smiles, self.Predictor.predict(possible_fragments_smiles).tolist()))
                
    #             # Merge, Sort and Select
    #             offsprings = sorted(offsprings + temp_offsprings, key = lambda x:x[1])[:2]
    #         # Here I should check that exist offsprings (there not None values as smiles). For now I will assume that we always get at least two. See on the future
    #         return Individual(smiles = offsprings[0][0]), Individual(smiles = offsprings[1][0])
    #     else:
    #         return individual1, individual2  

    def mutate(self, individual):
        # See the option max_replacment
        # Or select the mutant based on some criterion
        # try:
            # Here i will pick the molecules based on the model.
        # El problema de seleccionar asi los compuestos es que siempre seleccionamos los mismos. Siempre se esta entrando la misma estructura y terminamos con una pobalcion redundante
        # Esto tengo que pensarlo mejor
        # new_mols = list(mutate_mol(Chem.AddHs(individual.mol), self.crem_db_path, radius=3, min_size=1, max_size=8,min_inc=-3, max_inc=3, return_mol=True, ncores = ncores))
        # new_mols = [Chem.RemoveHs(i[1]) for i in new_mols]
        # best_mol, score = get_top(new_mols + [individual.mol], self.model)
        # smiles = Chem.MolToSmiles(best_mol)
        # mol = best_mol
        # print(score)
        # For now I am generating all the mutants and picking only one at random, this is very inefficient, should be better only generate one, but I am afraid that crem generate always the same or not generate any at all.
        # I think that what first crem does is randomly select on spot and find there all possible mutants. If this spot doesn't generate mutants, then you don't get nothing. But this is a supposition. 
if __name__ == '__main__':
    i1 = Individual(smiles = 'CCCO', cost = 10)
    i2 = Individual(smiles = 'CCC', cost = 2)
    i3 = Individual(smiles = 'CCCF', cost = 3)
    i4 = Individual(smiles = 'CCCF', cost = 69)
    i5 = Individual(smiles = 'CCCCCCF', cost = 69)
    # print(i1)
    # print(f"i1 == i2 :{i1 == i2}")
    # print(f"i1 != i2: {i1 != i2}")
    # print(f"i1 <= i2: {i1 <= i2}")
    # print(f"i1 >= i2: {i1 >= i2}")
    # print(f"i1 < i2: {i1 < i2}")
    # print(f"i1 > i2: {i1 > i2}")
    # print(f"i1 is i2: {i1 is i2}")
    # print(f"i1 + i2: {i2 + i1}")
    # print(f"i1 * i2: {i1 * i1}")
    print(f"i1 / i2: { i5//6}")
    print(-i1)
    # print(f"np.mean: {np.mean([i1, i2])}")
    # print(f"bool: {bool(i2)}")
    # print(f"min {min([i1,i2,i3])}")
    # print(f"i4 in [] { 'CCCF' in  [i1,i2,i3]}")
    
    # for i in [i4, i5, Individual('CClCC')]:
    #     if i not in [i1,i2,i3]:
    #         print(i)
    array1=[i1,i2,i3,i4]
    array1[0].pdbqt
    print(make_sdf(array1))
    array2=[i1,i2,i3,i4]
    print(sorted(array1), 55555555555)
    # Probabilities Selections
    costs = np.array(array1)
    probs = np.exp((2*costs).astype('float64'))# / np.sum(np.exp(36*costs))
    print(costs[0])