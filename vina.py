#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from lead import utils
import os, warnings
from datetime import datetime
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
import pickle as _pickle
from sklearn.ensemble import RandomForestRegressor

import os


# Alias for vina
vina_executable = os.path.abspath('vina')
#========Definitions of class and methods to get the output of Vina===========
class Atom:
    #https://userguide.mdanalysis.org/stable/formats/reference/pdbqt.html#writing-out
    def __init__(self, line):
        self.lineType = "ATOM"
        self.serial = int(line[6:11])
        self.name = line[12:16].strip()
        self.altLoc = line[16]
        self.resName = line[17:21].strip()
        self.chainID = line[21]
        self.resSeq = int(line[22:26])
        self.iCode = line[26]
        self.x = float(line[30:38])
        self.y = float(line[38:46])
        self.z = float(line[46:54])
        self.occupancy = line[54:60].strip()
        self.tempFactor = line[60:66].strip()
        self.partialChrg = line[66:76].strip()
        self.atomType = line[78:80].strip()
    def __getitem__(self, key):
        return self.__dict__[key]

class Hetatm:
    def __init__(self, line):
        self.lineType = "HETATM"
        self.serial = int(line[6:11])
        self.name = line[12:16].strip()
        self.altLoc = line[16]
        self.resName = line[17:21].strip()
        self.chainID = line[21]
        self.resSeq = int(line[22:26])
        self.iCode = line[26]
        self.x = float(line[30:38])
        self.y = float(line[38:46])
        self.z = float(line[46:54])
        self.occupancy = line[54:60].strip()
        self.tempFactor = line[60:66].strip()
        self.partialChrg = line[66:76].strip()
        self.atomType = line[78:80].strip()
    def __getitem__(self, key):
        return self.__dict__[key]

class Remark:
    def __init__(self, line):
        pass

class CHUNK_VINA_OUT:
    def __init__(self, chunk):
        self.chunk = chunk
        self.atoms = []
        self.run = None
        self.freeEnergy = None
        self.RMSD1 = None
        self.RMSD2 = None
        self.parse()

    def parse(self):
        for line in self.chunk:
            if line.startswith("MODEL"):
                self.run = int(line[5:])
            elif line.startswith("REMARK VINA RESULT:"):
                    (self.freeEnergy, self.RMSD1, self.RMSD2) = [float(number) for number in line.split(":")[-1].split()]
                    
            elif line.startswith("ATOM"):
                self.atoms.append(Atom(line))
            else:
                pass

    def get_atoms(self):
        """Return a list of all atoms.

        If to_dict is True, each atom is represented as a dictionary.
        Otherwise, a list of Atom objects is returned."""
        return [x.__dict__ for x in self.atoms]
        
    def write(self, name = None):
        if name:
            with open(name,"w") as f:
                f.writelines(self.chunk)
        else:
            with open(f"Run_{self.run}.pdbqt","w") as f:
                f.writelines(self.chunk)            

class VINA_OUT:
    """
    To acces the chunks you need to take into account that 
    """
    def __init__(self, file):
        self.file = file

        self.chunks = []
        self.parse()

    def parse(self):
        with open(self.file, "r") as input_file:
            lines = input_file.readlines()
        i = 0
        while i < len(lines):

            if lines[i].startswith("MODEL"):
                j = i
                tmp_chunk = []
                while (not lines[j].startswith("ENDMDL")) and (j < len(lines)):
                    tmp_chunk.append(lines[j])
                    j += 1
                    i += 1
                tmp_chunk.append("ENDMDL\n")
                self.chunks.append(CHUNK_VINA_OUT(tmp_chunk))

            i += 1
            
    def BestEnergy(self, write = False):
        min_chunk = min(self.chunks, key= lambda x: x.freeEnergy)
        if write: min_chunk.write("best_energy.pdbqt")
        return min_chunk


class VinaScoringPredictor:

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
        self.bfp = utils.rdkit_numpy_convert([AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(s), 2) for s in self.smiles_scoring])
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

            self.bfp = np.vstack((self.bfp, utils.rdkit_numpy_convert([AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(s), 2) for s in smiles_scoring2use])))
            self.model = None
            print(f"{len(smiles_scoring2use)} new structures will be incorporate to the model.")
        else:
            print("The introduced smiles are already in the data base. No need to update.")

    def predict(self, smiles:list):
        if self.model and smiles:
            if type(smiles) == str: smiles = [smiles]
            fps = [AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(s), 2) for s in smiles]
            return self.model.predict(utils.rdkit_numpy_convert(fps))
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
            utils.compressed_pickle(title, result)
        else:
            utils.full_pickle(title, result) 
#=================For compatibility=======================================
def VinaCost():
    return None
def VinaCostStar():
    return None
#==============================================================================

if __name__ == '__main__':
    pass
    # vina = VINA_OUT(".vina/0_out.pdbqt")
    # print(vina.BestEnergy().freeEnergy)
    #vina.PosNegConf(atom_1 = 17 , atom_2 = 1)
