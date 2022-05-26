#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import glob as glob
import multiprocessing as mp
import os, tempfile, tqdm
import numpy as np
import utility

#========Deffinitions of class and methosd to get the output of Vina===========
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
            
        
#===========================================================================



#==============================================================================

#=================Deffinition of the run=======================================
def VinaCost(idx, pdbqt, receptor_path, boxcenter, boxsize, exhaustiveness = 8, vina_cpus = 1,  num_modes = 1, wd = '.vina'):
    cmd = f"vina --receptor {receptor_path} --ligand {os.path.join(wd, f'{idx}.pdbqt')} "\
        f"--center_x {boxcenter[0]} --center_y {boxcenter[1]} --center_z {boxcenter[2]} "\
        f"--size_x {boxsize[0]} --size_y {boxsize[1]} --size_z {boxsize[2]} "\
        f"--out {os.path.join(wd, f'{idx}_out.pdbqt')} --cpu {vina_cpus} --exhaustiveness {exhaustiveness} --num_modes {num_modes}"
    print
    # Creating the ligand pdbqt
    with open(os.path.join(wd, f'{idx}.pdbqt'), 'w') as l:
        l.write(pdbqt)
    utility.run(cmd)

    # Getting the information
    best_energy = VINA_OUT(os.path.join(wd, f'{idx}_out.pdbqt')).BestEnergy()
    # Changing the xyz conformation by the conformation of the binding pose
    pdbqt = best_energy.chunk
    # Getting the Scoring function of Vina
    cost = best_energy.freeEnergy
    return pdbqt, cost
def VinaCostStar(args):
    return VinaCost(*args)
#==============================================================================



if __name__ == '__main__':
    pass
    
    vina = VINA_OUT(".vina/0_out.pdbqt")
    print(vina.BestEnergy().freeEnergy)
    #vina.PosNegConf(atom_1 = 17 , atom_2 = 1)
    #runvina("/home/ale/MY_PYTHON_PACKEGES/MDynamic/examples/Vina_Docking/ref_info/boxes.box")

