#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
For the functions: duplicate_conformers, get_mcs, generate_conformers, constraintconf
and the class ProteinLigandClashFilter:
    Code borrowed from Pat Walters
    https://github.com/PatWalters/fragment_expansion/blob/master/rdkit_eval/rd_gen_restricted_confs.py
    Which was borrowed from Joshua Meyers
    https://github.com/JoshuaMeyers/Snippets/blob/master/200405_constrained_conformers.ipynb
    and that code was adapted from Tim Dudgeon
    https://github.com/InformaticsMatters/pipelines/blob/master/src/python/pipelines/rdkit/constrained_conf_gen.py
    All I've done is change the command line wrapper, modify how to remove conformers that clash with the protein, add documentation and
    handle some possible run time exceptions.
"""
from copy import deepcopy
from rdkit import Chem
from rdkit.Chem import AllChem, rdFMCS
import os
from typing import  Optional
import warnings
from tqdm import tqdm
import Bio.PDB as PDB
from moldrug.utils import compressed_pickle
def duplicate_conformers(m: Chem.rdchem.Mol, new_conf_idx: int, rms_limit: float = 0.5) -> bool:
    """
    Check if a conformer with index new_onf_idx it is duplicated based on rms_limit

    Parameters
    ----------
    m : Chem.rdchem.Mol
        An RDKit molecule with conformers.
    new_conf_idx : int
        Index of the conformer to check.
    rms_limit : float, optional
        Threshold of rms to consider duplicate structure, by default 0.5.

    Returns
    -------
    bool
        True if the conformer new_conf_idx id is duplicate. false otherwise.
    """
    rmslist = []
    for i in range(m.GetNumConformers()):
        if i == new_conf_idx:
            continue
        rms = AllChem.GetConformerRMS(m, new_conf_idx, i, prealigned=True)
        rmslist.append(rms)
    return any(i < rms_limit for i in rmslist)

def get_mcs(mol_one: Chem.rdchem.Mol, mol_two: Chem.rdchem.Mol) -> str:
    """
    Code to find the maximum common substructure between two molecules.

    Parameters
    ----------
    mol_one : Chem.rdchem.Mol
        The first molecule.
    mol_two : Chem.rdchem.Mol
        The second molecule.

    Returns
    -------
    str
        The SMILES string of the Maximum Common Substructure (MCS).
    """
    mcs_smarts = Chem.MolFromSmarts(rdFMCS.FindMCS([mol_one, mol_two], completeRingsOnly=True, matchValences=True).smartsString)
    mcs_smi = Chem.MolToSmiles(mcs_smarts)
    # Workaround in case of fails
    if not Chem.MolFromSmiles(mcs_smi):
        valid_atoms = mol_one.GetSubstructMatch(mcs_smarts)
        mcs_smi = Chem.MolFragmentToSmiles(mol_one, atomsToUse=valid_atoms)
    return mcs_smi

def gen_aligned_conf(mol: Chem.rdchem.Mol, ref_mol: Chem.rdchem.Mol, ref_smi:str):
    """Generate a conformation of mol aligned to ref_mol

    Parameters
    ----------
    mol : Chem.rdchem.Mol
        Molecule to generate the conformation aligned to ref_mol
    ref_mol : Chem.rdchem.Mol
        Reference molecule with a conformation
    ref_smi : str
        The SMILES string for which the aligned will take place.

    Returns
    -------
    _type_
        _description_
    """
    AllChem.EmbedMolecule(mol)
    # AllChem.MMFFOptimizeMolecule(mol)
    ref_smi_mol = Chem.MolFromSmiles(ref_smi)
    Chem.rdMolAlign.AlignMol(
        mol,
        ref_mol,
        atomMap = list(zip(mol.GetSubstructMatch(ref_smi_mol), ref_mol.GetSubstructMatch(ref_smi_mol)))
        )
    return mol

def generate_conformers(mol: Chem.rdchem.Mol,
                        ref_mol: Chem.rdchem.Mol,
                        num_conf: int,
                        ref_smi: str = None,
                        minimum_conf_rms: Optional[float] = None,
                        ) -> Chem.rdchem.Mol:
    """
    Generate constrained conformers

    Parameters
    ----------
    mol : Chem.rdchem.Mol
        The molecule for which the conformer will be generated.
    ref_mol : Chem.rdchem.Mol
        The reference structure from which the reference coordinates will be taken.
    num_conf : int
        Maximum number of conformer to generate.
    ref_smi : str, optional
        part of the molecule to keep fix. If ref_smi is not given, assume to the MCS between mol and ref_mol, by default None
    minimum_conf_rms : Optional[float], optional
        Threshold of rms to consider duplicate structure, by default None

    Returns
    -------
    Chem.rdchem.Mol
        mol with with generated conformers. In case some Exception ocurred, it returns mol without conformers and write in the
        current directory generate_conformers_error.log with the nature of the Exception.
    """
    # if SMILES to be fixed are not given, assume to the MCS
    if ref_smi:
        if not Chem.MolFromSmiles(ref_smi):
            raise ValueError("The provided ref_smi is not valid.")
    else:
        ref_smi = get_mcs(mol, ref_mol)
        if not Chem.MolFromSmiles(ref_smi):
            raise ValueError("generate_conformers fails generating ref_smi based on the MCS between mol and ref_mol")

    try:
        # Creating core of reference ligand #
        core_with_wildcards = AllChem.ReplaceSidechains(ref_mol, Chem.MolFromSmiles(ref_smi))
        core1 = AllChem.DeleteSubstructs(core_with_wildcards, Chem.MolFromSmiles('*'))
        core1.UpdatePropertyCache()

        # Add Hs so that conf gen is improved
        mol.RemoveAllConformers()
        outmol = deepcopy(mol)
        mol_wh = Chem.AddHs(mol)

        # Generate conformers with constrained embed
        dup_count = 0
        for i in range(num_conf):
            temp_mol = Chem.Mol(mol_wh)  # copy to avoid inplace changes
            try:
                AllChem.ConstrainedEmbed(temp_mol, core1, randomseed=i)
            except Exception:
                temp_mol = gen_aligned_conf(temp_mol, ref_mol, ref_smi)
            temp_mol = Chem.RemoveHs(temp_mol)
            conf_idx = outmol.AddConformer(temp_mol.GetConformer(0), assignId=True)
            if minimum_conf_rms is not None:
                if duplicate_conformers(outmol, conf_idx, rms_limit=minimum_conf_rms):
                    dup_count += 1
                    outmol.RemoveConformer(conf_idx)
        if dup_count:
            pass
            # print(f'removed {dup_count} duplicated conformations')
        return outmol
    except Exception as e:
        print(e)
        cwd = os.getcwd()
        warnings.warn(f"generate_conformers failed. Check the file {os.path.join(cwd, 'generate_conformers_error.pbz2')}")
        compressed_pickle('generate_conformers_error', e)
        mol.RemoveAllConformers()
        return mol

class ProteinLigandClashFilter:
    """
    Class used to eliminate clash between a ligand and a protein.
    """
    def __init__(self, protein_pdbpath: str, distance: float = 1.5):
        """
        This is the constructor of the class

        Parameters
        ----------
        protein_pdbpath : str
            Path top the protein pdb file
        distance : float, optional
            Threshold of distance to consider a clash, by default 1.5
        """
        parser = PDB.PDBParser(QUIET=True, PERMISSIVE=True)
        s = parser.get_structure('protein', protein_pdbpath)
        self.kd = PDB.NeighborSearch(list(s.get_atoms()))
        self.radius = distance

    def __call__(self, conf: Chem.rdchem.Conformer) -> bool:
        """
        Call deffinition

        Parameters
        ----------
        conf : Chem.rdchem.Conformer
            The conformer to be evaluated

        Returns
        -------
        bool
            True if there clash, False otherwise.
        """
        for coord in conf.GetPositions():
            res = self.kd.search(coord, radius=self.radius)
            if len(res):
                return True
        return False

def constraintconf(pdb:str, smi:str, fix:str, out:str, max_conf:int = 25, rms:float = 0.01, bump:float = 1.5):
    """
    It generates several conformations in the binding pocket with a specified constraint.

    Parameters
    ----------
    pdb : str
        Protein pdb file
    smi : str
        Input SMILES file name
    fix : str
        File with fixed piece of the molecule
    out : str
        Output file name
    max_conf : int, optional
        Maximum number of conformers to generate, by default 25
    rms : float, optional
        RMS cutoff, by default 0.01
    bump : float, optional
        Bump cutoff, by default 1.5
    """

    ref = Chem.MolFromMolFile(fix)
    suppl = Chem.SmilesMolSupplier(smi, titleLine=False)
    writer = Chem.SDWriter(out)

    clash_filter = ProteinLigandClashFilter(pdb, distance=bump)

    for mol in tqdm(suppl):
        # generate conformers
        out_mol = generate_conformers(mol, ref,
                                      max_conf,
                                      ref_smi=Chem.MolToSmiles(ref),
                                      minimum_conf_rms=rms)

        # remove conformers that clash with the protein
        clashIds = [conf.GetId() for conf in out_mol.GetConformers() if clash_filter(conf)]
        _ = [out_mol.RemoveConformer(clashId) for clashId in clashIds]

        # write out the surviving conformers
        for conf in out_mol.GetConformers():
            writer.write(out_mol, confId=conf.GetId())


if __name__ == '__main__':
    pass