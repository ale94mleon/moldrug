#!/usr/bin/env python
"""Code borrowed from Pat Walters
https://github.com/PatWalters/fragment_expansion/blob/master/rdkit_eval/rd_gen_restricted_confs.py
Which was borrowed from Joshua Meyers
https://github.com/JoshuaMeyers/Snippets/blob/master/200405_constrained_conformers.ipynb
and that code was adapted from Tim Dudgeon
https://github.com/InformaticsMatters/pipelines/blob/master/src/python/pipelines/rdkit/constrained_conf_gen.py
All I've done is change the commandline wrapper and modify how to remove conformers that clash with the protein"""

from typing import List, Optional

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
from tqdm import tqdm
import argparse
from copy import deepcopy
import Bio.PDB as PDB

def duplicate_conformers(m: Chem.rdchem.Mol, new_conf_idx: int, rms_limit: float = 0.5) -> bool:
    rmslist = []
    for i in range(m.GetNumConformers()):
        if i == new_conf_idx:
            continue
        rms = AllChem.GetConformerRMS(m, new_conf_idx, i, prealigned=True)
        rmslist.append(rms)
    return any(i < rms_limit for i in rmslist)


def get_mcs(mol_one: Chem.rdchem.Mol, mol_two: Chem.rdchem.Mol) -> str:
    """Code to find the maximum common substructure between two molecules."""
    return Chem.MolToSmiles(
        Chem.MolFromSmarts(
            rdFMCS.FindMCS([mol_one, mol_two], completeRingsOnly=True, matchValences=True).smartsString
        )
    )


def generate_conformers(mol: Chem.rdchem.Mol,
                        ref_mol: Chem.rdchem.Mol,
                        num_conf: int,
                        ref_smi: str = None,
                        minimum_conf_rms: Optional[float] = None,
                        ) -> List[Chem.rdchem.Mol]:
    # if SMILES to be fixed are not given, assume to the MCS
    if not ref_smi:
        ref_smi = get_mcs(mol, ref_mol)

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
        AllChem.ConstrainedEmbed(temp_mol, core1, randomseed=i)
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


class ProteinLigandClashFilter:
    def __init__(self, protein_pdbpath: str, distance: float = 1.5):
        parser = PDB.PDBParser(QUIET=True, PERMISSIVE=True)
        s = parser.get_structure('protein', protein_pdbpath)
        self.kd = PDB.NeighborSearch(list(s.get_atoms()))
        self.radius = distance

    def __call__(self, conf: Chem.rdchem.Conformer) -> bool:
        for coord in conf.GetPositions():
            res = self.kd.search(coord, radius=self.radius)
            if len(res):
                return True
        return False


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument(
        '--pdb',
        help = 'Protein pdb file',
        dest = 'pdb',
        type = str,
    )
    parser.add_argument(
        '--smi',
        help='Input SMILES file name',
        dest = 'smi',
        type = str,
    )
    parser.add_argument(
        '--fix',
        help = 'File with fixed piece of the molecule',
        dest = 'fix',
        type = str,
    )
    parser.add_argument(
        '--out',
        help = 'Output file name',
        dest = 'out',
        type = str,
    )
    parser.add_argument(
        '--max',
        help = 'Maximum number of conformers to generate; default %(default)s',
        dest = 'max',
        default = 25,
        type = int,
    )
    parser.add_argument(
        '--rms',
        help = 'RMS cutoff; default %(default)s',
        dest = 'rms',
        default = 0.01,
        type = float,
    )
    parser.add_argument(
        '--bump',
        help = 'Bump cutoff default %(default)s',
        dest = 'bump',
        default = 1.5,
        type = float,
    )
    args = parser.parse_args()
    ref = Chem.MolFromMolFile(args.fix)
    suppl = Chem.SmilesMolSupplier(args.smi, titleLine=False)
    writer = Chem.SDWriter(args.out)

    clash_filter = ProteinLigandClashFilter(args.pdb, distance=args.bump)

    for mol in tqdm(suppl):
        # generate conformers
        out_mol = generate_conformers(mol, ref,
                                      args.max,
                                      ref_smi=Chem.MolToSmiles(ref),
                                      minimum_conf_rms=args.rms)

        # remove conformers that clash with the protein
        clashIds = [conf.GetId() for conf in out_mol.GetConformers() if clash_filter(conf)]
        [out_mol.RemoveConformer(clashId) for clashId in clashIds]

        # write out the surviving conformers
        for conf in out_mol.GetConformers():
            writer.write(out_mol, confId=conf.GetId())


if __name__ == "__main__":
    main()