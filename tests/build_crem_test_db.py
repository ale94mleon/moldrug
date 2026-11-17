#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Build the small CReM fragment database used by the moldrug test suite.

The tests used to download ``replacements02_sc2.db.gz`` from qsar4u.com on
every run. That file is gone (the server now answers 403) and the precompiled
databases that replaced it, published on Zenodo (doi:10.5281/zenodo.16909329),
are 159 MB or more -- far too large to pull on every test run.

Instead the tests build their own tiny database from
``data/crem_test_mols.smi``. The database is *not* committed: CReM warns that a
change in how RDKit writes SMILES can invalidate a prebuilt fragment database,
so building it on the fly keeps it consistent with whatever RDKit version is
installed.

It is only meant to make the CReM operators return *something* for the test
seed molecules; it is far too small for any real work. For real problems use
one of the Zenodo databases, see ``docs/source/notes.md``.

Usage
-----
Build the database (done automatically by the test suite)::

    python build_crem_test_db.py --output crem.db

Regenerate the SMILES input from the NCI subset that ships with RDKit::

    python build_crem_test_db.py --regenerate-smiles
"""
import argparse
import os
import subprocess  # nosec
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SMILES_FILE = os.path.join(HERE, 'data', 'crem_test_mols.smi')

#: Name of the fragment set inside the database.
SET_NAME = 'moldrug_test'
#: Only the radii the tests actually ask for; more radii means a bigger database.
RADII = ('1', '2', '3')
#: Small cores keep the database small while still covering the test seeds.
MAX_HEAVY_ATOMS = '8'

# Parameters used to select the molecules in crem_test_mols.smi
NUM_MOLECULES = 1500
MIN_HEAVY_ATOMS = 8
MAX_MOL_HEAVY_ATOMS = 35
ALLOWED_ELEMENTS = frozenset(('C', 'N', 'O', 'S', 'F', 'Cl', 'Br'))


def build(output_db: str, ncpu: int = 1) -> str:
    """Build the test fragment database.

    Parameters
    ----------
    output_db : str
        Path of the SQLite database to create. It is overwritten if it exists.
    ncpu : int, optional
        Number of processes handed to ``cremdb_create``, by default 1

    Returns
    -------
    str
        ``output_db``, for convenience.

    Raises
    ------
    RuntimeError
        If the database could not be built.
    """
    if os.path.exists(output_db):
        os.remove(output_db)

    cmd = [
        sys.executable, '-m', 'crem.scripts.cremdb_create',
        '--input', SMILES_FILE,
        '--output', output_db,
        '--set-name', SET_NAME,
        '--radii', *RADII,
        '--max-heavy-atoms', MAX_HEAVY_ATOMS,
        '--ncpu', str(ncpu),
    ]
    process = subprocess.run(cmd, capture_output=True, text=True)  # nosec
    if process.returncode or not os.path.isfile(output_db):
        raise RuntimeError(
            f"Could not build the test CReM database with:\n{' '.join(cmd)}\n"
            f"stdout:\n{process.stdout}\nstderr:\n{process.stderr}")
    return output_db


def regenerate_smiles() -> str:
    """Rewrite :data:`SMILES_FILE` from the NCI subset distributed with RDKit.

    Only needed if the molecule selection has to change; the result is meant to
    be committed.

    Returns
    -------
    str
        Path of the written SMILES file.
    """
    from rdkit import Chem, RDConfig, RDLogger
    RDLogger.DisableLog('rdApp.*')

    source = os.path.join(RDConfig.RDDataDir, 'NCI', 'first_5K.smi')
    written = 0
    with open(source) as handler, open(SMILES_FILE, 'w') as out:
        for line in handler:
            if written >= NUM_MOLECULES:
                break
            smiles = line.split()[0]
            # Keep single components only, CReM fragments connected molecules
            if '.' in smiles:
                continue
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                continue
            if not MIN_HEAVY_ATOMS <= mol.GetNumHeavyAtoms() <= MAX_MOL_HEAVY_ATOMS:
                continue
            if any(atom.GetSymbol() not in ALLOWED_ELEMENTS for atom in mol.GetAtoms()):
                continue
            out.write(f"{Chem.MolToSmiles(mol)}\tm{written}\n")
            written += 1
    print(f"Wrote {written} molecules to {SMILES_FILE}")
    return SMILES_FILE


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-o', '--output', default='crem.db',
                        help='path of the database to build (default: %(default)s)')
    parser.add_argument('-c', '--ncpu', type=int, default=1,
                        help='number of processes used to build it (default: %(default)s)')
    parser.add_argument('--regenerate-smiles', action='store_true',
                        help='rewrite the committed SMILES file and exit')
    args = parser.parse_args()

    if args.regenerate_smiles:
        regenerate_smiles()
        return
    print(f"Built {build(args.output, ncpu=args.ncpu)}")


if __name__ == '__main__':
    main()
