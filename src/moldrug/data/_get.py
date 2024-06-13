import os

import yaml

from moldrug import home

_AVAILABLE_DATA = ['6lu7', 'x0161']


class DataNotFound(Exception):
    pass


def get_data(name: str) -> dict:
    """Retrieve moldrug data

    Parameters
    ----------
    name : str
        Name of the data directory of moldrug

    Returns
    -------
    dict
        A dictionary with keys:

        * box: boxcenter, boxsize
        * ligand_3D: absolute path
        * smiles: the ligans's SMILES
        * protein: pdb, pdbqt

    Raises
    ------
    DataNotFound
        moldrug does not contains the data.
    """

    root_data = home.home(dataDir=name)

    if not os.path.exists(root_data) or not name:
        raise DataNotFound(f"moldrug data does not have '{name}' data set. "
                           f"Choose from: {_AVAILABLE_DATA}")

    with open(os.path.join(root_data, 'box.yml'), 'r') as f:
        box = yaml.safe_load(f)

    with open(os.path.join(root_data, 'ligand.smi'), 'r') as f:
        smiles = f.readline()

    data = {
        'box': box,
        'ligand_3D': os.path.join(root_data, 'ligand_3D.mol'),
        'smiles': smiles,
        'protein': {
            'pdb': os.path.join(root_data, 'protein.pdb'),
            'pdbqt': os.path.join(root_data, 'protein.pdbqt'),
        }
    }
    return data
