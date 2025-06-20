#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# This module was taken from Meeko 0.6.1
# https://github.com/forlilab/meeko
# Only unnecessary functionalities were removed.

__meeko_original_version__ = "0.6.1"

from .preparation import MoleculePreparation
from .molecule_pdbqt import PDBQTMolecule
from .rdkit_mol_create import RDKitMolCreate
from .writer import PDBQTWriterLegacy

import logging
from rdkit import rdBase
rdkit_logger = logging.getLogger("rdkit")
rdkit_logger.handlers[0].setLevel("WARNING")
rdkit_logger.handlers[0].setFormatter(
    logging.Formatter('[RDKit] %(levelname)s:%(message)s'),
)
rdBase.LogToPythonLogger()

__all__ = ['MoleculePreparation',
           'PDBQTMolecule',
           'RDKitMolCreate',
           'PDBQTWriterLegacy',
           ]
