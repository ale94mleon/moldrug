#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
For information of MolDrug:
    Docs: https://moldrug.readthedocs.io/en/latest/
    Source Code: https://github.com/ale94mleon/moldrug
"""
import os
from moldrug._version import __version__

__author__ = "Alejandro Martínez León"
__email__ = "ale94mleon@gmail.com"

if "MOLDRUG_VERBOSE" in os.environ:
    if os.environ["MOLDRUG_VERBOSE"].lower() in [1, 'true']:
        verbose = True
    elif os.environ["MOLDRUG_VERBOSE"].lower() in [0, 'false']:
        verbose = False
    else:
        raise ValueError(f"MOLDRUG_VERBOSE = {os.environ['MOLDRUG_VERBOSE']} is invalid. Choose from: 1, true, false (case insensitive).")
else:
    verbose = False
