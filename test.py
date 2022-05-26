#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import ga
A = ga.Individual(smiles = 'O=C(C)Oc1ccccc1C(=O)O', cost = 56)
out = ga.GA(A, ga.VinaCost, 20, 30, crem_db_path = '/home/ale/GITLAB/bi_crem_database/replacements02_sc2.5.db')