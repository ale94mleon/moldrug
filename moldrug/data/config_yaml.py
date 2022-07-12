r_x0161 = """
type: GA
njobs: 3
seed_smiles: COC(=O)C=1C=CC(=CC1)S(=O)(=O)N
costfunc: Cost
costfunc_kwargs:
  vina_executable: vina
  receptor_path: x0161.pdbqt
  boxcenter:
    - 12.11
    - 1.84
    - 23.56
  boxsize:
    - 22.5
    - 22.5
    - 22.5
  exhaustiveness: 4
  ncores: 4
  num_modes: 1
crem_db_path: /home/ale/GITLAB/bi_crem_database/replacements02_sc2.5.db
maxiter: 10
popsize: 10
beta: 0.001
pc: 1
get_similar: False
mutate_crem_kwargs:
  radius: 3
  min_size: 1
  max_size: 8
  min_inc: -5
  max_inc: 3
  ncores: 12
"""