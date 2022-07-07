import lead.data as data
from lead import utils, fitness
from lead.data import receptors, ligands, boxes
import tempfile, os, gzip, shutil, requests
from multiprocessing import cpu_count

# Creating a temporal directory
tmp_path = tempfile.TemporaryDirectory()
# Getting receptor a ligand information
ligand = ligands.r_x0161
receptor = receptors.r_x0161
receptor_file = os.path.join(tmp_path.name, 'receptor.pdbqt')
with open(receptor_file, 'w') as r:
    r.write(receptor)
box = boxes.r_x0161["A"]
# Getting the crem data base
url = "http://www.qsar4u.com/files/cremdb/replacements02_sc2.db.gz"
r = requests.get(url, allow_redirects=True)
open(os.path.join(tmp_path,'crem.db.gz'), 'wb').write(r.content)

crem_dbgz_path = os.path.join(tmp_path,'crem.db.gz')
crem_db_path = '/home/ale/GITLAB/lead/row_data/replacements02_sc2.db'

with gzip.open(crem_dbgz_path, 'rb') as f_in:
    with open(crem_db_path, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)




def test_single_receptor():

    maxiter = 1
    popsize = 2
    njobs = 3
    NumbCalls = 1

    """For a local optimization we could use
        min_size=1, max_size=1, min_inc=-1, max_inc=1
        Or use just grow instead
        min_size=0, max_size=0,
        And if we would like to integrate all this option
        min_size=0, max_size=1, min_inc=-1, max_inc=1
        this will add, delate heavy mutate heavy atoms.
        or change Hydrogens for heavy atoms

        with
        min_size=1, max_size=8, min_inc=-5, max_inc=3
        I explore efficiently the chemical space, not good for cost functions with similarity included. 
    """
    
    out = utils.GA(
        seed_smiles=ligand,
        maxiter=maxiter,
        popsize=popsize,
        crem_db_path = '/home/ale/GITLAB/bi_crem_database/replacements02_sc2.db',
        pc = 1,
        get_similar = False,
        mutate_crem_kwargs = {
            'radius':3,
            'min_size':1,
            'max_size':8,
            'min_inc':-5,
            'max_inc':3,
            'ncores':cpu_count(),
        },
        costfunc = fitness.Cost,#__CostSimilarity,# __VinaCostLipinski, Cost, __VinaCost, __QedSasVinaCost, CostMultiReceptors
        costfunc_kwargs = {
            'vina_executable': 'vina',
            'receptor_path': receptor_file,
            'boxcenter' : box['boxcenter'] ,
            'boxsize': box['boxsize'],
            'exhaustiveness': 8,
            'ncores': int(cpu_count() / njobs),
            'num_modes': 1,
            #'ref_smiles': init_smiles,
        },
        save_pop_every_gen = 20,
        pop_file_name = os.path.join(tmp_path.name, 'pop')
        )

    for i in range(NumbCalls):
        out(njobs = njobs)

    for o in out.pop:
        print(o.smiles, o.cost)
    out.pickle(os.path.join(tmp_path.name, f"result_NumGen_{out.NumGen}_PopSize_{popsize}"), compress=True)
    print(out.to_dataframe())

os.remove(crem_db_path)