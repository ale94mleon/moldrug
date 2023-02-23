from moldrug import utils
from fitness_plus_models import Cost
from rdkit import Chem
import tempfile, os
from moldrug.data import ligands, boxes, receptor_pdbqt
tmp_path = tempfile.TemporaryDirectory()
ligand_mol = Chem.MolFromSmiles(ligands.r_x0161)
I = utils.Individual(ligand_mol)
receptor_path = os.path.join(tmp_path.name,'receptor.pdbqt')
with open(receptor_path, 'w') as r: r.write(receptor_pdbqt.r_x0161)
box = boxes.r_x0161['A']

# Using the default desirability
NewI = Cost(
    Individual = I,wd = tmp_path.name,
    receptor_pdbqt_path = receptor_path,boxcenter = box['boxcenter'],
    boxsize = box['boxsize'],exhaustiveness = 4,ncores = 4,
    models = {
        'egfr': 'egfr.jlib',
        'hppb':  'hppb.jlib',
        'hppb_copy':  'hppb.jlib',
        'clearance': 'clearance.jlib',
        'clearance_copy': 'clearance.jlib',

    }, 
    desirability = {
        'egfr': {
            'w': 1,
            'LargerTheBest': {
                'LowerLimit': 4,
                'Target':10,
                'r': 1
            }
        },
        'hppb': {
            'w': 1,
            'LargerTheBest': {
                'LowerLimit': 25,
                'Target':75,
                'r': 1
            }
        },
        'hppb_copy': {
            'w': 1,
            'LargerTheBest': {
                'LowerLimit': 25,
                'Target':75,
                'r': 1
            }
        },
        'clearance': {
            'w': 1,
            'SmallerTheBest': {
                'Target': 20,
                'UpperLimit': 125,
                'r': 1
            }
        },
        'clearance_copy': {
            'w': 1,
            'SmallerTheBest': {
                'Target': 20,
                'UpperLimit': 125,
                'r': 1
            }
        },
    }
    
    )
print(
    f"cost: {NewI.cost}",
    f"vina_score: {NewI.vina_score}",
    f"egfr: {NewI.egfr}",
    f"clearance: {NewI.clearance}",
    f"hppb: {NewI.hppb}",
    f"clearance_copy: {NewI.clearance_copy}",
    f"hppb_copy: {NewI.hppb_copy}",)
