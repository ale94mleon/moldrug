from rdkit import Chem
from featurize import Predictors, Featurizer

if __name__ == '__main__':
    mol = Chem.MolFromSmiles("c1ccccc1")
    adme_models = {
        "clearance":"clearance.jlib",
        "hppb":"hppb.jlib",
    }
    predictor = Predictors(featurizer=Featurizer(), models=adme_models.values())
    print(predictor.predict(mol))