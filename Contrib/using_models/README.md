# ADME models

The fitness module: `fitness_plus_models.py` can use the models `hppb.jlib` and `clearence.jlib` by default and any other `jlib` model must be specified through the keywords: `models` and `desirability`. The models will be evaluated as follows:

```python
# The classes Featurizer and Predictors are already inside fitness_plus_admes.py
predictor = Predictors(featurizer=Featurizer(), models=models.values())
predictions = predictor.predict(Individual.mol)
```

The scripts `train_models.py`, `featurize.py` are used for the construction of the models. `use_models.py` is a simple example of how to use the models. These three scripts were kindly provided by [Miha Skalic](https://github.com/miha-skalic). The class `Featurizer` and `Predictors` from `featurize.py` are integrated into `fitness_plus_admes.py`.

## Getting the models

Install dependencies:

```bash
pip install sklearn
```

The other dependencies are already inside **MolDrug**. You must have in the same folder the scripts: `train_models.py`, `featurize.py` and then just:

```bash
python train_models.py
```

You should get the files: `clearance.jlib` and `hppb.jlib`.

## Future models

Here is a detailed list of public databases to construct ADME models: [adme](https://tdcommons.ai/single_pred_tasks/adme/).

## EGFR data set

This data was retrieved using this [EGFR-kinase tutorial](https://projects.volkamerlab.org/teachopencadd/talktorials/T001_query_chembl.html#Get-target-data-(EGFR-kinase)). The mentioned tutorial should encourage the user to retrieve data for proteins for which the 3D structure is not available and use ML/AI methods like the ones presented in this example and design molecules with MolDrug using those models as fitness functions.

## Run MolDrug from the command line

```bash
moldrug config.yml --fitness /the/path/for/fitness_plus_admes.py
```
