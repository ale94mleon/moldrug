# ADME models

The fitness module: `fitness_plus_admes.py` is able to use the models `hppb.jlib` and `clearence.jlib` by default and any other `jlib` model must be specified through the keywords: `adme_models` and `desirability`. The models will be evaluated as follow:

```python
# The classes Featurizer and Predictors are already inside fitness_plus_admes.py
predictor = Predictors(featurizer=Featurizer(), models=adme_models.values())
predictions = predictor.predict(Individual.mol)
```

The scripts `train_models.py`, `featurize.py` are used for the construction sof the models. `use_models.py` is a simple example on how to use the models. These three scripts were kindly provided by [Miha Skalic](https://github.com/miha-skalic). The class `Featurizer` and `Predictors` from `featurize.py` are integrated in `fitness_plus_admes.py`.

## Getting the models

Install dependencies:

```bash
pip install sklearn
```

The other dependencies are already inside **MolDrug**. Yopu must have in the same folder the scripts: `train_models.py`, `featurize.py` and then just:

```bash
python train_models.py
```

You should get the files: `clearance.jlib` and `hppb.jlib`.

## Future models

Here  is a detail list of public data base to construct ADME models: [adme](https://tdcommons.ai/single_pred_tasks/adme/).
