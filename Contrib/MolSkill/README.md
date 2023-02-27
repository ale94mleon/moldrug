# Scoring function with MolSkill

Recently, a promised "drug-likness" estimator was created: [MolSkill](https://github.com/microsoft/molskill). The fitness function presented in this directory, is a version of the standard `moldrug.fitness.cost`, only with QED substituted by MolSkill. Other cost functions (using models, multi-receptors, etc..) could easily be obtained in a similar way.

To use it you just have to install MolSkill (check its [GitHub](https://github.com/microsoft/molskill))

And, as usual:

```bash
moldrug config.yml --fitness /the/path/for/fitness_plus_molskill.py
```

This is the desirability used by default:

```python
desirability = {
    'molskill_score': {
        'w': 1,
        'SmallerTheBest': {
            'Target': -15,
            'UpperLimit': 0,
            'r': 1
        }
    },
    'sa_score': {
        'w': 1,
        'SmallerTheBest': {
            'Target': 3,
            'UpperLimit': 7,
            'r': 1
        }
    },
    'vina_score': {
        'w': 1,
        'SmallerTheBest': {
            'Target': -12,
            'UpperLimit': -6,
            'r': 1
        }
    }
}
```

Of course, you could change it for your specific application.