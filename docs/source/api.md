# Summary

## Input data

Currently, **moldrug** only accepts valid RDKit molecules and a valid pdbqt file that will be processed in case AutoDock-Vina is used.

## The idea

The idea is simple: we want the best drug for some applications. The problem is that it is not as simple as it sounds. The chemical space is enormous and the estimation of the fitness is also a term of concern.

Here we address this problem using a [Genetic Algorithm](https://moldrug.readthedocs.io/en/latest/source/modules/utils.html#moldrug.utils.GA). The inputs are:

- Biological target (pdbqt format).
- RDKit molecule of some known (or guessing) drug. (it could get easily get from the SMILES string)
- Definition of the binding pocket.
- Definition of the {py:class}`moldrug.utils.GA` running variables.

With the molecule, a random population of `popsize` individuals will be generated through the [CReM](https://github.com/DrrDom/crem) operation `mutate_mol`. Every individual will be evaluated through the selected fitness function. The best individual will be used for crossover and mutation, generating `pc*popsize` offspring. The offspring will be merged with the current population, sorted with respect to the fitness function and the first `popsize` individuals will be kept for the next generation. This cycle will run for `maxiter` generations.

In a general way **moldrug** will try to minimize:

```{math}
\mathrm{cost} = f(\vec{X})
```

Where {math}`\vec{X}` is some description of the chemical structure of the molecule which is mapped to the
{math}`\mathbb{R}` numbers through the function {math}`f`. This function will be called cost function.

Alternatively, there is also the class {py:class}`moldrug.utils.Local`
In this case, only an initial population is created and exported.

## Fitness functions

The default fitness functions could be accessed through the module {py:mod}`moldrug.fitness`.
There are currently four implemented:

- {py:func}`moldrug.fitness.Cost` (standard). It uses all the response variables.
- {py:func}`moldrug.fitness.CostOnlyVina`. Only use the information of the Vina score.
- {py:func}`moldrug.fitness.CostMultiReceptors`. It uses all the response variables.
- {py:func}`moldrug.fitness.CostMultiReceptorsOnlyVina`. Only use the information of the Vina scores.

```python
from moldrug import fitness
cost = fitness.Cost
```

A molecule must present several properties to be considered a drug. Some of the most important are: be potent, reach the biological target (good ADME profile) and be real. The last obvious property could be a bottleneck for computer-assisted drug design. Because we want to optimize several response variables at the same time. These cost functions (with the exception of {py:func}`moldrug.fitness.CostOnlyVina`) use the concept of [desirability functions](https://www.sciencedirect.com/science/article/pii/S0169743911000797) which optimizes several response variables on the hub.

The following are the response variables:

- [Vina score](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3041641/)
- [Quantitative Estimation of Drug-likeness (QED)](https://www.nature.com/articles/nchem.1243)
- [Synthetic accessibility score](https://jcheminf.biomedcentral.com/articles/10.1186/1758-2946-1-8)

For each of this response variable we create their corresponding [Derringer-Suich desirability functions](https://www.tandfonline.com/doi/abs/10.1080/00224065.1980.11980968) (check [here](https://moldrug.readthedocs.io/en/latest/notebooks/desirability.html) some examples of how it is implemented in moldrug). And then we combine them as a geometric mean:

```{math}
D = {\left[\prod_{i = 1}^{3} d_i^{w_i}\right]}^{\frac{1}{\sum_i w_i}}
```

where {math}`w_i` are the weights of each variable; and {math}`d_i` the desirability functions. Each individual {math}`d_i` ranges from 0 to 1 and therefore also {math}`D`. Because we are looking for the minimum, the function `cost` return {math}`1 - D`.

## Multi Receptor

Could be that our receptor presents high flexibility or that we are interested in generating specific small molecules. In this case could be convenient to add more than one receptor to the cost function. In {py:mod}`moldrug.fitness` module the cost functions {py:func}`moldrug.fitness.CostMultiReceptors` and {py:func}`moldrug.fitness.CostMultiReceptorsOnlyVina` try to reach this goal. For the case of flexibility, we could perform docking in an ensemble of protein structures and just keep the lower scoring rather than include all of them in the final desirability function.
