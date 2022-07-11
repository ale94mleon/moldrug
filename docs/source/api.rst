Summary
=======

Input data
----------

Currently, MolDrug only accepts valids RDKit SMILES and valids pdbqt files that
will be procesed by Vina.

The idea
--------
The idea is simple: we want the best drug for some application. The problem is
that it is not simple as it sounds. The chemical space is enormous and the estimation
of the fitness is also a term of concern.

Here we address this problem using a GA. The inputs are:

#. Biological target (pdbqt format).
#. SMILES string of some known (or guessing) drug.
#. Definition of the binding pocket.
#. Definition of the GA running variables.

With the initial SMILES, a random population of ``popsize``
individuals will be generated through the `CReM <https://github.com/DrrDom/crem>`_ 
operation ``mutate_mol``. Every individual will be evaluated through the fitness function.
The best individual will be used for crossover and mutation, generating ``pc*popsize`` offsprings.
The offsprings will be merge with the current population, sorted respect to the fitness function
and the first ``popsize`` individuals will be kept for the next generation.
This cycle will run for ``maxiter`` generations.

Fitness functions
-----------------

The default fitness function could be access through::

```python
from moldrug import fitness
cost = fitness.Cost
```
A molecule must present several properties to be considered a drug. Some of the must important are: be potent, reach the biological target (good ADME profile) and be real. The last obvious property could be a bottleneck for computer assisted drug design. Because we want to optimize several response variables at the same time; this `cost` function use the concept of desirability functions ([see this paper](https://www.sciencedirect.com/science/article/pii/S0169743911000797)) which optimize several response variables on the hub.

The following are the response variables:
1. Vina Scoring Function. [paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3041641/)
1. Quantitative Estimation of Drug-likeness. [paper](https://www.nature.com/articles/nchem.1243)
2. Synthetic accessibility. [paper](https://jcheminf.biomedcentral.com/articles/10.1186/1758-2946-1-8)

For each of this response variables we create their corresponded Derringer-Suich desirability functions. And them we combine them as a geometric mean:

$$
D = {\left[\prod_{i = 1}^{3} d_i^{w_i}\right]}^{\frac{1}{\sum_i w_i}}
$$

where $w_i$ are the weights of each variable; and $d_i$ the desirability functions. Each individual $d_i$ ranges from 0 to 1 and therefore also $D$.
Because we are looking for the minimum, the function `cost` return $ 1 - D$.

### Multi Receptor
Could be that our receptor presents high flexibility or that we are interested in generate specific small molecules. In this case could be convenient to add more than one receptor to the cost function. In the `fitness` module the cost function `CostMultiReceptors` tries to reach this goal. For the case of flexibility, we could perform docking in an ensemble of protein structures, and just keep the lower scoring rather that included all of them in the final desirability function.