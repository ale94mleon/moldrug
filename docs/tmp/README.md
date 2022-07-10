![logo](https://github.com/ale94mleon/moldrug/blob/main/row_data/logo.png?raw=true)

[![Documentation Status](https://readthedocs.org/projects/moldrug/badge/?version=latest)](https://moldrug.readthedocs.io/en/latest/?badge=latest)
[![PyPi license](https://badgen.net/pypi/license/moldrug/)](https://pypi.python.org/pypi/moldrug/)
![Tests](https://github.com/ale94mleon/moldrug/actions/workflows/python-package-conda.yml/badge.svg)
[![PyPI version shields.io](https://img.shields.io/pypi/v/moldrug.svg)](https://pypi.python.org/pypi/moldrug/)
[![PyPI pyversions](https://img.shields.io/pypi/pyversions/moldrug.svg)](https://pypi.python.org/pypi/moldrug/)
[![Downloads](https://static.pepy.tech/personalized-badge/moldrug?period=month&units=international_system&left_color=grey&right_color=brightgreen&left_text=Downloads)](https://pepy.tech/project/moldrug)

# Table of content
1.  [moldrug](#moldrug)
    1.  [The idea](#The-idea)
    2.  [Fitness functions](#Fitness-functions)
        1.  [Multi Receptor](#Multi-Receptor)
    3.  [Example of use](#Example-of-use)
        1.  [Saving the data](#Saving-the-data)
            1.  [Saving intermediate solution](#Saving-intermediate-solution)
            2.  [Exporting a DataFrame](#Exporting-a-DataFrame)
    4.  [Global, local and local-customize optimization](#Global,-local-and-local-customize-optimization)         

# moldrug

moldrug is a python package for moldrug generation and optimization of small molecules. It use a Genetic Algorithm (GA) as searching engine in the chemical space and CReM library ([crem](https://github.com/DrrDom/crem)) as chemical structure generator.

## Installation instruction
```bash
conda create -y -n moldrug
conda activate moldrug
conda install -y -c conda-forge rdkit">=2022.0"
conda install -y -c conda-forge openbabel">=3.1.0"
conda install -y -c bioconda autodock-vina
# To get the version on developing
pip install git+https://github.com/ale94mleon/moldrug.git@main
# To get the last "stable" version. This project is still in beta state.
pip install moldrug
```
In this way you will have a completely functional moldrug environment. It is needed through conda in order to get RDKit and OpenBabel, which have a non tribal installation through pip. 

## The idea

The idea is simple: we want the best drug for some application. The problem is that it is not simple as it sounds. The chemical space is enormous and the estimation of the fitness is also a term of concern.

Here we address this problem using a GA. The inputs are:
1. Biological target (pdbqt format).
2. SMILES string of some known (or guessing) drug.
3. Definition of the binding pocket.
4. Definition of the GA running variables.

With the initial SMILES, a random population of `popsize` individuals will be generated through the [crem](https://github.com/DrrDom/crem) operation `mutate_mol`. Every individual will be evaluated through the fitness function. The best individual will be used for crossover and mutation, generating `pc*popsize` offsprings. The offsprings will be merge with the current population, sorted respect to the fitness function and the first `popsize` individuals will be kept for the next generation. This cycle will run for `maxiter` generations.

## Fitness functions

The default fitness function could be access through:

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

## Example of use
```python
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from moldrug import fitness, utils
import json
from multiprocessing import cpu_count

receptor = 'x0161'
maxiter = 2
popsize = 2
njobs = 2

with open('data/box.json', 'r') as f:
    grid_opt = json.load(f)[receptor]['A']
with open('data/smi.json', 'r') as f:
    init_smiles = json.load(f)[receptor]

out = utils.GA(
    smiles=init_smiles,
    maxiter=maxiter,
    popsize=popsize,
    crem_db_path = '/provide/the/path/to/replacements02_sc2.5.db',
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
    costfunc = fitness.Cost,
    costfunc_kwargs = {
        'receptor_path': f'data/{receptor}.pdbqt',
        'boxcenter' : grid_opt['boxcenter'],
        'boxsize': grid_opt['boxsize'],
        'exhaustiveness': 8,
        'ncores': int(cpu_count() / njobs),
        'num_modes': 1,
    },
    )
out(njobs = njobs)
```
### Saving the data
After completion we could save the whole GA object with the `pickle` method. By default the compress option is set to False, but you could change it. To reopen the saved file use `utils.decompress_pickle` or `utils.loosen`  if the saved object was compressed or not respectively.
*   Non-Compressed
```python
out.pickle(f'result')
out = utils.loosen('result.pkl')
```
*   Compressed
```python
out.pickle(f'result', compress = True)
out = utils.decompress_pickle('result.pbz2')
```

If you would like to run in this object four more generations, you just need to set `maxiter` to 4 and call the class again. Even you can modify some parameters used in the searching (do not modify the cost function for the next call, it will not give errors but the results will not be correct)
```python
out.maxiter = 4
out.mutate_crem_kwargs = {
    'radius':5,
    'min_size':2,
    'max_size':5,
    'min_inc':-1,
    'max_inc':2,
    'ncores':2,
},
out(njobs = njobs)
```
The `out` instance has several useful attributes:
*   `out.pop`, the last (and the best) population of the simulation sorted by cost function. `out.pop[0]` will return the best `Individual`.
    *   Each `Individual` is a `utils.Individual` instance that have several attributes:
        *   `idx`: A unique identifier, that give the order of when the current Individual was generated during the simulation.
        *   `cost`: The cost calculated by the cost function (or np.inf by default)
        *   `mol`: A `Chem.rdchem.Mol` molecule.
        *   Other attributes that the cost function adds. In the case of `fitness.Cost`: `pdbqt`, `vina_score`, `qed` and `sa_score`.
*   `out.NumCalls`: Number of time that you call this instance. In this case, because we called twice, it should be 2.
*   `out.NumGen`: The total number of generations. During the first call we set `maxiter = 2` and in the second one `maxiter = 4`. Therefore, in total we made 6 generations. So, the result of this attribute should be 6.
*   `out.SawIndividuals` is a list of every generated `Individual` over the simulation.

#### Saving intermediate solution
Could be that for some reason the job is killed. In order to prevent loose all the information you could set in the initialization of the `GA` class the keywords:
*   `save_pop_every_gen`: Every how many generations the population will be saved.
*   `pop_file_name`: The name of the file to save the population. By default the extension .pkl is added.

Then you just need to initialize the `GA` class and give as population the saved one.
```python
from moldrug import utils
generation, init_pop = utils.loosen('pop.pkl')
# Initialize GA
out = utils.GA(...)
out.pop = init_pop
out(njobs = 3)
```
As you could imagine `out.pop` accept whatever you put there; but a correct input will be a list of `Individual` objects for which the corresponded `cost` function was passed. In other words, must have the `cost` attribute.

#### Exporting a DataFrame
In case that the simulation is really big (many Individuals in the population and many generations). Could be possible to export a simplified result.
```python
dataframe = out.to_dataframe()
```
This basically return the DataFrame of `out.SawIndividuals`. Every row is an `Individual`, and every columns the corresponded attributes. The attribute `mol` is deleted.

## Global, local and local-customize optimization
Walking on the chemical space is done through the method `mutate` of `GA` class tha is just a warper around `mutate_mol` function of [crem](https://github.com/DrrDom/crem) package. Depending on what is the final aim, will be the parameters that we will provided to `mutate` through the variable `mutate_crem_kwargs` in the initialization of `GA`.
If we know that the input smiles is not optimal, it should be convenient to take wider steeps in the chemical space. This could be accomplished with:
`min_size=1, max_size=8, min_inc=-5, max_inc=3`. In the other hand, if our interest is explore close to the input smiles, we should be more conservative and use `min_size=0, max_size=1, min_inc=-1, max_inc=1`. In this case the operation will be only: replace or delate one heavy atom, add one or two heavy atoms. In addition we could set `get_similar = True`. This flag doesn't ensure get similar molecules but get the most similar from the generated through the operation `mutate_mol` and a more conservative first population. Also we could add to the cost function the similarity as another response variable.

As standard, one possibility could be: a first run with wide steeps, and them call again the class but with more conservative parameters for the `mutate` method. In other words, a first run aiming global optimization and a second one aiming local optimization (this strategy was shown on [Saving the data](#Saving-the-data)).

For the "local-customize" optimization a new class is build in `utils.Local`. This class accept a molecule with the explicit Hs, cost function and parameters for the grow operation of [crem](https://github.com/DrrDom/crem). The results could be access with the attribute `pop` and also have similar method to save the data as `GA`. This class also use parallelization in the `__call__` method. 
 
## Brainstorm 
Here we intend to implement a Genetic Algorithm procedure for the moldrug optimization of chemical structures. In this way we are actively looking for a solution on the optimization problem.

The general idea is we give a starting ligand-protein complex. From there the ligand will be submitted to successive GA runs. For the GA the cost function could be any desirable property (LogP, affinity, minimum clashes in the binding pocket, etc...). In this first attend will be the Vina Scoring Function. Therefore a full docking without any restraint will be our cost function.

This runs will be done with the grow, mutate link operations of CReM. In this way we still ensure the accessibility of the generated compounds; which is (in my opinion), a drawback for the method presented on https://doi.org/10.1186/s13321-021-00501-7

But basically the issue here is how to implement in an efficient way the crossover and mutations?

What I have till now is the following:

For the mutation: I will use the grow and mutate method of CReM, but this ones have to be selected randomly as well how big could be the mutation
One idea is to go with higher substituent at the beginning and smaller one in the end (global to local optimization).

I have to keep track of the fragments that the each individual has in order do the crossover (the initial molecule could be considered as a fragments). Teh problem is that over the mutations and crossover will be more difficult to get track of this fragments But if we figure it out, then:

My initial idea for the crossover:
Cross the fragments and link them through the link option of CReM.

In addition, we could give more info to the mutate option, if we know (for some method) that some specific region is needed to improve the potency we could give this info to grow in this direction (the original idea). But know is plugged with a pseudo-GA optimization. I set pseudo because I had not clear how to code the crossover in an efficient way. Basically, the info to where grow will perform a selected mutation (grow in this case) therefore the convergency should speed up.

scscore and sascore are similar and give the same information.Therefore, for now will use sascore because is easy to get from rdkit.

We could think in the future add more than one conformer to the receptor. 