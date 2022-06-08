# lead

lead is a python package for lead generation and optimization of small molecules. It use a Genetic Algorithm (GA) as searching engine in the chemical space and CReM library ([crem](https://github.com/DrrDom/crem)) as chemical structure generator.

## The idea

The idea is simple: we want the best drug for some application. The problem is that it is not simple as it sounds. The chemical space is enormous and the estimation of the fitness is also a term of concern.

Here we address this problem using a GA. The inputs are:
1. Biological target (pdbqt format).
2. SMILES string of some known (or guessing) drug.
3. Definition of the binding pocket.
4. Definition of the GA running variables.

With the initial SMILES, a random population of `popsize` individuals will be generated through the crem operation `mutate_mol`. Every individual will be evaluated through the fitness function. The best individual will be used for crossover and mutation, generating `pc*popsize` offsprings. The offsprings will be merge with the current population, sorted respect to the fitness function and the first `popsize` individuals will be kept for the next generation. This cycle will run for `maxiter` generations.

## Fitness functions

The default fitness function could be access through:

```python
import lead as pb
cost = pb.Cost
```
A molecule must present several properties to be considered a drug. Some of the must important are: be potent, reach the biological target (good ADME profile) and be real. The last obvious property could be a bottleneck for computer assisted drug design. Because we want to optimize several response variables at the same time; this `cost` function use the concept of desirability functions ([see this paper](https://www.sciencedirect.com/science/article/pii/S0169743911000797)) which optimize several response variables on the hub.

The following are the response variables:
1. Vina Scoring Function. [paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3041641/)
1. Quantitative Estimation of Drug-likeness. [paper](https://www.nature.com/articles/nchem.1243)
2. Synthetic accessibility. [paper](https://jcheminf.biomedcentral.com/articles/10.1186/1758-2946-1-8)

For each of this response variables we create their corresponded Derringer-Suich desirability functions. And them we combined as a geometric mean:

$$
D = {\left[\prod_{i = 1}^{3} d_i^{w_i}\right]}^{\frac{1}{\sum_i w_i}}
$$

where $w_i$ are the weights of each variable; and $d_i$ the desirability functions. Each individual $d_i$ ranges from 0 to 1 and therefore also $D$.
Because we are looking for the minimum, the function `cost` return $ 1 - D$.

## Example of use
```python
import lead as pb
import json

receptor = 'x0161'
maxiter = 3
popsize = 4

with open('data/box.json', 'r') as f:
    grid_opt = json.load(f)[receptor]['A']
with open('data/smi.json', 'r') as f:
    init_smiles = json.load(f)[receptor]

out = pb.ga.GA(
    smiles=init_smiles,
    maxiter=maxiter,
    popsize=popsize,
    crem_db_path = '/provide/the/path/to/replacements02_sc2.5.db',
    pc = 1,
    costfunc = pb.fitness.Cost,
    receptor_path =f'data/{receptor}.pdbqt',
    boxcenter = grid_opt['boxcenter'],
    boxsize = grid_opt['boxsize'],
    exhaustiveness = 8,
    vina_cpus = 3,
    num_modes = 1,
    )
out(njobs = 4)
for o in out.pop:
    print(o.smiles, o.cost)
out.pickle(f'model.pkl')
```

## Git commands
```bash
git config --global init.defaultBranch main
git init
git add README.md
git commit -m "first commit"
git branch -M main
git remote add origin https://github.com/ale94mleon/lead.git
git commit -m "first commit"
```


# Brainstorm
Here we intend to implement a Genetic Algorithm procedure for the lead optimization of chemical structures. In this way we are actively looking for a solution on the optimization problem.

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

scscore and sascore are similar and give the same information.Therefore, for now will use sascore because is easy to get from rdkit