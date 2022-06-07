# lead
```
git config --global init.defaultBranch main
git init
git add README.md
git commit -m "first commit"
git branch -M main
git remote add origin https://github.com/ale94mleon/lead.git
git commit -m "first commit"
```


# The idea
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