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