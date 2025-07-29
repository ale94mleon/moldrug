# 📝 Notes

## CReM database

```{warning}
For the tutorials, we will use the CReM database `replacements02_sc2.db.gz`. This is just because is a smaller one and therefore it has a fast download. For real problems consider reading first the [CReM paper](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-020-00431-w) and also see the discussion on [moldrug repo](https://github.com/ale94mleon/moldrug/discussions/6)
```

## Explicit hydrogens and RDKit library

```{warning}
Many RDKit functions, including the `sascore`, operate optimally when applied to molecules without explicit hydrogens. The `sascore` assigns a score ranging from 1 to 13, where 1 indicates an easily synthesizable molecule, and 13 signifies a challenging one. However, when explicit hydrogens are present, the function tends to yield higher scores, incorrectly suggesting a greater difficulty of synthesis. Detailed discussions on this issue can be found at [SA_Score gives different values depending if the molecule has explicit hydrogens](https://github.com/rdkit/rdkit/discussions/7047).

The moldrug-3.7.0 release has addressed this concern for all built-in fitness functions. Refer to [moldrug's CHANGELOG](https://moldrug.readthedocs.io/en/latest/source/CHANGELOG.html) for a comprehensive overview of the fixes and improvements.
```
