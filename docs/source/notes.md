# 📝 Notes

## CReM database

```{warning}
The precompiled CReM databases are distributed on [Zenodo](https://doi.org/10.5281/zenodo.16909329) and listed on the [CReM page](http://www.qsar4u.com/pages/crem.php). The older `replacements02_sc2.db.gz` from `qsar4u.com/files/cremdb` is no longer available. For the tutorials, the smallest one, `chembl22_sa2_hac12.db.gz` (159 MB), is a reasonable choice. For real problems consider reading first the [CReM paper](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-020-00431-w) and also see the discussion on [moldrug repo](https://github.com/ale94mleon/moldrug/discussions/6)

You can also build your own database from a SMILES file with CReM's `cremdb_create`. That is what moldrug's test suite does, see `tests/build_crem_test_db.py`; note that such a small database is only useful for testing.
```

## Explicit hydrogens and RDKit library

```{warning}
Many RDKit functions, including the `sascore`, operate optimally when applied to molecules without explicit hydrogens. The `sascore` assigns a score ranging from 1 to 13, where 1 indicates an easily synthesizable molecule, and 13 signifies a challenging one. However, when explicit hydrogens are present, the function tends to yield higher scores, incorrectly suggesting a greater difficulty of synthesis. Detailed discussions on this issue can be found at [SA_Score gives different values depending if the molecule has explicit hydrogens](https://github.com/rdkit/rdkit/discussions/7047).

The moldrug-3.7.0 release has addressed this concern for all built-in fitness functions. Refer to [moldrug's CHANGELOG](https://moldrug.readthedocs.io/en/latest/source/CHANGELOG.html) for a comprehensive overview of the fixes and improvements.
```
