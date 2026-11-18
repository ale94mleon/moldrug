# CReM sampling

`crem_sampling.py` measures what CReM's `max_replacements` does to the speed of a mutation
and to the way the chemical space is explored, on whichever fragment database you give it.

Capping `max_replacements` saves the work of building mutants that are then discarded, but
it is not only a speed
setting: CReM gives each fragment of the molecule an equal quota, however many replacements
the database holds for it, so sampling re-weights which part of the molecule gets modified.
Run this before changing it.

```bash
python ../../tests/build_crem_test_db.py --output crem.db   # or use a real database
python crem_sampling.py --db crem.db
```

The results, and what they mean for `mutate_crem_kwargs`, are written up in the CReM
sampling page of the [documentation](https://moldrug.readthedocs.io).
