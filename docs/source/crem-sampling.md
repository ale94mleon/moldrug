# 🎲 CReM sampling and exploration

{py:meth}`moldrug.opt.GA.mutate` asks CReM for the mutants of one molecule and keeps a
single one of them. One CReM keyword decides how that is done, and it is worth understanding
before changing it, because it is not only a speed setting.

```python
mutate_crem_kwargs = {
    'radius': 3,
    'max_replacements': None,   # moldrug's default: build every mutant, then pick one
}
```

| | what CReM does | what moldrug does with it |
|---|---|---|
| `max_replacements = None` | builds **every** possible mutant | picks one at random, discards the rest |
| `max_replacements = N` | builds **at most N**, sampled | picks one at random from those |

Capping is the faster of the two, and the bigger the database the more it saves. It also
searches differently. This page explains how, with numbers you can reproduce.

## How CReM builds the mutants of a molecule

CReM **cuts the molecule into fragments**: every piece that could be removed, each with the
chemical context around it. For each cut, it looks up in the database every fragment that has
been seen in the same context, and swaps them in. Each swap is one product.

The cuts of a molecule are wildly unequal in how much the database has to offer for them. For
the x0161 ligand against the small database the test suite builds, counting which cut each of
the 151 products came from:

| fragment removed | products it yields | share |
|---|---|---|
| `NS(=O)(=O)[*:1]` | 120 | 79.5% |
| `OCc1ccccc1[*:1]` | 14 | 9.3% |
| `O[*:1]` | 9 | 6.0% |
| `OC[*:1]` | 5 | 3.3% |
| `c1ccc([*:1])c(C[*:2])c1` | 3 | 2.0% |

Replacing the sulfonamide accounts for four out of five of all possible mutants, simply
because that group is common in the database and has many known substitutes. The ring cut at
the bottom accounts for three mutants in total.

```{warning}
That table is **not** the list of places this molecule can be cut. CReM keeps a set of the
products it has already returned and skips any repeat, so a molecule reachable by cutting in
two different places is credited to whichever cut was processed first. Enumeration walks the
cuts in a fixed order, so a cut whose products were all already emitted by an earlier one is
credited with nothing and never appears.

This molecule has **at least 8** cut points, not 5. Three of them, `O=S(=O)([*:2])[*:1]`,
`N[*:1]` and `c1ccc([*:2])c([*:1])c1`, are credited with 0 products here and surface only
once sampling shuffles the order. They are smaller cuts inside the sulfonamide and the ring,
reaching overlapping chemistry.

The product set itself is unambiguous, all 151 of them, and so is the share of mutations a
strategy spends on each cut. Only the attribution of one molecule to one cut is
order-dependent. This is why the tables on this page do not all list the same cuts.
```

## What sampling changes

With `max_replacements = N`, CReM shuffles the cuts, gives each one a quota of
`N // number_of_cuts` (at least 1), fills that quota with replacements drawn at random from
the database, and stops once it has `N` products. Leftovers are pooled and topped up from at
the end.

**The quota is per cut, and it ignores how many replacements that cut actually has.** A cut
with 120 available replacements and a cut with 3 get the same quota, so the products of one
call are spread evenly over the cuts rather than following the database. `N = 1` is the
extreme: a call stops at the first cut of the shuffle, so one randomly chosen cut supplies
the single product. The shuffle is redrawn per call, so over many mutations the cuts still
come up about equally often.

Where the mutations land, as a share of 400 mutations of the same molecule:

| fragment removed | enumerate all | cap 1 | cap 10 | cap 100 |
|---|---|---|---|---|
| `NS(=O)(=O)[*:1]` | **77.2%** | **14.3%** | 29.5% | 59.0% |
| `OCc1ccccc1[*:1]` | 10.8% | 12.6% | 14.0% | 10.5% |
| `O[*:1]` | 7.8% | 11.1% | 9.5% | 9.5% |
| `O=S(=O)([*:2])[*:1]` | 0.0% | 13.8% | 12.0% | 7.8% |
| `OC[*:1]` | 1.8% | 10.3% | 11.2% | 5.8% |
| `c1ccc([*:1])c(C[*:2])c1` | **2.5%** | **13.8%** | 7.5% | 1.2% |
| `N[*:1]` | 0.0% | 10.8% | 9.2% | 5.0% |
| `c1ccc([*:2])c([*:1])c1` | 0.0% | 13.3% | 7.0% | 1.2% |


Enumerating follows the database: the sulfonamide gets 77% of the attention because that is
where 79.5% of the possible mutants are. Capping flattens the column: at `N = 1` every cut
gets roughly an eighth, whatever it has to offer.

## What that costs

### Exploration

The experiment: take the x0161 ligand and mutate it **400 separate times**, which is what the
GA does whenever the same parent keeps being selected across generations. Each mutation
returns **one** child. Against this database there are **151 possible children**.

So the question is: over those 400 mutations, how much of the available chemical space did we
visit, and how often did we get the same molecule back?

| `max_replacements` | different children seen | the single most common child | the 3 most common together | mutations that returned nothing |
|---|---|---|---|---|
| `None` (enumerate all) | **140** of 151 | 1.8% of the 400 | 5.2% | 0 |
| 100 (cap 100) | 142 of 151 | 2.2% | 6.0% | 0 |
| 10 (cap 10) | 111 of 151 | 8.0% | 15.5% | 0 |
| 1 (cap 1) | **92** of 151 | **13.8%** | **29.6%** | 2 |

Reading the two extremes:

- **Enumerating everything**, 400 mutations produced **140 different molecules** out of the
  151 that exist, and no molecule came up more than 7 times. The search spread over nearly
  the whole space. Similar result is obtained at cap 100, but this is close to the total max of 151.
- **`max_replacements = 1`**, the same 400 mutations produced only **92 different molecules**,
  one of them came back 55 times, and three of them together accounted for 118 of the 398
  mutations that returned anything.

Columns 3 and 4 measure concentration. If all 151 children were equally likely, any one
would come up 0.7% of the time and any three 2.0%. The `None` row sits close to that; the `1`
row is twenty times above it in the third column and fifteen times in the fourth.

Note also the last column: `max_replacements = 1` returned nothing at all in 2 of 400 calls,
the single sampled replacement having been filtered out downstream.
{py:meth}`moldrug.opt.GA.mutate` then returns the parent unchanged and the offspring slot is
wasted. Any value above 1 avoids this.

### Where the repeats come from

The repeated molecules are the children of the cuts that have almost nothing to offer.
Grouping those same 400 mutations at `max_replacements = 1` by the cut CReM reported:

| fragment removed | mutations from it | different children it gave | its most repeated child |
|---|---|---|---|
| `NS(=O)(=O)[*:1]` | 57 | 48 | 3 times |
| `O=S(=O)([*:2])[*:1]` | 55 | 13 | 7 times |
| `c1ccc([*:1])c(C[*:2])c1` | 55 | 3 | 20 times |
| `c1ccc([*:2])c([*:1])c1` | 53 | 2 | 29 times |
| `OCc1ccccc1[*:1]` | 50 | 17 | 7 times |
| `O[*:1]` | 44 | 9 | 9 times |
| `N[*:1]` | 43 | 6 | 13 times |
| `OC[*:1]` | 41 | 7 | 9 times |

Every cut gets a similar number of mutations, between 41 and 57, which is the quota working
as designed. But one of them turns its 57 mutations into 48 different molecules, while
another turns its 53 into **two**. A cut that can only give two children has no choice but to
hand the same two back, 29 times and 24 times.

That is the whole mechanism:

> even attention across cuts **+** uneven richness per cut **=** repetition concentrated on
> the poorest cuts.

## What this means for moldrug

Read the other way round, CReM's sampling is a **novelty-seeking** strategy, and a sensible
one. Giving every cut the same quota deliberately rescues the rare chemistry: the positions
your database barely covers get as much attention as the ones it knows inside out. For
enumerating a diverse set of ideas from a single molecule, that is exactly what you want.

moldrug's loop is a different situation, and the way it pays for repeats is worth being
precise about. {py:meth}`moldrug.opt.GA.mutate` is called `nc` times per generation, and a
child is kept only if it has not been seen before:

```python
if children not in self.SawIndividuals and children not in popc and children.pdbqt:
    popc.append(children)
```

**A repeated molecule is not costed twice.** It is recognised as already seen and dropped
before the cost function runs, so no docking is repeated. What is lost is the attempt
itself: that slot of the generation produced nothing, and the generation ends with fewer new
candidates than it asked for.

That is the real price. The more often a small `max_replacements` makes CReM return a
molecule the run has already seen, the more generations end with fewer new children than
`nc`, and the fewer distinct molecules the run has explored after the same amount of work.
The search advances more slowly per generation, and a population fed by a thinner stream of
new candidates is more exposed to **premature convergence**: less novelty enters, so the
population has less to recombine away from whatever it has already found.

The second effect compounds it. The well-covered positions are squeezed at the same time:
replacing the sulfonamide falls from 77% of mutations to 14%, so the hundred-odd plausible,
database-supported substitutions there are reached far more slowly, while attempts pile up on
cuts that have only two or three children to give.

So the trade is not "faster, slightly different". It is a choice between spending the mutation
budget evenly across the molecule, which wastes attempts on positions that have nothing left
to offer, or spending it in proportion to what the database actually supports. Neither is
wrong. For a one-shot enumeration of ideas from one molecule the first is attractive; for an
iterative search that needs a steady supply of new candidates, the second keeps the
generations full.

**That is why moldrug leaves `max_replacements` unset by default.** Set it when the cost of
enumerating is what is hurting you, and prefer a value large enough that generations stay
full, rather than the smallest one that runs fast.

```{warning}
**Every number on this page was measured against the small database the test suite builds,
made from 1500 molecules. None of it has been measured against a production database.**

How much a given cap flattens the search depends on how many replacements the database holds
per cut, and that is a property of the database, not of moldrug. In the tables above
`max_replacements = 100` looks almost identical to enumerating everything, but only because
the cap barely binds here: the molecule has 151 mutants in total, so 100 is hardly a cap.

The *expectation*, not a measurement, is that a database built from ChEMBL or Enamine holds
far more replacements per cut, so the same cap would bind much harder and behave closer to
the bottom rows of these tables. That is an inference from their size, hundreds of megabytes
compressed against about ten here, and it may be wrong: the number of replacements per cut
and how unevenly they are spread are both empirical questions, and the answer may differ by
molecule as well as by database.

**Measure on the database you actually use**, with the script below, before choosing a value.
```

## Speed

Same molecule, same database, median of 5 calls:

| `max_replacements` | time (s) | mutants built | speedup |
|---|---|---|---|
| None (enumerate all) | 0.150 | 151 | 1.0x |
| 1 | 0.073 | 1 | 2.1x |
| 10 | 0.078 | 9 | 1.9x |
| 100 | 0.118 | 84 | 1.3x |

The gains here are modest because the database is tiny and 151 molecules are cheap to build.
The setting is meant for the case this table cannot show, a database where each mutation
builds far more molecules than are needed to keep one, but how much is saved there has not
been measured and will depend on the database.

## Reproducing these numbers

Every table on this page comes from
[`Contrib/crem_sampling/crem_sampling.py`](https://github.com/ale94mleon/moldrug/blob/main/Contrib/crem_sampling/crem_sampling.py):

```bash
cd Contrib/crem_sampling
python ../../tests/build_crem_test_db.py --output crem.db
python crem_sampling.py --db crem.db
```

Against a real database, from [Zenodo](https://doi.org/10.5281/zenodo.16909329):

```bash
python crem_sampling.py --db chembl33_sa2_f5.db --draws 200
```

It takes `--smiles` to change the molecule and `--caps` to change which values are compared,
and prints every table as Markdown, ready to paste back here.

## If you do set it: seeding

CReM samples with a random generator of its own, `random.Random(seed)`, which does **not**
read the global random state that moldrug seeds through `randomseed`. moldrug therefore draws
a fresh seed per call from the global generator. There are three ways `seed` can reach CReM,
and two of them look reasonable but are wrong.

1. **A fixed value on every call.** The generator is rebuilt identically each time, so the same
parent yields the same child on every mutation. A parent selected across several generations
is mutated into one molecule, over and over, and the mutation operator stops exploring.
Reproducible, but not random.

2. **`seed` left as `None`.** CReM seeds itself from the operating system. The mutations do vary
within a run, so nothing looks wrong, but the run cannot be reproduced: the same
`randomseed` gives different molecules every time. Random, but not reproducible.

3. **A fresh seed per call, drawn from the global generator.** This is what moldrug does from
version 4.0.0. Earlier versions never passed a `seed`, so anyone who set
`max_replacements` themselves got the second behaviour. The
seed changes from one mutation to the next, so the children vary, and the sequence of seeds
is itself determined by `randomseed`, so the whole run repeats. Random *and* reproducible.

If you put your own `seed` in `mutate_crem_kwargs`, moldrug leaves it alone and you get the
first behaviour.

Sampling with a seed is only available from CReM 0.3.0. It closed
[CReM issue #12](https://github.com/DrrDom/crem/issues/12), *"Reproducible output"*, open
since 2022.
