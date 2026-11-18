#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Measure what CReM's ``max_replacements`` does to speed and to exploration.

:meth:`moldrug.opt.GA.mutate` asks CReM for the mutants of one molecule and keeps a single
one. CReM can either enumerate every possible mutant, which is what moldrug does by
default, or sample at most ``max_replacements`` of them, which saves the work of building
mutants that are then discarded.

The two are not equivalent. CReM cuts the molecule into fragments and gives each one a quota
of ``max_replacements // number_of_fragments``, however many replacements the database holds
for it. Enumerating everything and picking one at random makes every product equally likely,
so the fragment that gets modified follows how many replacements it has; sampling makes the
fragments roughly equally likely instead. It re-weights which part of the molecule changes.

This script measures both: how long each strategy takes, and how concentrated the molecules
it returns are. Run it before changing ``max_replacements`` in ``mutate_crem_kwargs``, on
the database you actually use, because the effect depends on how uneven that database is.

Usage
-----
Against the small database the test suite builds::

    python ../../tests/build_crem_test_db.py --output crem.db
    python crem_sampling.py --db crem.db

Against a real one, from https://doi.org/10.5281/zenodo.16909329::

    python crem_sampling.py --db chembl33_sa2_f5.db --draws 200
"""
import argparse
import collections
import random
import statistics
import time

from crem.crem import mutate_mol
from rdkit import Chem

#: x0161, one of the ligands moldrug ships as example data
DEFAULT_SMILES = 'OCc1ccccc1COC(=O)C=1C=CC(=CC1)S(=O)(=O)N'

#: The defaults of moldrug.opt.GA
CREM_KWARGS = dict(radius=3, min_size=1, max_size=8, min_inc=-5, max_inc=3,
                   return_mol=True, ncores=1)


def time_one_call(mol, db, repeats=5, **extra):
    """Median wall time of one ``mutate_mol`` call, and how many mutants it built."""
    timings, built = [], 0
    for _ in range(repeats):
        start = time.perf_counter()
        built = len(list(mutate_mol(mol, db, **CREM_KWARGS, **extra)))
        timings.append(time.perf_counter() - start)
    return statistics.median(timings), built


def draw_many(mol, db, draws, max_replacements, seed=1):
    """Mutate the same molecule ``draws`` times, the way the GA does.

    Returns
    -------
    tuple[collections.Counter, collections.Counter, int]
        How often each product was drawn, how often each source fragment was the one
        modified, and how many calls returned nothing at all.
    """
    kwargs = dict(CREM_KWARGS, return_rxn=True)
    random.seed(seed)
    counts, by_fragment, empty = collections.Counter(), collections.Counter(), 0
    for _ in range(draws):
        if max_replacements is None:
            products = list(mutate_mol(mol, db, **kwargs))
        else:
            # A seed per call: CReM samples with a generator of its own, so a fixed seed
            # would return the same molecule every time and no seed would make the run
            # unreproducible. See https://github.com/DrrDom/crem/issues/12
            products = list(mutate_mol(mol, db, max_replacements=max_replacements,
                                       seed=random.randint(0, 2**31 - 1), **kwargs))
        if not products:
            empty += 1
            continue
        _, rxn, product = random.choice(products)
        counts[Chem.MolToSmiles(product)] += 1
        by_fragment[rxn.split('>>')[0]] += 1
    return counts, by_fragment, empty


def fragment_breakdown(mol, db):
    """How many products each removable fragment of the molecule accounts for.

    Returns
    -------
    collections.Counter
        Keyed by the SMARTS of the fragment that was cut out.
    """
    return collections.Counter(
        rxn.split('>>')[0]
        for _, rxn in mutate_mol(mol, db, return_rxn=True,
                                 **{k: v for k, v in CREM_KWARGS.items() if k != 'return_mol'}))


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--db', required=True, help='path to a CReM fragment database')
    parser.add_argument('--smiles', default=DEFAULT_SMILES, help='molecule to mutate (default: x0161)')
    parser.add_argument('--draws', type=int, default=400, help='mutations per strategy (default: %(default)s)')
    parser.add_argument('--caps', type=int, nargs='+', default=[1, 10, 100],
                        help='values of max_replacements to compare (default: %(default)s)')
    args = parser.parse_args()

    mol = Chem.AddHs(Chem.MolFromSmiles(args.smiles))

    everything = {Chem.MolToSmiles(m) for _, m in mutate_mol(mol, args.db, **CREM_KWARGS)}
    total = len(everything)
    print(f"{args.db}: {total} distinct mutants of the seed molecule\n")

    fragments = fragment_breakdown(mol, args.db)
    print(f"Enumerating them credits {len(fragments)} cuts of the molecule, very unequally.")
    print("CReM skips a product it has already returned, so a molecule reachable by cutting in")
    print("two places is credited to whichever cut came first: this is not the full list of")
    print("cuts, and sampling below will shuffle that order and credit more of them.\n")
    print("| fragment removed | products it yields | share |")
    print("|---|---|---|")
    for frag, n in fragments.most_common():
        print(f"| `{frag}` | {n} | {100 * n / sum(fragments.values()):.1f}% |")
    print()

    print("Speed, median of 5 calls")
    print("| max_replacements | time (s) | mutants built | speedup |")
    print("|---|---|---|---|")
    full_time, full_built = time_one_call(mol, args.db)
    print(f"| None (enumerate all) | {full_time:.3f} | {full_built} | 1.0x |")
    for cap in args.caps:
        seconds, built = time_one_call(mol, args.db, max_replacements=cap, seed=1)
        print(f"| {cap} | {seconds:.3f} | {built} | {full_time / seconds:.1f}x |")

    print(f"\nExploration: the same molecule mutated {args.draws} separate times, one child each,")
    print(f"out of {total} possible children. If every child were equally likely, the single most")
    print(f"common would be {100 / total:.1f}% and the three most common together {100 * 3 / total:.1f}%.\n")
    print("| max_replacements | different children seen | most common child "
          "| 3 most common together | returned nothing |")
    print("|---|---|---|---|---|")
    caps = [None] + args.caps
    per_fragment = {}
    for cap in caps:
        counts, by_fragment, empty = draw_many(mol, args.db, args.draws, cap)
        per_fragment[cap] = (counts, by_fragment)
        drawn = sum(counts.values())
        top3 = counts.most_common(3)
        label = 'None (enumerate all)' if cap is None else str(cap)
        print(f"| {label} | {len(counts)} of {total} "
              f"| {100 * top3[0][1] / drawn:.1f}% "
              f"| {100 * sum(n for _, n in top3) / drawn:.1f}% "
              f"| {empty} |")

    # Which fragment CReM reported for each draw. Note that CReM keeps a set of products
    # and skips one it has already emitted, so a product reachable by cutting two different
    # fragments is credited to whichever was processed first, and sampling shuffles that
    # order. These shares are therefore about the strategy as a whole, not a canonical
    # attribution of a molecule to a fragment.
    # Rows of the table: every fragment seen under any cap, most modified overall first.
    total_by_fragment = collections.Counter()
    for _, by_fragment in per_fragment.values():
        total_by_fragment += by_fragment
    seen_fragments = [frag for frag, _ in total_by_fragment.most_common()]

    print("\nWhich part of the molecule got modified, as a share of all mutations")
    print("| fragment removed | " + " | ".join(('enumerate all' if c is None else f'cap {c}') for c in caps) + " |")
    print("|---|" + "---|" * len(caps))
    for frag in seen_fragments:
        row = []
        for cap in caps:
            _, by_fragment = per_fragment[cap]
            drawn = sum(by_fragment.values())
            row.append(f"{100 * by_fragment.get(frag, 0) / drawn:.1f}%")
        print(f"| `{frag}` | " + " | ".join(row) + " |")

    smallest = min(args.caps)
    counts, by_fragment = per_fragment[smallest]
    print(f"\nWhere the repeats come from, with max_replacements = {smallest}")
    print("| fragment removed | mutations from it | different children it gave | its most repeated child |")
    print("|---|---|---|---|")
    children = collections.defaultdict(collections.Counter)
    random.seed(1)
    for _ in range(args.draws):
        products = list(mutate_mol(mol, args.db, max_replacements=smallest,
                                   seed=random.randint(0, 2**31 - 1),
                                   **dict(CREM_KWARGS, return_rxn=True)))
        if not products:
            continue
        _, rxn, product = random.choice(products)
        children[rxn.split('>>')[0]][Chem.MolToSmiles(product)] += 1
    for frag in sorted(children, key=lambda f: -sum(children[f].values())):
        c = children[frag]
        print(f"| `{frag}` | {sum(c.values())} | {len(c)} | {c.most_common(1)[0][1]} times |")
    print("\nA fragment that can only give a handful of children has to repeat them.")

if __name__ == '__main__':
    main()
