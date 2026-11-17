# ⚡ Parallelization

moldrug parallelizes at two independent levels, and it is worth being explicit about which
one you are tuning, because they multiply.

| Level | What it spreads | How you set it |
|---|---|---|
| Inside one cost evaluation | the work of a single docking or mutation | `ncores`, inside `costfunc_kwargs` and `mutate_crem_kwargs` |
| Across the population | the cost evaluations of a generation | `njobs`, or a {py:class}`moldrug.runner.Runner` |

The second level is where the {py:mod}`moldrug.runner` module comes in, and it is the one
that can reach beyond a single machine.

## Inside one cost evaluation

The built-in cost functions in {py:mod}`moldrug.fitness` pass `ncores` straight to
AutoDock-Vina as `--cpu`, so Vina threads a single docking across that many cores:

```python
costfunc_kwargs = {
    'receptor_pdbqt_path': 'receptor.pdbqt',
    'boxcenter': [0, 0, 0],
    'boxsize': [20, 20, 20],
    'ncores': 4,          # Vina uses 4 cores for each individual
}
```

CReM takes a separate `ncores` when it generates mutants:

```python
mutate_crem_kwargs = {
    'radius': 3,
    'ncores': 4,          # CReM uses 4 cores when mutating
}
```

These are unrelated to each other and to `njobs`.

## Across the population

{py:meth}`moldrug.opt.GA.__call__` and {py:meth}`moldrug.opt.Local.__call__` evaluate the
individuals of a generation independently, which is the natural place to parallelize.

The simple way is `njobs`:

```python
ga(njobs=12)
```

That is shorthand: internally moldrug builds
`Runner(RunnerMode.MULTIPROCESSING, process_count=njobs)` and runs the generation through a
multiprocessing pool. With the default `njobs=1` you get a pool of one, so the evaluations
run one after another.

For anything else, pass a {py:class}`moldrug.runner.Runner` instead.

```{warning}
`njobs` and `runner` are mutually exclusive. Passing `njobs > 1` together with a `runner`
raises an `AssertionError`; choose one.
```

### The execution modes

{py:class}`moldrug.runner.RunnerMode` offers four modes, each taking different parameters:

| Mode | Parameters | Use it for |
|---|---|---|
| `SERIAL` | none | debugging, and profiling a single cost evaluation |
| `MULTIPROCESSING` | `process_count` | one machine; this is what `njobs` gives you |
| `DASK_LOCAL` | `process_count` and/or `thread_count` | one machine, through the Dask scheduler |
| `DASK_JOB_QUEUE_CLUSTER` | `dask_cluster` | many machines, through a job queue |

```python
from moldrug.runner import Runner, RunnerMode

# Serial, for debugging
ga(runner=Runner(RunnerMode.SERIAL))

# The explicit form of ga(njobs=12)
ga(runner=Runner(RunnerMode.MULTIPROCESSING, process_count=12))

# Dask on this machine: 4 workers, 2 threads each
ga(runner=Runner(RunnerMode.DASK_LOCAL, process_count=4, thread_count=2))
```

### Distributed execution

Reaching more than one machine means handing a Dask cluster object to the runner. The
cluster describes the jobs the scheduler will submit on your behalf:

```python
from dask_jobqueue import SLURMCluster
from moldrug.runner import Runner, RunnerMode

cluster = SLURMCluster(
    queue='short',
    cores=12,
    processes=1,
    memory='8GB',
    walltime='00:30:00',
    job_extra_directives=[]
)
cluster.scale(4)          # ask for four jobs

runner = Runner(RunnerMode.DASK_JOB_QUEUE_CLUSTER, dask_cluster=cluster)

for _ in range(NumbCalls):
    ga(runner=runner)
```

The Dask modes need optional dependencies, installed as the `cluster` extra:

```bash
pip install moldrug[cluster]
```

Without them, constructing a `Runner` in either Dask mode fails with an `AssertionError`
pointing at this package set.

## Choosing the numbers

The two levels multiply: `njobs=12` with `ncores=4` asks for 48 cores. Going past what the
machine has does not fail, it just makes everything slower, as the processes fight over the
same cores.

A reasonable starting point on a single machine is to spend your cores on the population
rather than on each docking, since the population level scales more cleanly:

```python
from multiprocessing import cpu_count

njobs = cpu_count()
costfunc_kwargs['ncores'] = 1
```

Two things to keep in mind when you tune away from that:

- `njobs` above `popsize` buys nothing. There are only ever `popsize` cost evaluations in a
  generation, so the extra workers sit idle.
- Vina's threading does not scale linearly. Raising `ncores` helps a single docking far less
  than running more dockings at once, so prefer increasing `njobs` while there are still
  individuals left to evaluate.

On a cluster the same arithmetic applies per job: the `cores` you give `SLURMCluster` has to
cover the `ncores` each evaluation asks for.
