"""Execution backends used to evaluate a population.

:meth:`moldrug.opt.GA.__call__` and :meth:`moldrug.opt.Local.__call__` evaluate the cost
function of every individual of a generation independently, which is where moldrug
parallelizes across the population. This module holds the backends that carry out that
evaluation, from a plain serial loop to a Dask cluster spanning several machines.

The usual entry point is the ``njobs`` parameter of those two methods, which is shorthand
for :class:`Runner` in :attr:`RunnerMode.MULTIPROCESSING` mode. Pass a :class:`Runner`
explicitly to use any of the other modes.

The Dask modes need optional dependencies, installed as the ``cluster`` extra::

    pip install moldrug[cluster]

See the parallelization section of the documentation for how this level relates to the
``ncores`` given to AutoDock-Vina and CReM inside a single cost evaluation.
"""
import multiprocessing as mp
from enum import Enum
from typing import Callable, List, Optional

import tqdm

try:
    from dask.distributed import Client

    dask_available = True
except ImportError:
    dask_available = False


def serial_execution(func: Callable, entries: List):
    """Evaluate every entry one after another, in the current process.

    Parameters
    ----------
    func : Callable
        The function to apply, typically the cost function of an optimization.
    entries : List
        The arguments to call ``func`` with, one per individual.

    Returns
    -------
    list
        The results, in the order of ``entries``.
    """
    return [func(args) for args in tqdm.tqdm(entries, total=len(entries))]


def parallel_execution_multiprocessing(func: Callable, entries: List, process_count: int):
    """Evaluate the entries on a pool of local processes.

    Parameters
    ----------
    func : Callable
        The function to apply, typically the cost function of an optimization. It has to be
        picklable, as it is sent to the worker processes.
    entries : List
        The arguments to call ``func`` with, one per individual.
    process_count : int
        The number of worker processes in the pool.

    Returns
    -------
    list
        The results, in the order of ``entries``.
    """
    pool = mp.Pool(process_count)
    results = [entry for entry in tqdm.tqdm(pool.imap(func, entries), total=len(entries))]
    pool.close()

    return results


if dask_available:
    def parallel_execution_dask_local(func: Callable, entries: List, client: Client):
        """Evaluate the entries through a Dask client running on this machine.

        Parameters
        ----------
        func : Callable
            The function to apply, typically the cost function of an optimization.
        entries : List
            The arguments to call ``func`` with, one per individual.
        client : dask.distributed.Client
            The client to submit the work to.

        Returns
        -------
        list
            The gathered results, in the order of ``entries``.
        """
        return client.gather(client.map(func, entries))

    def parallel_execution_dask_cluster(func: Callable, entries: List, client: Client):
        """Evaluate the entries on the workers of a Dask cluster.

        Parameters
        ----------
        func : Callable
            The function to apply, typically the cost function of an optimization. Whatever
            it needs, AutoDock-Vina included, has to be available on every worker node.
        entries : List
            The arguments to call ``func`` with, one per individual.
        client : dask.distributed.Client
            The client connected to the cluster.

        Returns
        -------
        list
            The gathered results, in the order of ``entries``.
        """
        return client.gather(client.map(func, entries))


class RunnerMode(Enum):
    """How a :class:`Runner` evaluates a generation.

    Each mode takes its own parameters, which :class:`Runner` checks on construction.

    Attributes
    ----------
    SERIAL
        One entry after another, in the current process. Takes no parameters. Useful to
        debug a cost function or to profile a single evaluation, where the workers of the
        other modes would hide the traceback or the timings.
    MULTIPROCESSING
        A pool of local processes, sized by ``process_count``. This is what ``njobs``
        selects, and the right choice on a single machine.
    DASK_LOCAL
        A Dask client on this machine, with ``process_count`` workers of ``thread_count``
        threads each. Needs the ``cluster`` extra.
    DASK_JOB_QUEUE_CLUSTER
        The workers of the ``dask_cluster`` handed to the runner, so the generation is
        spread over the nodes the cluster submits jobs to. Needs the ``cluster`` extra.
    """
    SERIAL = "serial"
    MULTIPROCESSING = "multiprocessing"
    DASK_LOCAL = "dask_local"
    DASK_JOB_QUEUE_CLUSTER = "dask_job_queue_cluster"


class Runner:
    """Evaluates a generation with one of the :class:`RunnerMode` backends.

    Pass an instance to :meth:`moldrug.opt.GA.__call__` or
    :meth:`moldrug.opt.Local.__call__` through their ``runner`` parameter, instead of
    ``njobs``. The two are mutually exclusive.

    Examples
    --------
    .. code-block:: python

        from moldrug.runner import Runner, RunnerMode

        # The explicit form of ga(njobs=12)
        ga(runner=Runner(RunnerMode.MULTIPROCESSING, process_count=12))

        # Serial, to debug a cost function
        ga(runner=Runner(RunnerMode.SERIAL))

    On a cluster, build the cluster object first and hand it over:

    .. code-block:: python

        from dask_jobqueue import SLURMCluster

        cluster = SLURMCluster(queue='short', cores=12, processes=1,
                               memory='8GB', walltime='00:30:00')
        cluster.scale(4)
        ga(runner=Runner(RunnerMode.DASK_JOB_QUEUE_CLUSTER, dask_cluster=cluster))
    """
    mode: RunnerMode
    thread_count: Optional[int]
    process_count: Optional[int]
    dask_cluster = None
    client = None

    def __init__(self, mode: RunnerMode, thread_count: Optional[int] = None, process_count: Optional[int] = None, dask_cluster = None):
        """Constructor

        Parameters
        ----------
        mode : RunnerMode
            The backend to evaluate with.
        thread_count : int, optional
            Threads per worker, only for :attr:`RunnerMode.DASK_LOCAL`, by default None.
        process_count : int, optional
            Worker processes, for :attr:`RunnerMode.MULTIPROCESSING` and
            :attr:`RunnerMode.DASK_LOCAL`, by default None.
        dask_cluster : optional
            A Dask cluster object, for example a ``dask_jobqueue.SLURMCluster``. Only for
            :attr:`RunnerMode.DASK_JOB_QUEUE_CLUSTER`, by default None.

        Raises
        ------
        ImportError
            If a Dask mode is asked for without the ``cluster`` extra installed.
        ValueError
            If the parameters do not match the mode, or the mode is unknown.
        """
        if mode in [RunnerMode.DASK_LOCAL, RunnerMode.DASK_JOB_QUEUE_CLUSTER] and not dask_available:
            raise ImportError("Execution using Dask requires installation of optional dependencies. "
                              "The optional pip package group is called 'cluster'")

        if mode is RunnerMode.SERIAL:
            if not (thread_count is None and process_count is None and dask_cluster is None):
                raise ValueError("Serial execution doesn't take any parameters.")
        elif mode is RunnerMode.MULTIPROCESSING:
            if not (thread_count is None and process_count is not None and dask_cluster is None):
                raise ValueError("Only process count is needed.")
        elif mode is RunnerMode.DASK_LOCAL:
            if not ((thread_count is not None or process_count is not None) and dask_cluster is None):
                raise ValueError("Only process count and/or thread count are needed.")
        elif mode is RunnerMode.DASK_JOB_QUEUE_CLUSTER:
            if not (thread_count is None and process_count is None and dask_cluster is not None):
                raise ValueError("Dask execution takes only a Dask cluster object.")
        else:
            raise ValueError(f"Unknown execution mode: {mode}")

        self.mode = mode
        self.thread_count = thread_count
        self.process_count = process_count
        self.cluster = dask_cluster

        if mode in [RunnerMode.DASK_LOCAL]:
            self.client = Client(n_workers=process_count, threads_per_worker=thread_count)

        if mode in [RunnerMode.DASK_JOB_QUEUE_CLUSTER]:
            self.client = Client(dask_cluster)

    def run(self, func: Callable, entries: List):
        """Evaluate the entries with the configured backend.

        Parameters
        ----------
        func : Callable
            The function to apply, typically the cost function of an optimization.
        entries : List
            The arguments to call ``func`` with, one per individual.

        Returns
        -------
        list
            The results, in the order of ``entries``.
        """
        if self.mode == RunnerMode.SERIAL:
            return serial_execution(func, entries)
        elif self.mode == RunnerMode.MULTIPROCESSING:
            return parallel_execution_multiprocessing(func, entries, process_count=self.process_count)
        elif self.mode == RunnerMode.DASK_LOCAL:
            return parallel_execution_dask_local(func, entries, self.client)
        elif self.mode == RunnerMode.DASK_JOB_QUEUE_CLUSTER:
            return parallel_execution_dask_cluster(func, entries, self.client)
        else:
            raise ValueError(f"Unknown execution mode: {self.mode}")
