import multiprocessing as mp
from enum import Enum
from typing import Callable, List, Optional

import tqdm

try:
    from dask.distributed import Client
    from dask_jobqueue import JobQueueCluster

    dask_available = True
except ImportError:
    dask_available = False


def serial_execution(func: Callable, entries: List):
    return [func(args) for args in tqdm.tqdm(entries, total=len(entries))]


def parallel_execution_multiprocessing(func: Callable, entries: List, process_count: int):
    pool = mp.Pool(process_count)
    results = [entry for entry in tqdm.tqdm(pool.imap(func, entries), total=len(entries))]
    pool.close()

    return results


if dask_available:
    def parallel_execution_dask_local(func: Callable, entries: List, process_count: int, process_threads_count: int):
        client = Client(n_workers=process_count, threads_per_worker=process_threads_count)

        return client.gather(client.map(func, entries))


    def parallel_execution_dask_cluster(func: Callable, entries: List, cluster: JobQueueCluster):
        client = Client(cluster)

        return client.gather(client.map(func, entries))


class RunnerMode(Enum):
    SERIAL = "serial"
    MULTIPROCESSING = "multiprocessing"
    DASK_LOCAL = "dask_local"
    DASK_JOB_QUEUE_CLUSTER = "dask_slurm"


class Runner:
    mode: RunnerMode
    thread_count: Optional[int]
    process_count: Optional[int]
    dask_cluster = None

    def __init__(self, mode: RunnerMode, thread_count: Optional[int] = None, process_count: Optional[int] = None, dask_cluster = None):
        if mode in [RunnerMode.DASK_LOCAL, RunnerMode.DASK_JOB_QUEUE_CLUSTER]:
            assert dask_cluster is not None, "Execution using Dask requires installation of optional dependencies. The optional pip package group is called 'dask'"

        if mode is RunnerMode.SERIAL:
            assert thread_count is None and process_count is None and dask_cluster is None, "Serial execution doesn't take any parameters."
        elif mode in [RunnerMode.MULTIPROCESSING]:
            assert thread_count is None and process_count is not None and dask_cluster is None, "Only process count is needed."
        elif mode in [RunnerMode.DASK_LOCAL]:
            assert (thread_count is not None or process_count is None) and dask_cluster is None, "Only process count and/or thread count are needed."
        elif mode in [RunnerMode.DASK_JOB_QUEUE_CLUSTER]:
            assert thread_count is None and process_count is None and dask_cluster is not None, "Dask execution takes only a Dask cluster object."
        else:
            assert False

        self.mode = mode
        self.thread_count = thread_count
        self.process_count = process_count
        self.cluster = dask_cluster

    def run(self, func: Callable, entries: List):
        if self.mode == RunnerMode.SERIAL:
            return serial_execution(func, entries)
        elif self.mode == RunnerMode.MULTIPROCESSING:
            return parallel_execution_multiprocessing(func, entries, process_count=self.process_count)
        elif self.mode == RunnerMode.DASK_LOCAL:
            return parallel_execution_dask_local(func, entries, process_count=self.process_count, process_threads_count=self.thread_count)
        if self.mode == RunnerMode.DASK_JOB_QUEUE_CLUSTER:
            return parallel_execution_dask_cluster(func, entries, cluster=self.cluster)
        else:
            assert False
