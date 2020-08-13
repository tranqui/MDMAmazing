#!/usr/bin/env python3
"""
Module giving helper functions for parallel processing with MPI.

The main functions are:
- synchronise: forces all MPI processes to wait, to ensure e.g. file operations are all in sync.
- parallel_map: runs a function on multiple data points taking advantage of the multiple MPI processes.
"""

import itertools, numpy, math

from .correlation import chunk_indices, chunk, chunk_pairs, flatten, takespread

from mpi4py import MPI
cpus = MPI.COMM_WORLD
num_cpus = cpus.Get_size()
rank = cpus.Get_rank()

def synchronise():
    """ Makes all CPUs wait for the others to get to the same line.

    Useful for making CPUs wait for file operations."""
    cpus.alltoall(None)

def parallel_map(func, data, progress=None):
    """ Extension of the standard map function in python, optimised for
    parallel processing.

    Args:
        func: function to run on each data point
        data: container of data points
        progress: optional progressbar widget to display progress
    Returns:
        list: results of function evaluated on each data point (returned to first cpu)
    """
    iterator = chunk(data, num_cpus, rank)
    if progress and rank is 0: iterator = progress(list(iterator))
    return flatten(cpus.gather(list(map(func, iterator))))

if __name__ == '__main__':
    data = numpy.arange(20)
    def square(x): return x**2
    out = parallel_map(square, data)

    if rank is 0:
        print('A trivial example: compute square of first %d natural numbers using %d cpus\n' % (len(data), num_cpus))
        print()

    synchronise()
    print('result: cpu=%d return=%r' % (rank, out))
    synchronise()

    if rank is 0:
        print('\nTry running with a different number of cpus (via -n argument to mpirun)!')
