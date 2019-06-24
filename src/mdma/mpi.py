#!/usr/bin/env python3
"""
Module giving helper functions for parallel processing with MPI.

The main functions are:
- synchronise: forces all MPI processes to wait, to ensure e.g. file operations are all in sync.
- parallel_map: runs a function on multiple data points taking advantage of the multiple MPI processes.
- pair_chunk: chunks data points into pairs for calculations involving multiple data points, e.g. correlations between snapshots in a trajectory.
"""

import itertools, numpy, math

from mpi4py import MPI
cpus = MPI.COMM_WORLD
num_cpus = cpus.Get_size()
rank = cpus.Get_rank()

def synchronise():
    """Makes all CPUs wait for the others to get to the same line.

    Useful for making CPUs wait for file operations."""
    cpus.alltoall(None)

def chunk_indices(size, n=num_cpus, i=rank):
    """Indices of data to divide a container into (approximately) equal sized chunks for processing by separate cpus.

    The chunks are made as equal as possible but this is only possible exactly when the number of
    elements can be evenly divided by the number of cpus.

    Args:
        size: length of data to break up
        n (int): number of cpus to divide the data amongst
        i (int): which cpu to chunk for
    Returns:
        indices: numpy array containing indices for this chunk

    Examples
    --------

    >>> for i in range(3): print(i, chunk_indices(10, 3, i))
    0 [0, 1, 2, 3]
    1 [4, 5, 6]
    2 [7, 8, 9]

    >>> for i in range(3): print(i, chunk_indices(12, 3, i))
    0 [0, 1, 2, 3]
    1 [4, 5, 6, 7]
    2 [8, 9, 10, 11]
    """
    k, m = size // n, size % n
    return numpy.arange(i*k + min(i, m), (i+1)*k + min(i+1, m))

def chunk(data, n=num_cpus, i=rank):
    """Break data up into (approximately) equal sized chunks for processing by separate cpus.

    The chunks are made as equal as possible but this is only possible exactly when the number of
    elements can be evenly divided by the number of cpus.

    Args:
        data (list or similar): data to divide
        n (int): number of cpus to divide the data amongst
        i (int): which cpu to chunk for
    Returns:
        chunk (iterator): chunk for the ith cpu

    Examples
    --------

    >>> for i in range(3): print(i, list(chunk(range(10), 3, i)))
    0 [0, 1, 2, 3]
    1 [4, 5, 6]
    2 [7, 8, 9]

    >>> for i in range(3): print(i, list(chunk(range(12), 3, i)))
    0 [0, 1, 2, 3]
    1 [4, 5, 6, 7]
    2 [8, 9, 10, 11]
    """
    k, m = len(data) // n, len(data) % n
    for l in range(i*k + min(i, m), (i+1)*k + min(i+1, m)): yield data[l]

def flatten(data):
    """Flatten a list of lists into a single list.

    Args:
        data: a list of lists
    Returns:
        list: flattened list

    Examples
    --------

    >>> flatten([[1, 2, 3], [4, 5, 6]])
    [1, 2, 3, 4, 5, 6]
    """
    if data is None: return
    else: return list(itertools.chain(*data))

def parallel_map(func, data, progress=None):
    """Extension of the standard map function in python, optimised for
    parallel processing.

    Args:
        func: function to run on each data point
        data: container of data points
        progress: optional progressbar widget to display progress
    Returns:
        list: results of function evaluated on each data point (returned to first cpu)
    """
    iterator = chunk(data)
    if progress and rank is 0: iterator = progress(list(iterator))
    return flatten(cpus.gather(list(map(func, iterator))))

def chunk_pairs(data, dt=1, max_samples=None, min_samples=1):
    """Return pairs between data points in a sorted container.

    A helper function for computing properties which depend on multiple
    data points, e.g. correlations between organised data.

    If a maximum number of sample pairs is specified then the chosen
    pairs will be spread out as much as possible: this is intended to
    minimise the correlations between samples. For example, if the data
    corresponds to a trajectory in time and we are computing a short-time
    correlation function, we may specify a maximum number of samples to
    speed up calculation over a long trajectory. However, if all of the
    short time samples are taken from the start of the trajectory, these
    may not be independent due to the presence longer-time correlations.
    By spreading the pairs out we maximise the independence of the measured
    correlations.

    Args:
        data: container of samples with some sort of linear order e.g. a trajectory composed of snapshots of equal time separation.
        dt: distance between pairs inside the data
        max_samples: maximum number of pairs to return (for details see description above)
    Returns:
        pairs (iterator): collection of pairs of data points for analysis

    Examples
    --------

    >>> list(chunk_pairs(range(5)))
    [(0, 1), (1, 2), (2, 3), (3, 4)]

    >>> list(chunk_pairs(range(5), 2))
    [(0, 2), (1, 3), (2, 4)]

    >>> list(chunk_pairs(range(5), 2, min_samples=4))
    []

    >>> list(chunk_pairs(range(10), 2, 3))
    [(0, 2), (3, 5), (6, 8)]
    """
    usable_length = len(data)-dt
    if max_samples is None: use_samples = usable_length
    else: use_samples = min(usable_length, max_samples)

    if use_samples >= min_samples:
        skip = max((len(data)+1-dt)//use_samples, 1)
        for i in range(0, skip*use_samples, skip):
            yield data[i], data[i+dt]

def takespread(sequence, num):
    """Evenly sample data points from a sequence.

    Args:
        sequence: data points to sample
        num: number of points to sample from sequence.
    Returns:
        points (iterator): sampled data points

    Examples
    --------

    >>> list(takespread(range(10),2))
    [0, 5]

    >>> list(takespread(range(10),3))
    [0, 4, 7]
    """
    for i in range(num):
        yield sequence[int(math.ceil(i * len(sequence) / num))]

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
