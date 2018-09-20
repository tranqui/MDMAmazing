#!/usr/bin/env python3
"""To write documentation..."""


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

    for i in range(3): print(i, chunk_indices(10, 3, i))
    >>> 0 [0, 1, 2, 3]
    >>> 1 [4, 5, 6]
    >>> 2 [7, 8, 9]

    for i in range(3): print(i, chunk_indices(12, 3, i))
    >>> 0 [0, 1, 2, 3]
    >>> 1 [4, 5, 6, 7]
    >>> 2 [8, 9, 10, 11]

    Args:
        size: length of data to break up
        n (int): number of cpus to divide the data amongst
        i (int): which cpu to chunk for
    Returns:
        indices: numpy array containing indices for this chunk
    """
    k, m = size // n, size % n
    return numpy.arange(i*k + min(i, m), (i+1)*k + min(i+1, m))

def chunk(data, n=num_cpus, i=rank):
    """Break data up into (approximately) equal sized chunks for processing by separate cpus.

    The chunks are made as equal as possible but this is only possible exactly when the number of
    elements can be evenly divided by the number of cpus.

    for i in range(3): print(i, list(chunk(range(10), 3, i)))
    >>> 0 [0, 1, 2, 3]
    >>> 1 [4, 5, 6]
    >>> 2 [7, 8, 9]

    for i in range(3): print(i, list(chunk(range(12), 3, i)))
    >>> 0 [0, 1, 2, 3]
    >>> 1 [4, 5, 6, 7]
    >>> 2 [8, 9, 10, 11]

    Args:
        data (list or similar): data to divide
        n (int): number of cpus to divide the data amongst
        i (int): which cpu to chunk for
    Returns:
        chunk (iterator): chunk for the ith cpu
    """
    k, m = len(data) // n, len(data) % n
    for l in range(i*k + min(i, m), (i+1)*k + min(i+1, m)): yield data[l]

def flatten(data):
    """Flatten a list of lists into a single list."""
    if data is None: return
    else: return list(itertools.chain(*data))

def parallel_map(func, data, progress=None):
    iterator = chunk(data)
    if progress and rank is 0: iterator = progress(list(iterator))
    return flatten(cpus.gather(list(map(func, iterator))))

def pairs(trajectory, dt=1, max_samples=None, min_samples=1):
    usable_length = len(trajectory)-dt
    if max_samples is None: use_samples = usable_length
    else: use_samples = min(usable_length, max_samples)

    if use_samples >= min_samples:
        skip = max((len(trajectory)+1-dt)/use_samples, 1)
        for i in range(0, skip*use_samples, skip):
            yield trajectory[i], trajectory[i+dt]

def takespread(sequence, num):
    length = float(len(sequence))
    for i in range(num): yield sequence[int(math.ceil(i * length / num))]

if __name__ == '__main__':
    data = numpy.arange(20)
    def square(x): return x**2
    out = parallel_map(square, data)
    print('final:', rank, out)
