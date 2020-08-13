#!/usr/bin/env python3
"""
Module giving helper functions for evaluating correlation functions.

The main functions are:
- pair_chunk: chunks data points into pairs for calculations involving multiple data points, e.g. correlations between snapshots in a trajectory.
"""

import numpy

def chunk_indices(size, n, i):
    """ Indices of data to divide a container into (approximately) equal sized chunks for separate processing.

    The chunks are made as equal as possible but this is only possible exactly when the number of
    elements can be evenly divided by n.

    Args:
        size (int): length of data to break up
        n (int): number of chunks to divide the data into
        i (int): which entity to chunk for
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

def chunk(data, n, i):
    """ Break data up into (approximately) equal sized chunks for separate processing.

    The chunks are made as equal as possible but this is only possible exactly when the number of
    elements can be evenly divided by n.

    Args:
        data (container): data to divide
        n (int): number of chunks to divide the data into
        i (int): which entity to chunk for
    Returns:
        chunk (iterator): iterates through the different chunks

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

def chunk_pairs(data, dt=1, max_samples=None, min_samples=1):
    """ Return pairs between data points in a sorted container.

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
        data (container): samples with some sort of linear order e.g. a trajectory composed of snapshots of equal time separation.
        dt (int): distance between pairs inside the data
        max_samples (int): maximum number of pairs to return (for details see description above)
        min_samples (int): minimum number of pairs to return or else None will be returned (default=1)
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

def flatten(data):
    """ Flatten a list of lists into a single list.

    Args:
        data (container): a list of lists
    Returns:
        list: flattened list

    Examples
    --------

    >>> flatten([[1, 2, 3], [4, 5, 6]])
    [1, 2, 3, 4, 5, 6]
    """
    if data is None: return
    else: return list(itertools.chain(*data))

def takespread(sequence, num):
    """ Evenly sample data points from a sequence.

    Args:
        sequence (container): data points to sample
        num (int): number of points to sample from sequence.
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

def one_point_time_average(trajectory, func, max_samples=None, **func_kwargs):
    """ Average a one-point function over a trajectory (e.g. an observable or a spatial correlation function).

    Args:
        trajectory (container): container of snapshots forming the trajectory.
        func (function): function to evaluate at each point in time (taking a single snapshots as its argument).
        max_samples (int): maximum snapshots to sample.
        **func_kwargs: Additional keyword arguments are passed to the func call.
    Returns:
        The time averaged observable.
    """
    if max_samples is None: max_samples = len(trajectory)
    return numpy.average([func(snap, **func_kwargs) for snap in trajectory])

def two_point_time_average(trajectory, func, max_samples=None, initial_value=1, **func_kwargs):
    """ Average a two-point time correlation function over a trajectory.

    Args:
        trajectory (container): ordered container of snapshots assumed to be arranged have equal time separation between adjacent snapshots (linear time).
        func (function): correlation function to evaluate (taking two snapshots as its arguments), some general function :math:`G(t' - t)`.
        max_samples (int): maximum number of pairs to sample over.
        initial_value (scalar): initial value for correlation function :math:`G(0)` (default=1).
        **func_kwargs: Additional keyword arguments are passed to the func call.
    Returns:
        The time-averaged correlation function :math:`G(t' - t)`.
        Each element corresponds to a distinct :math:`\delta t = t' - t` in steps of the time separation between adjacent snapshots.
    """
    G = numpy.zeros(len(trajectory))
    G[0] = initial_value

    for dt in range(1, len(trajectory)):
        count = 0
        for snap1, snap2 in chunk_pairs(trajectory, dt, max_samples=20):
            G[dt] += func(snap1, snap2, **func_kwargs)
            count += 1
        G[dt] /= count

    return G
