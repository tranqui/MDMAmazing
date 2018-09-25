#!/usr/bin/env python3
"""Unit tests for distance calculations."""

import pytest

import numpy
from mdma.spatial.distance import pdist, pdist_histogram
from scipy.spatial.distance import pdist as npdist

@pytest.mark.parametrize("n", [10, 100])
@pytest.mark.parametrize("d", [2, 3, 4])
def test_pdist(n, d, eps=1e-14):
    """Test pairwise distances calculated correctly.

    Args:
        n: number of points
        d: dimensionality of space
        eps: tolerance in the floating point discrepancy to consider the test a success
    """
    x = numpy.random.random((n,d))
    assert numpy.linalg.norm(pdist(x) - npdist(x)) < eps

@pytest.mark.parametrize("n", [10, 100])
@pytest.mark.parametrize("d", [2, 3, 4])
def test_pdist_histogram(n, d, rmin=0, rmax=0.5, nbins=10):
    """Test histogram of pairwise distances calculated correctly.

    Args:
        n: number of points
        d: dimensionality of space
        rmin: start point included in histogram bins
        rmax: final point included in histogram bins
        nbins: number of bins across range
    """
    x = numpy.random.random((n,d))
    bins = numpy.linspace(rmin, rmax, nbins+1)
    histA = pdist_histogram(x, rmin, rmax, nbins).astype(numpy.float)
    histB = numpy.histogram(pdist(x), bins)[0]
    assert numpy.all(histA == histB)

# import timeit
# v = numpy.random.random((1000,3))
# print(timeit.timeit('npdist(v)', globals=globals(), number=5000))
# print(timeit.timeit('pdist(v)', globals=globals(), number=5000))
# #print(numpy.linalg.norm(npdist(v) - pdist(v)))
