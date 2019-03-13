#!/usr/bin/env python3
"""
Module defining helper functions for distance calculations in periodic boundary
conditions.
"""

import numpy
from scipy.spatial.distance import cdist, squareform

def delta(x1, x2, box):
    """ Displacement between two sets of points in a d-dimensional periodic space.

    Args:
        x1 (numpy.ndarray):
            An n by d array for first set of n points in a d-dimensional space.
        x2 (numpy.ndarray):
            An n by d array for second set of n points in a d-dimensional space.
        box (numpy.ndarray):
            A d-dimensional array giving the box width in each dimension.
    Returns:
        delta (numpy.ndarray):
            x2 - x1 but with displacements wrapped inside the box boundaries
            (i.e. the nearest image convention).
    """
    n,d = x1.shape
    delta = x2 - x1
    for c in range(d):
        delta[:,c][delta[:,c] >  0.5*box[c]] -= box[c]
        delta[:,c][delta[:,c] < -0.5*box[c]] += box[c]
    return delta

def distance(x1, x2, box, correct_drift=False):
    """ Distance between two sets of points in a d-dimensional periodic space.

    This functions analagously to numpy.spatial.cdist.

    Args:
        x1 (numpy.ndarray):
            An n by d array for first set of n points in a d-dimensional space.
        x2 (numpy.ndarray):
            An n by d array for second set of n points in a d-dimensional space.
        box (numpy.ndarray):
            A d-dimensional array giving the box width in each dimension.
        correct_drift (bool):
            Whether to correct for drift in the centre of mass.
    Returns:
        delta (numpy.ndarray):
            The n distances i.e.
            :math:`\{|\mathbf{x}_2^{(i)} - \mathbf{x}_1^{(i)}|\}_{i=1}^n`,
            where displacements are wrapped inside the box boundaries
            (i.e. the nearest image convention).
    """
    dx = delta(x1,x2,box)
    if correct_drift:
        dX = numpy.average(dx, axis=0)
        dx -= dX

    return numpy.linalg.norm(dx, axis=1)

def pdist(x, box):
    """ Pairwise distances between points in a d-dimensional periodic space.

    Args:
        x (numpy.ndarray):
            An n by d array for the n points in a d-dimensional space.
        box (numpy.ndarray):
            A d-dimensional array giving the box width in each dimension.
    Returns:
        delta (numpy.ndarray):
            The n(n-1)/2 distances between points.
    """
    n,d = x.shape
    delta = numpy.empty((n,n,d))
    for c in range(d):
        delta[:,:,c] = cdist(x[:,c].reshape(-1,1),x[:,c].reshape(-1,1))
        delta[:,:,c][delta[:,:,c] >  0.5*box[c]] -= box[c]
        delta[:,:,c][delta[:,:,c] < -0.5*box[c]] += box[c]
    return squareform(numpy.linalg.norm(delta, axis=2))

def self_overlap(x1, x2, box, tol=0.3):
    """ Self overlap between two configurations of points.

    Overlap is defined as:

        :math:`Q(t_1, t_2) = \\frac{1}{N} \sum_{i=1}^N \\Theta\Big(|\mathbf{x}_i(t_1) - \mathbf{x}_i(t_2)| - a\Big)`

    where :math:`\\Theta(\cdots)` is the
    `Heaviside theta function <https://en.wikipedia.org/wiki/Heaviside_step_function>`_
    and :math:`a` is the tolerance to count the pair as overlapping.

    Args:
        x1 (numpy.ndarray):
            An n by d array for first set of n points in a d-dimensional space.
        x2 (numpy.ndarray):
            An n by d array for second set of n points in a d-dimensional space.
        box (numpy.ndarray):
            A d-dimensional array giving the box width in each dimension.
        tol (bool):
            Tolerance to count an overlap (:math:`a` in the above equation).
    Returns:
        overlap (scalar):
            The self-overlap (:math:`Q` in the above equation).
    """
    return numpy.average(distance(x1,x2,box,True) < tol)

def self_intermediate_scattering_function(x1, x2, box, q=2*numpy.pi):
    """ Self intermediate scattering function between two configurations of points.

    The correlation function is defined as the Fourier transform of the
    self part of the
    `van Hove function <https://en.wikipedia.org/wiki/Dynamic_structure_factor#The_van_Hove_Function>`_
    (also known as the dynamic structure factor):

        :math:`F(t_1, t_2; \mathbf{q}) = \\frac{1}{N} \Big\langle \sum_{i=1}^N \exp{(i \mathbf{q} \cdot (\mathbf{x}_i(t_1) - \mathbf{x}_i(t_2)))} \Big\\rangle`

    where :math:`\mathbf{q}` is the wave vector of the Fourier transform.
    This implementation assumes isotropy, so the exponential reduces to a
    `sinc <https://en.wikipedia.org/wiki/Sinc_function>`_ function
    and only the magnitude of the wavevector matters so :math:`q` is a scalar.

    Args:
        x1 (numpy.ndarray):
            An n by d array for first set of n points in a d-dimensional space.
        x2 (numpy.ndarray):
            An n by d array for second set of n points in a d-dimensional space.
        box (numpy.ndarray):
            A d-dimensional array giving the box width in each dimension.
        q (scalar):
            Magnitude of the wavevector, typically taken as :math:`2\pi / \sigma`
            where :math:`\sigma` is the typical particle size.
    Returns:
        F (scalar):
            The self-intermediate scattering function
            (:math:`F` in the above equation).
    """
    return numpy.average(numpy.sinc((q/numpy.pi)*distance(x1,x2,box,True)))
