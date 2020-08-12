#!/usr/bin/env python3
"""
Module defining helper functions for distance calculations in periodic boundary
conditions.
"""

import numpy
from scipy.spatial.distance import cdist, squareform

def delta(x1, x2, box_dimensions):
    """ Displacement between two sets of points in a d-dimensional periodic space.

    Args:
        x1 (numpy.ndarray):
            An n by d array for first set of n points in a d-dimensional space.
        x2 (numpy.ndarray):
            An n by d array for second set of n points in a d-dimensional space.
        box_dimensions (numpy.ndarray):
            A d-dimensional array giving the box width in each dimension.
    Returns:
        delta (numpy.ndarray):
            x2 - x1 but with displacements wrapped inside the box boundaries
            (i.e. the nearest image convention).
    """
    n,d = x1.shape
    delta = x2 - x1
    for c in range(d):
        delta[:,c][delta[:,c] >  0.5*box_dimensions[c]] -= box_dimensions[c]
        delta[:,c][delta[:,c] < -0.5*box_dimensions[c]] += box_dimensions[c]
    return delta

def distance(x1, x2, box_dimensions, correct_drift=False):
    """ Distance between two sets of points in a d-dimensional periodic space.

    This functions analagously to numpy.spatial.cdist.

    Args:
        x1 (numpy.ndarray):
            An n by d array for first set of n points in a d-dimensional space.
        x2 (numpy.ndarray):
            An n by d array for second set of n points in a d-dimensional space.
        box_dimensions (numpy.ndarray):
            A d-dimensional array giving the box width in each dimension.
        correct_drift (bool):
            Whether to correct for drift in the centre of mass.
    Returns:
        delta (numpy.ndarray):
            The n distances i.e.
            :math:`\{|\\vec{x}_2^{(k)} - \\vec{x}_1^{(k)}|\}_{k=1}^n`,
            where displacements are wrapped inside the box boundaries
            (i.e. the nearest image convention).
    """
    dx = delta(x1,x2,box_dimensions)
    if correct_drift:
        dX = numpy.average(dx, axis=0)
        dx -= dX

    return numpy.linalg.norm(dx, axis=1)

def pdist(x, box_dimensions):
    """ Pairwise distances between points in a d-dimensional periodic space.

    Args:
        x (numpy.ndarray):
            An n by d array for the n points in a d-dimensional space.
        box_dimensions (numpy.ndarray):
            A d-dimensional array giving the box width in each dimension.
    Returns:
        delta (numpy.ndarray):
            The n(n-1)/2 distances between points.
    """
    n,d = x.shape
    delta = numpy.empty((n,n,d))
    for c in range(d):
        delta[:,:,c] = cdist(x[:,c].reshape(-1,1),x[:,c].reshape(-1,1))
        delta[:,:,c][delta[:,:,c] >  0.5*box_dimensions[c]] -= box_dimensions[c]
        delta[:,:,c][delta[:,:,c] < -0.5*box_dimensions[c]] += box_dimensions[c]
    return squareform(numpy.linalg.norm(delta, axis=2))

def self_overlap(x1, x2, box_dimensions, tol=0.3):
    """ Self overlap between two configurations of points.

    Overlap is defined as:

     .. math:: Q(\\vec{x}_1, \\vec{x}_2; \\delta) = \\frac{1}{N} \\sum_{k=1}^N \\Theta\\left( \\left| \\vec{x}_1^{(k)} - \\vec{x}_2^{(k)} \\right| - \\delta \\right)

    where :math:`\\Theta(\cdots)` is the `Heaviside step function
    <https://en.wikipedia.org/wiki/Heaviside_step_function>`_,
    :math:`\\vec{x}_{\\{1,2\\}}^{(k)}` indicates the kth particle position in the each system
    and :math:`\\delta` is a small parameter that determines whether particles are sufficiently
    close to be considered to overlap.
    :math:`\delta` is typically taken to be :math:`0.3\sigma` where :math:`\sigma` is
    the (effective) particle diameter.

    Args:
        x1 (numpy.ndarray):
            An n by d array for first set of n points in a d-dimensional space.
        x2 (numpy.ndarray):
            An n by d array for second set of n points in a d-dimensional space.
        box_dimensions (numpy.ndarray):
            A d-dimensional array giving the box width in each dimension.
        tol (bool):
            Tolerance to count an overlap (:math:`\delta` in the above equation).
    Returns:
        overlap (scalar):
            The self-overlap (:math:`Q` in the above equation).
    """
    return numpy.average(distance(x1,x2,box_dimensions,True) < tol)

def self_intermediate_scattering_function(x1, x2, box_dimensions, q=2*numpy.pi):
    """ Self intermediate scattering function between two configurations of points.

    The correlation function is defined as the Fourier transform of the self part of the
    `van Hove function <https://en.wikipedia.org/wiki/Dynamic_structure_factor#The_van_Hove_Function>`_:

    .. math:: F(\\vec{x}_1, \\vec{x}_2; \\vec{q}) = \\frac{1}{N} \\left\\langle \sum_{k=1}^N \exp{\\left(i \\vec{q} \cdot \\left( \\vec{x}_1^{(k)} - \\vec{x}_2^{(k)} \\right) \\right)} \\right\\rangle

    where :math:`\mathbf{q}` is the wave vector of the Fourier transform.
    This implementation assumes isotropy, so the exponential reduces to a
    `sinc <https://en.wikipedia.org/wiki/Sinc_function>`_ function
    and only the magnitude of the wavevector matters so :math:`q` is a scalar.
    :math:`|\\vec{q}|` is typically taken to be :math:`2\pi / \sigma`.

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
    return numpy.average(numpy.sinc((q/numpy.pi)*distance(x1,x2,box_dimensions,True)))
