#!/usr/bin/env python3

import numpy
from scipy.spatial.distance import cdist, squareform

def delta(x1, x2, box):
    n,d = x1.shape
    delta = x2 - x1
    for c in range(d):
        delta[:,c][delta[:,c] >  0.5*box[c]] -= box[c]
        delta[:,c][delta[:,c] < -0.5*box[c]] += box[c]
    return delta

def distance(x1, x2, box, correct_drift=False):
    dx = delta(x1,x2,box)
    if correct_drift:
        dX = numpy.average(dx, axis=0)
        dx -= dX

    return numpy.linalg.norm(dx, axis=1)

def self_overlap(x1, x2, box, tol=0.3):
    return numpy.sum(distance(x1,x2,box,True) < tol) / len(x1)

def self_intermediate_scattering_function(x1, x2, box, q=2*numpy.pi):
    return numpy.sum(numpy.sinc((q/numpy.pi)*distance(x1,x2,box,True))) / len(x1)

def pdist(x, box):
    n,d = x.shape
    delta = numpy.empty((n,n,d))
    for c in range(d):
        delta[:,:,c] = cdist(x[:,c].reshape(-1,1),x[:,c].reshape(-1,1))
        delta[:,:,c][delta[:,:,c] >  0.5*box[c]] -= box[c]
        delta[:,:,c][delta[:,:,c] < -0.5*box[c]] += box[c]
    return squareform(numpy.linalg.norm(delta, axis=2))
