#!/usr/bin/env python3
"""
Module for reading and writing snapshots from and to Protein Data Bank (.pdb)
file formats. The main class is PDBSnapshot, but some additional functions are
defined to provide a simplified interface to this class.

The module defines:
  - PDBSnapshot: the class the defining the file interface to this file format
  - read: shorthand for PDBSnapshot.read_single
  - write: create a snapshot from coordinates to write to disk
"""

import sys, io, numpy, pandas
from .snapshot import stream_safe_open, NoSnapshotError, Snapshot

__author__ = "Joshua Robinson and Christopher Brasnett"
__maintainer__ = "Joshua Robinson"

class PDBSnapshot(Snapshot):
    """Snapshot of a system of particles in Protein Data Bank (.pdb) file format.

    Interface defined in parent class Snapshot. Further documentation can be found there.
    """

    def read(self, path_or_file):
        """Read a snapshot from a file, overwriting any existing data.

        Args:
            path_or_file: file stream or path to read snapshot from
        Raises:
            NoSnapshotError: if could not read from file
            RuntimeException: if did not recognise file format
        """
        self.x = []
        self.species = []
        self.molecules_by_atom = []

        with stream_safe_open(path_or_file) as f:
            for line in f.readlines():
                line = line.split()
                tag = line[0]

                # Box comes from lattice parameter of the crystal unit cell.
                if tag == 'CRYST1':
                    self.box = numpy.array([list(map(float, line[1:4]))]).T
                    self.box = numpy.hstack((numpy.zeros((3,1)),
                                             self.box.reshape(3,1)))

                # Atom information.
                elif tag == 'ATOM':
                    # Chemical species information.
                    self.species += [line[2]]
                    self.molecules_by_atom += [int(line[4])]
                    # Atomic coordinates.
                    self.x += [list(map(float, line[5:8]))]

        self.x = numpy.array(self.x)
        self.species = numpy.array(self.species)
        self.molecules_by_atom = numpy.array(self.molecules_by_atom)

        self.molecule_labels = numpy.unique(self.molecules_by_atom)
        self.atoms_by_molecule = {m: [] for m in self.molecule_labels}
        for i,molecule in enumerate(self.molecules_by_atom):
            self.atoms_by_molecule[molecule] += [i]

    @property
    def molecule_sizes(self):
        """Sizes of each molecule in the snapshot."""
        return numpy.array([len(atoms) for molecule,atoms in self.atoms_by_molecule.items()])

    def molecule(self, molecule_number):
        """Summarise all atoms within a molecule."""
        atoms = self.atoms_by_molecule[molecule_number]
        x, species = self.x[atoms], self.species[atoms]

        summary = pandas.DataFrame({'molecule': molecule_number,
                                    'species': species,
                                    'x': x[:,0], 'y': x[:,1], 'z': x[:,2]})

        return summary

    @property
    def __str__(self):
        """String representation of the snapshot in Protein Data Bank (.pdb) format."""
        raise NotImplementedError
        f = io.StringIO()
        #f.write('?\n')
        return f.getvalue()

def read(*args, **kwargs):
    """Read a single snapshot from the disk."""
    return PDBSnapshot.read_single(*args, **kwargs)

def write(x, box, out=sys.stdout, species=None, **kwargs):
    """Write a single configuration to the disk."""
    snapshot = PDBSnapshot(x, box=box, species=species, **kwargs)
    snapshot.write(out)
