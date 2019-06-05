#!/usr/bin/env python3
"""
Module for reading and writing snapshots from and to LAMMPS (.atom) file formats.
The main class is AtomSnapshot, but some additional functions are defined to provide
a simplified interface to this class.

The module defines:
  - AtomSnapshot: the class the defining the file interface to this file format
  - read: shorthand for AtomSnapshot.read_single
  - read_trajectory: shorthand for AtomSnapshot.read_trajectory
  - write: create a snapshot from coordinates to write to disk
  - write_trajectory: shorthand to write a series of coordinates to disk
"""

import sys, io, numpy, pandas
from .snapshot import stream_safe_open, NoSnapshotError, Snapshot

class AtomSnapshot(Snapshot):
    """Snapshot of a system of particles in LAMMPS (.atom) file format.

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
        with stream_safe_open(path_or_file) as f:
            self.x = numpy.empty((0,0))
            self.time = self.box = None

            while True:
                item = f.readline().split()
                if not item: raise NoSnapshotError
                assert item[0] == 'ITEM:'

                # Timestep within a trajectory.
                if item[1] == 'TIMESTEP':
                    line = f.readline()
                    try: self.time = int(line)
                    except ValueError: self.time = float(line)

                # Number of atoms in the header
                elif ' '.join(item[1:4]) == 'NUMBER OF ATOMS':
                    n = int(f.readline())
                    self.x = numpy.empty((n,self.d))

                # Item containing the bounding box.
                elif ' '.join(item[1:3]) == 'BOX BOUNDS':
                    d = len(item[3:])
                    self.x = numpy.empty((self.n,d))
                    boundaries = [f.readline().split() for _ in range(d)]

                    try:
                        self.box = numpy.zeros((d,2), dtype=numpy.longdouble)

                        for c,boundary in enumerate(boundaries):
                            self.box[c][:] = [float(b) for b in boundary]
                    except:
                        self.box = numpy.zeros((d,3), dtype=numpy.longdouble)

                        for c,boundary in enumerate(boundaries):
                            self.box[c][:] = [float(b) for b in boundary]

                # Main table contains the per-atom data. Should come at the end.
                elif item[1] == 'ATOMS':
                    assert self.n > 0
                    assert self.d >= 1 and self.d <= 3
                    assert self.box is not None

                    headings = item[2:]
                    assert 'id' in headings
                    assert 'x' or 'xs' in headings

                    c = io.StringIO()
                    for i in range(n): c.write(f.readline())
                    c.seek(0)
                    table = pandas.read_table(c, index_col=0, sep='\s+', names=headings, nrows=n)
                    #try: table = table.sort_values('id')
                    #except: table = table.sort('id')

                    if 'xs' in headings:
                        cols = ['xs','ys','zs'][:self.d]
                        self.x = table[cols].values.copy('c').astype(numpy.longdouble)
                        for c in range(d): self.x[:,c] *= self.box_dimensions[c]
                    else:
                        cols = ['x','y','z'][:self.d]
                        self.x = table[cols].values.copy('c').astype(numpy.longdouble)

                    self.species = numpy.array(table['type'])
                    return

                else: raise RuntimeError('unknown header: %s' % item)

    def __str__(self):
        """String representation of the snapshot in LAMMPS (.atom) format."""
        f = io.StringIO()
        f.write('ITEM: TIMESTEP\n%r\n' % self.time)
        f.write('ITEM: NUMBER OF ATOMS\n%r\n' % self.n)
        f.write('ITEM: BOX BOUNDS')
        for _ in self.box: f.write(' pp')
        f.write('\n')

        # Full left boundary and right boundary box conditions are specified.
        try:
            for c in range(self.d):
                for b in self.box[c]:
                    f.write('%.8f' % b)
                f.write('\n')
        # Only dimensions of box are specified
        except TypeError:
            for c in range(self.d):
                f.write('0 %.8f\n' % self.box[c])

        f.write('ITEM: ATOMS id type x y z')
        for i,(name,x) in enumerate(zip(self.species,self.x)):
            f.write('\n')
            f.write('%r %s' % (i,name))
            for coord in x: f.write(' %.8f' % coord)
        return f.getvalue()

def read(*args, **kwargs):
    """Read a single snapshot from the disk."""
    return AtomSnapshot.read_single(*args, **kwargs)
def read_trajectory(*args, **kwargs):
    """Read a trajectory (i.e. multiple snapshots) from the disk."""
    return AtomSnapshot.read_trajectory(*args, **kwargs)

def write(x, box, out=sys.stdout, species=None, **kwargs):
    """Write a single configuration to the disk."""
    snapshot = AtomSnapshot(x, box=box, species=species, **kwargs)
    snapshot.write(out)
def write_trajectory(trajectory, out=sys.stdout, **kwargs):
    """Write a configuration container to the disk."""
    for x in trajectory: write(x, out, **kwargs)
