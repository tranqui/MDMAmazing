Examples
########

Creating and running a LAMMPS simulation
========================================

In this section we will create a LAMMPS simulation of the Kob-Anderson binary Lennard-Jones potential::

  from mdma.lammps import potentials, wrapper

  # Define the system parameters to simulate.
  natoms, temperature, density = 1024, 1, 1.2
  system = potentials.KobAnderson(natoms, temperature, density=density)

  # Create the simulation.
  sim = wrapper.LammpsExecutable()
  sim.initialise_system(system)

To run the simulation for 100 timesteps we perform::

  sim.run(100)

To write a snapshot to a file we can use the atom module::

  atom_columns = ['id', 'type', 'x', 'y', 'z']
  sim.write_dump('all', 'custom', 'dump.atom', *atom_columns, 'modify', 'pbc', 'yes', 'sort', 'id')

This will create a file in the current directory named `dump.atom` storing the snapshot in LAMMPS'
atom format.

Sometimes it is desirable to save the coordinates in rescaled
coordinates which are bounded between 0 and 1, in which case we can do::

  atom_columns = ['id', 'type', 'xs', 'ys', 'zs']
  sim.write_dump('all', 'custom', 'dump.atom', *atom_columns, 'modify', 'pbc', 'yes', 'sort', 'id')

Alternatively, we can use the atom module to write the coordinates. This has the advantage of
being more transparent and flexible due to being written in python, at the cost of slower
performance::

  from mdma import atom
  with open('dump.atom','w') as f:
      atom.write(sim.coordinates, sim.box, f, species=sim.species)

Created atom files can be inspected with visualisation software (e.g. `ovito <https://ovito.org/>`_).
To read a snapshot previously stored as an atom file we can do::

  from mdma import atom
  snap = atom.read('dump.atom')

Comparison of this snapshot with the simulation will confirm that the snapshot is identical (up to rounding errors)::

  import numpy
  assert numpy.allclose(snap.box, sim.box)
  assert numpy.allclose(snap.box_dimensions, sim.box_dimensions)
  assert numpy.allclose(snap.x, sim.coordinates)


Computing intermediate scattering functions
===========================================

Coming soon.
