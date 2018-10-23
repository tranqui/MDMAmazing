Examples
########

Creating and running a LAMMPS simulation
========================================

Initialising and running the simulations
----------------------------------------

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

Processing the results
----------------------

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

Creating and running a DynamO simulation
========================================

Initialising and running the simulations
----------------------------------------

Coming soon.

Processing the results
----------------------

Suppose we have simulated a hard sphere system and produced a configuration file `config.end.xml` in the current directory. To read this file and convert it into the more flexible LAMMPS format we can try::

  from mdma import dynamo, atom
  snap = dynamo.read('config.end.xml')
  with open('dump.atom','w') as f:
      atom.write(snap.x, snap.box, f, species=snap.species)

The file `dump.atom` will be created, ready for analysis or visualisation (with e.g. `ovito <https://ovito.org/>`_).

Computing intermediate scattering functions
===========================================

Coming soon.

Parallel analysis with MPI
==========================

The following program gives an example of how to chunk data for speed up with mpi:

.. code-block:: python
   :caption: example.py
   :name: mpi-example

   from mdma import mpi
   from glob import glob

   files_to_process = glob('*.txt')
   print(mpi.rank, list(mpi.chunk(files_to_process)))

   def analysis_function(path):
       with open(path) as f:
           data = int(f.read())
       print(mpi.rank, 'file %s contains:' % path, data)

       return data

   data = mpi.parallel_map(analysis_function, files_to_process)
   print(mpi.rank, 'return:', data)

Suppose we have created 8 files each containing a number, which we can create on the Linux/Mac command line via::

  for i in $(seq 8); do echo $i > $i.txt; done

Then, running the program above with 3 cores produces the following output::

  >>> mpirun -n 3 python3 example.py
  0 ['6.txt', '1.txt', '8.txt']
  1 ['3.txt', '4.txt', '7.txt']
  2 ['5.txt', '2.txt']
  0 file 6.txt contains: 6
  1 file 3.txt contains: 3
  2 file 5.txt contains: 5
  0 file 1.txt contains: 1
  1 file 4.txt contains: 4
  0 file 8.txt contains: 8
  2 file 2.txt contains: 2
  1 file 7.txt contains: 7
  2 return: None
  1 return: None
  0 return: [6, 1, 8, 3, 4, 7, 5, 2]

Note that running the program without the :code:`mpirun` command will use normal serial analysis::

  >>> python3 example.py
  0 ['1.txt', '2.txt', '3.txt', '4.txt', '5.txt', '6.txt', '7.txt', '8.txt']
  0 file 1.txt contains: 1
  0 file 2.txt contains: 2
  0 file 3.txt contains: 3
  0 file 4.txt contains: 4
  0 file 5.txt contains: 5
  0 file 6.txt contains: 6
  0 file 7.txt contains: 7
  0 file 8.txt contains: 8
  0 return: [1, 2, 3, 4, 5, 6, 7, 8]
