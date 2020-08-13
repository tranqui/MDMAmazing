Examples
########


Reading and writing snapshots
=============================


Snapshots are the fundamental units in :mod:`mdma` for analysis of molecular dynamics simulations.
We provides a number of snapshot classes corresponding to different common file formats used by various molecular dynamics software.
The main modules for this are:

* :mod:`mdma.xyz` provides :class:`mdma.xyz.XYZSnapshot` for reading XYZ format (.xyz extension)
* :mod:`mdma.atom` provides :class:`mdma.atom.AtomSnapshot` for reading LAMMPS dump files (.atom extension)
* :mod:`mdma.pdb` provides :class:`mdma.pdb.PDBSnapshot` for reading Protein Database format (.pdb extension)
* :mod:`mdma.dynamo` provides :class:`mdma.dynamo.xmlsnapshot.DynamoSnapshot` for reading DynamO format (.xml extension)

.. note:: If your desired file format is not in the above list, then you can create a new class by deriving from :class:`mdma.snapshot.Snapshot` following the templates above.
          If you are unsure how to do this, then feel free to submit a request that support for your file format be added to the `GitHub issue tracker <https://github.com/tranqui/MDMAmazing/issues>`_ or contact the `author <index.html#author>`_ directly with your request.

Each of the classes contain different information depending on the file format.

Reading snapshots from the disk
-------------------------------

Each of the above modules provides a :code:`read` function to read a single snapshot from a file (or an open filestream).

To read a single snapshot from the LAMMPS file `my_snapshot.atom` we run::

  from mdma import atom
  snap = atom.read('my_snapshot.atom')

Equivalently we could pass an open filestream as in::

  from mdma import atom
  with open('my_snapshot.atom') as f:
     snap = atom.read(f)

Whichever method used, the above examples will load the data into the variable :code:`snap` which can then be used. For example, we can obtain the number of particles, the time (in simulation time units) and the particle coordinates via::

  >>> snap.natoms
  10976
  >>> snap.time
  1234
  >>> snap.coordinates
  array([[52.631  , 54.3353 , 56.8596 ],
         [ 6.94414, 54.7188 , 50.713  ],
         [48.832  ,  1.72831,  2.68473],
         ...,
         [14.4879 , 44.3192 , 56.926  ],
         [21.4345 , 44.3192 , 31.5277 ],
         [10.4642 , 20.0305 , 31.3072 ]], dtype=float128)

Or we can use the following equivalent short-hand names for these variables::

  >>> snap.n
  10976
  >>> snap.t
  1234
  >>> snap.x
  array([[52.631  , 54.3353 , 56.8596 ],
         [ 6.94414, 54.7188 , 50.713  ],
         [48.832  ,  1.72831,  2.68473],
         ...,
         [14.4879 , 44.3192 , 56.926  ],
         [21.4345 , 44.3192 , 31.5277 ],
         [10.4642 , 20.0305 , 31.3072 ]], dtype=float128)

The simulation box and its length in each dimension can be obtained via::

  >>> snap.box
  array([[ 0.    , 58.2363],
         [ 0.    , 58.2363],
         [ 0.    , 58.2363]], dtype=float128)
  >>> snap.box_dimensions
  array([58.2363, 58.2363, 58.2363], dtype=float128)

Each row of the member variable :code:`snap.box` describes the left and right positions of a rectangular box in that dimension. Each entry of :code:`snap.box_dimensions` gives the length (right - left). This convention is used throughout the :mod:`mdma` package.

The chemical species are obtained with the :code:`snap.species` member variable, i.e.::

  >>> snap.species
  array([2, 1, 1, ..., 1, 2, 1])

.. note:: Which member variables are contained in the snapshot object will depend on the file format.
          Consult the documentation for the specific module to learn what data is contained in the particular snapshot.
          
          For example, the XYZ format does not state the simulation time so this will be absent from an :class:`mdma.xyz.XYZSnapshot` object, and may produce an error if you try to access it. The XYZ file format also does not define a box, so it will try to approximate a box from the coordinates which can lead to errors in the analysis. In general, it is better to use other file formats that contain more simulation information.

If you know the file contains many snapshots (e.g. in a trajectory), then you have to open it as a filestream so that the file will not be closed upon reading each snapshot. To read in 10 snapshots from a file try::

  with open('my_trajectory.atom') as f:
      for i in range(10):
          snap = atom.read(f)
          # your code processing the snapshot goes here

This example hard-codes the number of snapshots to read. If you do not know the number in advance, we provide a convenience function `read_trajectory` for each file type which returns a generator that can be looped over until the end of the file. For example::

  for snap in atom.read_trajectory('my_trajectory.atom'):
      # your code processing the snapshot goes here

The previous two examples are convenient for processing large trajectories because only a single snapshot is loaded into memory at one time. Sometimes it is necessary to load the entire trajectory into memory, which can be done as follows::

  trajectory = list(atom.read_trajectory('my_trajectory.atom')
  for snap in trajectory:
      # your code processing the snapshot goes here

.. warning:: Be careful when reading an entire trajectory into memory, as this can easily consume a large portion of available resources for large systems and/or long trajectories.

Another common situation is for the snapshots forming a trajectory to be stored in separate files.
In this example we assume that there are files in the current directory labelled `1.atom, 2.atom, 3.atom` etc that sequentially describe a complete trajectory.
We can obtain a list of all files in the current directory with the `.atom` extension using `glob <https://docs.python.org/3/library/glob.html>`_, but these will not necessarily be in the correct order (the snapshot `10.atom` would come before `2.atom`) so to obtain the correct order we make use of the `natsort <https://pypi.org/project/natsort/>`_ module::

  from glob import glob
  from natsort import natsorted
  paths_in_order = natsorted(glob('*.atom'))

Now we can load the trajectory into memory via::

  trajectory = [atom.read(path) for path in paths_in_order]

For other fileformats replace :code:`atom` in the above examples with any of the other modules listed at the start of this section. For example, try::

  from mdma import xyz
  trajectory = list(xyz.read_trajectory('my_trajectory.xyz'))

to read in a trajectory in XYZ format.

Printing snapshots to the console or writing them to the disk
-------------------------------------------------------------

To print the snapshot to the console in its native format we can use::

  print(snap)

which for simple applications can be combined with `BASH redirects <https://www.gnu.org/software/bash/manual/html_node/Redirections.html>`_ on Linux to output snapshot files.
For more control over where the snapshot is written you can use the method :meth:`mdma.snapshot.Snapshot.write`, e.g.::

  snap.write('my_snapshot.atom')

Or equivalently we can pass an open file handle::

  with open('my_snapshot.atom', 'w') as f:
      snap.write(f)

The latter example is particularly useful for writing entire trajectories to a single file, because we can chain calls to :meth:`mdma.snapshot.Snapshot.write`, i.e.::

  with open('my_trajectory.atom', 'w') as f:
      for snap in trajectory:
          snap.write(f)

If your data is not contained in a snapshot object (e.g. if you have the raw coordinates/box in numpy arrays) then you can use the functions :code:`write` or :code:`write_trajectory` inside the relevant snapshot module.
Refer to the documentation inside your module for how to use these, e.g. for :mod:`mdma.atom` refer to :func:`mdma.atom.write` and :meth:`mdma.atom.write_trajectory`.


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

To write a snapshot to a file we can use the native LAMMPS dump files::

  atom_columns = ['id', 'type', 'x', 'y', 'z']
  sim.write_dump('all', 'custom', 'dump.atom', *atom_columns, 'modify', 'pbc', 'yes', 'sort', 'id')

This will create a file in the current directory named `dump.atom` storing the snapshot in LAMMPS' atom format.

Sometimes it is desirable to save the coordinates in rescaled coordinates which are bounded between 0 and 1, in which case we can do::

  atom_columns = ['id', 'type', 'xs', 'ys', 'zs']
  sim.write_dump('all', 'custom', 'dump.atom', *atom_columns, 'modify', 'pbc', 'yes', 'sort', 'id')

Alternatively, we can use the :mod:`mdma.atom` module to write the coordinates (see also `above <#reading-and-writing-snapshots>`_ for more comprehensive overview of reading/writing snapshots).
This has the advantage of being more transparent and flexible due to being written in python, at the cost of slower performance::

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

Suppose we have simulated a hard sphere system and produced a configuration file `config.end.xml` in the current directory.
To read this file and convert it into the more flexible LAMMPS format we can try::

  from mdma import dynamo, atom
  snap = dynamo.read('config.end.xml')
  with open('dump.atom','w') as f:
      atom.write(snap.x, snap.box, f, species=snap.species)

The file `dump.atom` will be created, ready for analysis or visualisation (with e.g. `ovito <https://ovito.org/>`_).


Common two-point correlation functions
======================================

We are going to explore how to obtain some two-point functions commonly used to analyse molecular dynamics simulations, of the general form

.. math:: F = F(\vec{x}_1, \vec{x}_2)

where :math:`\vec{x}_{\{1,2\}}` are two sets of coordinates that could correspond to two different systems, or the same system at two different times.
We require that the number of particles :math:`N` (and dimensions) are the same for both systems.
We will run through examples of how to calculate some common correlation functions, then we will show how to average these correlation functions over a trajectory to obtain correlation functions for two-points in time.

Spatial correlation functions are defined within the submodule :mod:`mdma.spatial`.
Currently only simulations in periodic boxes are supported, so the only module there is :mod:`mdma.spatial.periodic`.
In all of the examples in subsequent sections this is assumed to have been imported via::

  from mdma.spatial import periodic

We assume that two snapshots are loaded called :code:`snap1` and :code:`snap2`, that correspond to the two systems above.
See `Reading and writing snapshots <#reading-and-writing-snapshots>`_ for examples showing how to read snapshots.

.. todo:: Generalise correlation functions for NPT simulations where the box size will fluctuate between two points in time.

Displacements
-------------

One of the simplest quantities for spatial correlations is the displacement between two sets of coordinates, which is crucial for calculating other more interesting quantities.
In the absence of periodic boundary conditions, displacements are extremely simple, i.e. we take the difference::

  dx = snap1.x - snap2.x

However, with periodic boundary conditions we have to take into account the wrapping at the boundaries.
To this we provide the following functions:

* :func:`mdma.spatial.periodic.delta`: calculates the displacement between particles with the same indices in the two systems
* :func:`mdma.spatial.periodic.distance`: calculates the :math:`N` distances between particles with the same indices in the two systems
* :func:`mdma.spatial.periodic.pdist`: calculates all :math:`N(N-1)/2` distances between all particles within a *single* system. This is the periodic equivalent to `scipy.spatial.distance.pdist <https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.pdist.html>`_.

The equivalent of the above example for periodic systems would be::

  dx = periodic.delta(snap1.x, snap2.x, snap1.box_dimensions)

.. note:: We assume the box dimensions are the same in both systems. We also make this assumption in all subsequent examples.
          If the box size differs then the correlation functions will produce erroneous results.

Self-overlap
------------

The self-overlap is defined as

.. math:: Q(\vec{x}_1, \vec{x}_2; \delta) = \frac{1}{N} \sum_{k=1}^N \Theta\left( \left| \vec{x}_1^{(k)} - \vec{x}_2^{(k)} \right| - \delta \right)

where :math:`\Theta(\cdots)` is the `Heaviside step function <https://en.wikipedia.org/wiki/Heaviside_step_function>`_, :math:`\vec{x}_{\{1,2\}}^{(k)}` indicates the kth particle position in the each system and :math:`\delta` is a small parameter that determines whether particles are sufficiently close to be considered to overlap.
:math:`\delta` is typically taken to be :math:`0.3\sigma` where :math:`\sigma` is the (effective) particle diameter.

To compute this quantity we provide :func:`mdma.spatial.periodic.self_overlap`, which can be used via::

  Q = periodic.self_overlap(snap1, snap2, tol=0.3)

This can be called equivalently with the raw coordinate data, i.e.::

  Q = periodic.self_overlap(snap1, snap2.x, snap1.box_dimensions, tol=0.3)

which is useful when you have data outside of a snapshot instance.
Refer to the documentation of :func:`mdma.spatial.periodic.self_overlap` for descriptions of the arguments.

.. todo:: Show how to calculate the overlap between two clusters (not periodic), which requires finding the optimal alignment.

Intermediate scattering function
--------------------------------

The self intermediate scattering function (ISF) is defined as the Fourier transform of the self part of the `van Hove function <https://en.wikipedia.org/wiki/Dynamic_structure_factor#The_van_Hove_Function>`_:

.. math:: F(\vec{x}_1, \vec{x}_2; \vec{q}) = \frac{1}{N} \sum_{k=1}^N \exp{\left(i \vec{q} \cdot \left( \vec{x}_1^{(k)} - \vec{x}_2^{(k)} \right) \right)}

:math:`|\vec{q}|` is typically taken to be :math:`2\pi / \sigma`.

.. note:: We assume isotropy and :math:`d=3` so :math:`\vec{q} \to |\vec{q}|`, and the exponential reduces to a `sinc <https://en.wikipedia.org/wiki/Sinc_function>`_ function.
.. todo:: Replace sinc implementation with a spherical Bessel function for arbitrary :math:`d`.

To compute this quantity we provide :func:`mdma.spatial.periodic.self_intermediate_scattering_function`, which can be used via::

  F = periodic.self_intermediate_scattering_function(snap1, snap2)

As before, this can be called equivalently with the raw coordinate data, i.e.::

  F = periodic.self_overlap(snap1, snap2.x, snap1.box_dimensions, tol=0.3)

which is useful when you have data outside of a snapshot instance.

.. warning:: The above examples will give erroneous results in general, because the self-ISF function takes a fourth argument :math:`|\vec{q}|` which we have ignored.
             By default this function sets :math:`|\vec{q}| = 2\pi` if this is not specified, which implicitly assumes the effective diameter :math:`\sigma = 1`.
             In general you must pass the wavenumber explicitly to get reasonable results.
             Refer to the documentation :func:`mdma.spatial.periodic.self_intermediate_scattering_function` for descriptions of the additional arguments.

Averaging temporal correlation functions over trajectories
----------------------------------------------------------

A common operation is to take some two-point correlation function, and find its average value in equilibrium (where time-translation invariance is recovered) i.e.\

.. math:: \lim_{t \to \infty} \left \langle G(t, t') \rangle = \langle G(\delta t = t' - t) \right \rangle

for some correlation function :math:`G(t, t')` and where :math:`\langle \cdots \rangle` indicates an ensemble average.
The ensemble average equals the long-time average in equilibrium, so we can evaluate this for a long trajectory by sampling over all the pairs of snapshots in a trajectory (although in practice a subset usually suffices).

.. note:: In the literature on supercooled liquids a trajectory is conventionally taken to be *long-enough* for this procedure if it contains several decays of the time-correlation functions, e.g. if the ISF decays :math:`\mathcal{O}(10)` times.
          In the latter example, by a "decay" we mean that the correlation function reaches :math:`F(\delta t) \le 1/e` from an initial value of :math:`F(\delta t = 0) = 1`, and the reference time is reset when this event occurs.

Assuming we have loaded a trajectory into the variable :code:`trajectory` (see `Reading and writing snapshots <#reading-and-writing-snapshots>`_ for examples showing how to do this).
We can obtain a quick estimate of what the correlation function looks like by only comparing with the first snapshot, i.e.::

  import numpy as np
  F = np.empty(len(trajectory))
  snap1 = trajectory[0]
  for dt in range(len(trajectory)):
      snap2 = trajectory[dt]
      F[dt] = periodic.self_intermediate_scattering_function(snap1, snap2)

Again, we have assumed that the box dimensions and number of particles are constant throughout the trajectory.
Plotting the variable :math:`F` at this point can give a rough idea of how it is varying.

The above example will typically feature a lot of noise because each value of :math:`\delta t` only contains a single sample.
In general, it is much better to perform the average via::

  import numpy as np
  F = np.zeros(len(trajectory))
  F[0] = 1
  for dt in range(1, len(trajectory)):
      for i in range(len(trajectory) - dt):
          j = i + dt
          snap1 = trajectory[i]
          snap2 = trajectory[j]
          F[dt] += periodic.self_intermediate_scattering_function(snap1, snap2)
      F[dt] /= len(trajectory) - dt

This code snippet evaluates :math:`F(\delta t)` over all :math:`m(m-1)/2` pairs of snapshots, where :math:`m` is the number of snapshots, and can be quite slow.
A better approach is to take a subset of reference times for each :math:`\delta t` so that this average becomes :math:`\mathcal{O}(m)` rather than :math:`\mathcal{O}(m^2)`.
When sampling a subset of reference points it is best to space them as far apart as much as possible, because adjacent snapshots (in time) will be highly correlated; to do this we use the function :func:`mdma.correlation.chunk_pairs` ::

  import numpy as np
  from mdma import correlation
  F = np.zeros(len(trajectory))
  F[0] = 1
  for dt in range(1, len(trajectory)):
      count = 0
      for snap1, snap2 in correlation.chunk_pairs(trajectory, dt, max_samples=25):
          F[dt] += periodic.self_intermediate_scattering_function(snap1, snap2)
          count += 1
      F[dt] /= count

The above example will take a maximum of 25 samples for each value of :math:`\delta t` and thus involves :math:`\sim 25m` operations (though slightly fewer than this because there will not be 25 samples for the largest values of :math:`\delta t`).
This example shows how you can start to build more sophisticated correlation functions, but if this is all you wish to do we provide :func:`mdma.correlation.two_point_time_average` with the same functionality::

  from mdma import correlation
  F = correlation.two_point_time_average(trajectory, periodic.self_intermediate_scattering_function, max_samples=25)

is equivalent to the previous example and much more succinct.
You can still pass additional arguments to the correlation function by name, e.g. if you want to set the wavenumber to :math:`2\pi / 5` (for effective particle diameters of :math:`\sigma = 5`) then you would change this to::

  from mdma import correlation
  F = correlation.two_point_time_average(trajectory, periodic.self_intermediate_scattering_function, max_samples=25, q=2*np.pi/5)

The following is a minimal working example that calculates and plots an ISF for a trajectory contained in the file `my_trajectory.atom` ::

  import numpy as np
  from mdma import atom, correlation
  from mdma.spatial import periodic

  # We have to load the entire trajectory into memory by casting to list(...)
  # because the two-point correlators need to reference multiple snapshots.
  trajectory = list(atom.read_trajectory('my_trajectory.atom'))

  # Extract the allowable values of dt.
  # (NB: these should be linear increments or the time-averaging in the next step will fail)
  t0 = trajectory[0].time
  dt = [snap.time - t0 for snap in trajectory]

  # Calculate the ISF with a time-average to reduce noise.
  sigma = 1 # set to the effective diameter for your system
  F = correlation.two_point_time_average(trajectory, periodic.self_intermediate_scattering_function, max_samples=25, q=2*np.pi/sigma)

  # Plot the function.
  import matplotlib.pyplot as plt
  plt.semilogx(dt, F)
  plt.xlabel('t')
  plt.ylabel('F(t)')
  plt.show()

This assumes you have `matplotlib <https://matplotlib.org>`_ installed for plotting the function.

.. note:: Replace :code:`self_intermediate_scattering_function` with the correlation function of your choice in the examples above.

.. todo:: Vectorise the above.


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

  >>> mpirun -n 3 python example.py
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

  >>> python example.py
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
