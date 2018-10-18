Examples
########

Creating and a LAMMPS simulation
================================

In this section we will create::

  natoms, temperature, density = 1024, 1, 1.2
  system = KobAnderson(natoms, temperature, density=density)

  sim = LammpsExecutable()
  sim.initialise_system(system)

To run the simulation for 100 timesteps we perform::

  sim.run(100)

To write a snapshot to a file we must use the atom module::

  from mdma import atom
  with open('dump.atom','w') as f:
      atom.write(sim.coordinates, sim.box, f, species=sim.species)

The created file can be inspected with visualisation software (e.g. `ovito <https://ovito.org/>`_).

Computing intermediate scattering functions
===========================================

Coming soon.
