#!/usr/bin/env python3
"""
Module defining a wrapper on top of the python interface to LAMMPS in order to
streamline the process of creating and running molecular dynamics simulations.
"""

import numpy, pandas
import os, ast, json, string, ctypes
from glob import glob
from natsort import natsorted

import lammps
from mdma import atom

def random_seed():
    """Generate a random seed for use by LAMMPS."""
    return numpy.random.randint(numpy.iinfo(ctypes.c_int).max)

# Default columns when dumping snapshots in .atom format.
atom_columns = ['id', 'type', 'x', 'y', 'z', 'fx', 'fy', 'fz'] #, 'vx', 'vy', 'vz']

class Ensemble:
    """Numerical identifiers for the possible ensembles.

    This is intended as a simple enumeration implementation.
    """
    NVT = 1 # canonical ensemble
    NPT = 2 # isobaric ensemble

class LammpsExecutable(lammps.PyLammps):
    """
    Executable for LAMMPS simulations.

    This is the main entry point for starting new simulations or resuming old ones.
    """

    def __init__(self, *args, **kwargs):
        """Constructor defaults to PyLammps base class.

        For most applications this should be called without arguments.

        See `PyLammps documentation <https://lammps.sandia.gov/doc/Howto_pylammps.html#creating-a-new-instance-of-pylammps>`_
        for information on possible arguments.
        """
        super().__init__(*args, **kwargs)
        self.initialised = False
        # exec_id keeps track of whether this is a new instance of a trajectory
        # or a continuation of a previous one (i.e. resumed from disk).
        self.exec_id = 0

    def initialise_system(self, system, initial_snapshot=None):
        """Initialise a new LAMMPS simulation from parameters.

        Args:
            system (mdma.lammps.PotentialBase):
                The parameters of the model system to simulate, as well as 
                simulation parameters e.g. the thermostat or neighbour detection.
            initial_snapshot (mdma.Snapshot or None):
                The starting configuration for the simulation.
                If None then one will be generated from a random seed.
        Raises:
            RuntimeError:
                When attempting to initialise an existing simulation.
            ValueError:
                When system contains a partial description of a simulation.
                This happens when:

                    * Neither a density, pressure or an initial configuration is given so a starting configuration cannot be generated.
                    * Pressure is specified but not an initial configuration.
                    * An unknown ensemble is specified (this should never happen as the ensemble is chosen from parsing the incoming data).
        """
        if self.initialised:
            raise RuntimeError('cannot overwrite a previous simulation!')
        self.initialised = True

        self.units(system.units)
        self.atom_style('atomic')
        self.atom_modify('map', 'array')
        self.boundary(' '.join(['p']*system.d))

        # Determine ensemble based on system parameters.
        try:
            pressure = system.pressure
            ensemble = Ensemble.NPT
        except AttributeError:
            try:
                final_density = system.density
                ensemble = Ensemble.NVT
            except AttributeError:
                if initial_snapshot is None:
                    raise ValueError('without an initial snapshot (one of) density or pressure must be set (for NVT or NPT simulations respectively)')
                else: ensemble = Ensemble.NVT

        # Create the box (for now only cubic boxes are implemented)
        initial_box = self.initial_box(system, initial_snapshot)
        self.region('box', 'block', *initial_box.reshape(-1))
        self.create_box(system.num_species, 'box')

        # Initialise coordinates and species.
        self.create_atoms(1, 'random', system.natoms, random_seed(), 'NULL')
        if initial_snapshot is None:
            self.composition = system.composition
        else:
            self.coordinates = initial_snapshot.x.astype(numpy.double)
            try: self.species = initial_snapshot.species
            except:
                species = [1+string.ascii_uppercase.index(s) for s in initial_snapshot.species]
                self.species = numpy.array(species)

        # Initialise interactions.
        system.initialise_potential(self)
        self.timestep(system.timestep)
        self.dt = system.timestep

        # Minimise in case there are (high energy) overlaps from the random generation.
        if initial_snapshot is None:
            energy_tol, force_tol, max_iters, max_evals = 1e-4, 1e-6, 100, 1000
            self.minimize(energy_tol, force_tol, max_iters, max_evals)
            self.reset_timestep(0)

        self.velocity('all', 'create', system.temperature, random_seed(), 'loop', 'geom')
        # set neighbour options here? should maybe include an option to customise
        self.run_style('verlet')

        if ensemble == Ensemble.NVT:
            self.fix(1, 'all', 'nvt',
                     'temp', system.temperature, system.temperature,
                     system.thermostat_damping_interval)

            if initial_snapshot is None:
                initial_box_dimensions = initial_box[:,1] - initial_box[:,0]
                initial_volume = numpy.prod(initial_box_dimensions)
                desired_volume = system.natoms / system.density

                rescale = (desired_volume / initial_volume)**(1./system.d)
                final_box = initial_box * rescale

                self.fix(2, 'all', 'deform', 1,
                        'x', 'final', final_box[0,0], final_box[0,1],
                        'y', 'final', final_box[1,0], final_box[1,1],
                        'z', 'final', final_box[2,0], final_box[2,1])
                self.run(system.initial_thermalisation_steps)
                self.unfix(2)
                self.run(system.initial_thermalisation_steps)
                self.reset_timestep(0)
                self.runs = []

        elif ensemble == Ensemble.NPT:
            self.fix(1, 'all', 'npt',
                     'temp', system.temperature, system.temperature,
                     system.thermostat_damping_interval,
                     'iso', system.pressure, system.pressure,
                     system.barostat_damping_interval)

            if initial_snapshot is None:
                raise ValueError('must give initial snapshot for NPT run!')

        else:
            raise ValueError('unknown ensemble encountered (ensemble=%d)!' % ensemble)
        self.ensemble = ensemble

        # Running for zero steps does not actually run the simulation, but it
        # forces lammps to compute thermodynamic properties for initial state which
        # can be retrieved (otherwise an error occurs).
        self.run(0)

    def initial_box(self, system, initial_snapshot, seed_density=1e-3):
        """Determine size of initial box for system.

        This is an internal helper function for initialise_system, so the typical
        end user should not need to touch this.

        There are three (valid) scenarios:

            1. An initial snapshot is specified but not a density: the snapshots box is returned.
            2. An initial snapshot is specified as well as a density: the snapshots box is rescaled to match the target density and returned.
            3. A snapshot is not specified but a density is: a cubic box of the target density is returned.

        Args:
            system (mdma.lammps.PotentialBase):
                The parameters of the model system to simulate, as well as 
                simulation parameters e.g. the density which will be used to set
                the box dimensions if specified.
            initial_snapshot (mdma.Snapshot or None):
                The starting configuration for the simulation.
                If None then that means one should be generated from a random seed,
                and a cubic box will be chosen for this purpose.
        Returns:
            initial_box (numpy.ndarray):
                The box boundaries in each dimension.
        """
        try:
            initial_box = initial_snapshot.box
            initial_box_dimensions = (initial_box[:,1] - initial_box[:,0])**(1./system.d)

            # If density is specified we need to rescale the box appropriately
            try:
                desired_volume = system.natoms / system.density
                snapshot_volume = numpy.prod(initial_box_dimensions)
                rescale = (desired_volume / snapshot_volume)**(1./system.d)

                initial_box *= rescale
                initial_box_dimensions *= rescale

            # Density was not specified by system: ignore (we *could* be in NPT ensemble so not necessarily an error).
            except AttributeError:
                pass

            return initial_box

        # Either no snapshot was given, or one was given without a box specified:
        # we will have to generate coordinates from random so we must start with a low
        # density box (which will later be shrunk to achieve the desired density).
        except AttributeError:
            unit_box = numpy.array([0, 1.]*system.d).reshape(system.d,2)
            initial_volume = system.natoms / seed_density
            initial_box = unit_box * initial_volume**(1./system.d)

            return initial_box

    def save(self, path, overwrite=True):
        """Save the simulation to the disk.

        If the simulation is a continuation of a previous trajectory, then this
        will only append to existing data.

        Args:
            path (str): Path to write to.
            overwrite (bool):
                Whether it is okay to overwrite an existing trajectory.
                Should be set to false the first time writing a new trajectory to
                disk to avoid deleting old data.
        """
        os.makedirs(path, exist_ok=overwrite)
        os.makedirs('%s/dumps' % path, exist_ok=overwrite)

        # run_path = '%s/%d.runs' % (path,self.exec_id)
        # history_path = '%s/%d.history' % (path,self.exec_id)
        # if not overwrite and (os.path.exists(run_path) or os.path.exists(history_path)):
        #     raise RuntimeError('saving to %s will overwrite previous simulation!' % path)

        with open('%s/%d.runs' % (path,self.exec_id), 'w') as f:
            for run in self.runs:
                thermo_data = run.thermo.__dict__
                json.dump(thermo_data, f)
                f.write('\n')

        self.write_restart('%s/dumps/*.restart' % path)
        self.write_dump('all', 'custom', '%s/dumps/*.atom' % path, *atom_columns, 'modify', 'pbc', 'yes', 'sort', 'id')

        history = '\n'.join(self._cmd_history)
        history = history.replace('%s/' % path, '')
        with open('%s/%d.history' % (path,self.exec_id), 'w') as f:
            f.write(history)

    def resume(self, path):
        """Resume a saved trajectory from the disk.

        Args:
            path (str): Location of the trajectory to load.
        Raises:
            RuntimeError:
               If attempting to load a trajectory over an existing simulation;
               this should only be called on new LammpsExecutable objects to
               prevent loss of data.
        """
        if self.initialised:
            raise RuntimeError('cannot overwrite a previous simulation!')
        self.initialised = True

        latest_restart = natsorted(glob('%s/dumps/*.restart' % path))[-1]
        self.read_restart(latest_restart)

        for run_path in natsorted(glob('%s/*.runs' % path)):
            with open(run_path) as f:
                for line in f.readlines():
                    dictionary = ast.literal_eval(line)
                    del dictionary['_name']
                    self.runs += [{'thermo': dictionary}]

        previous_logs = natsorted(glob('%s/*.history' % path))
        history = []
        for log in previous_logs: history += self.read_script(log)

        # Neighbour list and fix commands need to be reissued
        self.restore_fixes(history)
        for command in history:
            if 'pair' in command or 'neigh' in command:
                self.command(command)

        # Integration timestep dt needs to be restored
        timestep_commands = [command for command in history if command.split()[0] == 'timestep']
        self.dt = float(timestep_commands[-1].split()[1])
        self.timestep(self.dt)

        latest_log = previous_logs[-1]
        self.exec_id = latest_log.split('/')[-1]
        self.exec_id = self.exec_id.split('.history')[0]
        self.exec_id = 1 + int(self.exec_id)

    def restore_fixes(self, history):
        """Implementation detail when resuming from the disk: restores thermostat
        fixes from the command history.

        Args:
            history (list):
                Command history from the trajectory.
                Active fixes are determined and applied from this.
        """
        fix_history = [command for command in history if 'fix' in command.split()[0]]
        fix_labels = numpy.array([command.split()[1] for command in fix_history])

        unique_fixes = numpy.unique(fix_labels)
        for label in unique_fixes:
            fixes = numpy.where(fix_labels == label)[0]
            latest_command = fix_history[fixes[-1]]

            # Don't do anything if the last command with this label was an 'unfix'
            active = latest_command.split()[0] == 'fix'
            if active:
                self.command(latest_command)

                fix_type = latest_command.split()[3]
                if  fix_type == 'nvt':
                    self.ensemble = Ensemble.NVT
                elif fix_type == 'npt':
                    self.ensemble = Ensemble.NPT

    def read_script(self, path):
        """Implementation detail for resume: reads e.g. history files.

        Args:
            path (str): Location of the file to read.
        Returns:
            lines (list): List of the lines (strings) in the file.
        """
        with open(path) as f:
            return f.read().splitlines()

    @property
    def thermodynamics(self):
        """Return table of summary thermodynamic data from trajectory's history.

        Returns:
            A pandas.DataFrame object summarising the thermodynamics.
        """
        # Headings of selected data columns.
        headings = ['Step', 'Temp', 'Press', 'TotEng']
        if self.ensemble == Ensemble.NPT: headings += ['Volume']

        # Populate the table with data.
        table = {key: [] for key in headings}
        for run in self.runs:
            if type(run) is not dict:
                run_data = run.thermo.__dict__
            else: run_data = run['thermo']

            for key in headings:
                val = run_data[key]
                if len(val) > 1: val = val[-1]
                else: val = val[0]
                table[key] += [val]

        # Extra inferred data which is handy for summarising.
        if self.ensemble == Ensemble.NPT:
            headings += ['Density']
            table['Density'] = self.natoms / numpy.array(table['Volume'])
        headings.insert(1, 'Time')
        table['Time'] = numpy.array(table['Step']) * self.dt

        # Return in an easily printable pandas format.
        return pandas.DataFrame(data=table, columns=headings)

    @property
    def snapshot(self):
        """Current snapshot of simulation state.
        
        Returns:
            snapshot (mdma.atom.AtomSnapshot): The current state.
        """
        return atom.AtomSnapshot(self.coordinates, self.box, self.species)

    @property
    def natoms(self):
        """Total number of atoms.
        
        Returns:
            natoms (int): The number of atoms.
        """
        return self.system.natoms

    @property
    def coordinates(self):
        """Atom coordinates.
        
        Returns:
            x (numpy.ndarray): n by d array giving coordinates of the atoms.
        """
        x = self.lmp.gather_atoms('x', 1, 3)
        return numpy.array(x).reshape(-1, 3)

    @coordinates.setter
    def coordinates(self, x):
        """Atom coordinates.
        
        Args:
            x (numpy.ndarray): n by d array giving coordinates of the atoms.
        """
        self.lmp.scatter_atoms('x', 1, 3, x.ctypes)

    @property
    def species(self):
        """Species of each atom.
        
        Returns:
            x (numpy.ndarray): n-dimensional array giving species of the atoms.
        """
        species = self.lmp.gather_atoms('type', 0, 1)
        return numpy.array(species)

    @species.setter
    def species(self, species):
        """Species of each atom.
        
        Args:
            x (numpy.ndarray): n-dimensional array giving species of the atoms.
        """
        species = species.astype(ctypes.c_int)
        self.lmp.scatter_atoms('type', 0, 1, species.ctypes)

    @property
    def box(self):
        """The simulation box.
        
        Returns:
            box (numpy.ndarray): a 2 by d array giving the boundaries in each dimension.
        """
        return numpy.array([[self.system.xlo, self.system.xhi],
                            [self.system.ylo, self.system.yhi],
                            [self.system.zlo, self.system.zhi]])

    @property
    def box_dimensions(self):
        """The simulation box's dimensions.
        
        Returns:
            box (numpy.ndarray): a d-dimensional array giving the width of each dimension.
        """
        return numpy.array([self.system.xhi - self.system.xlo,
                            self.system.yhi - self.system.ylo,
                            self.system.zhi - self.system.zlo])
    @property
    def volume(self):
        """Total system volume.
        
        Returns:
            volume (scalar): total system volume.
        """
        return numpy.prod(self.box_dimensions)

    @property
    def density(self):
        """The number density.
        
        Returns:
            density (scalar): atoms per unit volume.
        """
        return self.natoms / self.volume

    @property
    def composition(self):
        """Get composition of system, i.e. the number of atoms of each species as
        a fraction of the total number of particles.
        
        Returns:
            composition (numpy.ndarray):
                An m dimensional array where m in the total number of unique species,
                giving the compositional make-up of each species.
                The array is ordered by species number.
        """
        species, proportions = numpy.unique(self.species, return_counts=True)
        order = numpy.argsort(species)
        return proportions[order] / self.natoms

    @composition.setter
    def composition(self, composition):
        """Set composition of system, i.e. the number of atoms of each species as
        a fraction of the total number of particles.

        Particles are randomly assigned species to achieve the desired composition.

        Args:
            composition (numpy.ndarray):
                An m dimensional array where m in the total number of unique species,
                giving the compositional make-up of each species.
                The array should be ordered by species number.
        """
        assert sum(composition) == 1

        if len(composition) > 1:
            species = (self.natoms*ctypes.c_int)()
            species = numpy.ones(self.natoms, dtype=ctypes.c_int)

            # Maintain a list of atoms which have not been assigned a species so we can ensure the exact requested composition is achieved.
            remaining_atoms = numpy.arange(self.natoms, dtype=ctypes.c_int)

            for i,proportion in enumerate(composition[1:]):
                # Atoms to convert to the new type chosen randomly without replacement.
                to_convert = numpy.random.choice(len(remaining_atoms),
                                                 int(proportion*self.natoms),
                                                 replace=False)

                species[to_convert] = i+2
                #for j in to_convert: species[remaining_atoms[j]] = i+2
                # Update the list of atoms which have not been assigned a type yet.
                remaining_atoms = numpy.delete(remaining_atoms, to_convert)

            species[remaining_atoms] = 1
            #for i in remaining_atoms: species[i] = 1
            self.species = species
