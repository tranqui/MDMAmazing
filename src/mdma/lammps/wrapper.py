#!/usr/bin/env python3

import numpy, pandas
import os, ast, json, string, ctypes
from glob import glob
from natsort import natsorted

import lammps
from mdma import atom

def random_seed():
    return numpy.random.randint(numpy.iinfo(ctypes.c_int).max)

atom_columns = ['id', 'type', 'x', 'y', 'z', 'fx', 'fy', 'fz'] #, 'vx', 'vy', 'vz']

class Ensemble:
    NVT = 1
    NPT = 2

class LammpsExecutable(lammps.PyLammps):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.exec_id = 0
        self.initialised = False

    def initialise_system(self, system, initial_snapshot=None):
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
                raise RuntimeError('must give initial snapshot for NPT run!')

        else:
            raise RuntimeError('unknown ensemble encountered (ensemble=%d)!' % ensemble)
        self.ensemble = ensemble

        # Run for zero steps so lammps will compute thermodynamic properties for initial state.
        self.run(0)

    def initial_box(self, system, initial_snapshot, seed_density=1e-3):
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

            # Density was not specified by system: ignore (we must be in NPT ensemble)
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

        latest_log = previous_logs[-1]
        self.exec_id = latest_log.split('/')[-1]
        self.exec_id = self.exec_id.split('.history')[0]
        self.exec_id = 1 + int(self.exec_id)

    def restore_fixes(self, history):
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
        with open(path) as f:
            return f.read().splitlines()

    @property
    def thermodynamics(self):
        headings = ['Step', 'Temp', 'Press', 'TotEng']
        if self.ensemble == Ensemble.NPT: headings += ['Volume']

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

        if self.ensemble == Ensemble.NPT:
            headings += ['Density']
            table['Density'] = self.natoms / numpy.array(table['Volume'])

        return pandas.DataFrame(data=table, columns=headings)

    @property
    def snapshot(self):
        return atom.AtomSnapshot(self.coordinates, self.box, self.species)

    @property
    def natoms(self):
        return self.system.natoms

    @property
    def coordinates(self):
        x = self.lmp.gather_atoms('x', 1, 3)
        return numpy.array(x).reshape(-1, 3)

    @coordinates.setter
    def coordinates(self, x):
        self.lmp.scatter_atoms('x', 1, 3, x.ctypes)

    @property
    def species(self):
        species = self.lmp.gather_atoms('type', 0, 1)
        return numpy.array(species)

    @species.setter
    def species(self, species):
        species = species.astype(ctypes.c_int)
        self.lmp.scatter_atoms('type', 0, 1, species.ctypes)

    @property
    def box(self):
        return numpy.array([[self.system.xlo, self.system.xhi],
                            [self.system.ylo, self.system.yhi],
                            [self.system.zlo, self.system.zhi]])

    @property
    def box_dimensions(self):
        return numpy.array([self.system.xhi - self.system.xlo,
                            self.system.yhi - self.system.ylo,
                            self.system.zhi - self.system.zlo])
    @property
    def volume(self):
        return numpy.prod(self.box_dimensions)

    @property
    def density(self):
        return self.natoms / self.volume

    @property
    def composition(self):
        species, proportions = numpy.unique(self.species, return_counts=True)
        order = numpy.argsort(species)
        return proportions[order] / self.natoms

    @composition.setter
    def composition(self, composition):
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
