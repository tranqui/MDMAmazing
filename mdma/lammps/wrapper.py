#!/usr/bin/env python3

import numpy
import os, json, ctypes
from glob import glob
from natsort import natsorted

import lammps

def random_seed():
    return numpy.random.randint(numpy.iinfo(ctypes.c_int).max)

atom_columns = ['id', 'type', 'x', 'y', 'z', 'fx', 'fy', 'fz'] #, 'vx', 'vy', 'vz']

class LammpsExecutable(lammps.PyLammps):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, *kwargs)
        self.exec_id = 0

    def save(self, path, overwrite=True):
        os.makedirs(path, exist_ok=overwrite)
        os.makedirs('%s/dumps' % path, exist_ok=overwrite)

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
        latest_restart = natsorted(glob('%s/dumps/*.restart' % path))[-1]
        self.read_restart(latest_restart)

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
            if active: self.command(latest_command)

    def read_script(self, path):
        with open(path) as f:
            return f.read().splitlines()

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
    def box_size(self):
        return numpy.array([self.system.xhi - self.system.xlo,
                            self.system.yhi - self.system.ylo,
                            self.system.zhi - self.system.zlo])
    @property
    def volume(self):
        return numpy.prod(self.box_size)

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
